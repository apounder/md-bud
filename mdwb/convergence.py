"""Observable-specific, descriptive convergence curves; no automatic pass/fail."""
from __future__ import annotations

import itertools
from .figure_names import metric_figure_stem
import csv
import json
from pathlib import Path
import numpy as np
from .convergence_config import validate_references
from .diagnostics import _write, js_distance, sampling_diagnostic

DEFAULTS = {"enabled": True, "checkpoints": 20, "window_frames": 50,
            "reference_mode": "pooled", "tail_ns": 10.0, "reference_path": None,
            "references": {}, "metrics": []}


def split_rhat(arrays):
    """Classical split R-hat, not rank-normalized; equal-length prefix halves."""
    n = min(map(len, arrays)) // 2
    if len(arrays) < 2 or n < 2:
        return None
    chains = np.array([half for a in arrays for half in (a[:n], a[n:2*n])])
    if not np.isfinite(chains).all():
        return None
    variances = np.var(chains, axis=1, ddof=1)
    if np.any(variances == 0):
        return None  # Constant/stuck halves cannot establish mixing.
    within = variances.mean()
    between = n * np.var(chains.mean(axis=1), ddof=1)
    return float(np.sqrt(((n - 1) * within / n + between / n) / within))


def moments(values, circular=False):
    if not circular:
        return float(np.mean(values)), float(np.std(values))
    z = np.mean(np.exp(1j * np.deg2rad(values)))
    r = abs(z)
    if r < 1e-12:
        return None, None
    return float(np.rad2deg(np.angle(z))), float(np.rad2deg(np.sqrt(-2 * np.log(min(1, r)))))


def write_convergence(config, series_list, output):
    opt = {**DEFAULTS, **config.get("convergence", {})}
    if not opt["enabled"]:
        return {"enabled": False}
    folder = Path(output) / "reports" / "convergence"
    folder.mkdir(parents=True, exist_ok=True)
    refs = opt["references"]
    if opt["reference_mode"] == "file":
        refs = json.loads(Path(opt["reference_path"]).read_text())
    validate_references(refs)
    groups, omissions = {}, []
    limits = config.get("diagnostics", {})
    for s in series_list:
        identity = '/'.join((s.analysis, s.section, s.name))
        if s.kind != "timeseries" or (opt["metrics"] and identity not in opt["metrics"] and s.analysis not in opt["metrics"]):
            continue
        groups.setdefault(s.key, []).append(s)
    within, between, populations, reference_rows, histogram_rows = [], [], [], [], []
    used_series = used_pairs = 0
    for key, group in groups.items():
        s = group[0]
        identity = '/'.join(key[:3])
        arrays = [np.asarray(t.values, float) for t in group]
        pair_count = len(group) * (len(group)-1) // 2
        reason = None
        if used_series + len(group) > limits.get("max_series", 100) or sum(map(len, arrays)) > limits.get("max_samples", 1000000) or used_pairs + pair_count > limits.get("max_pairs", 10000):
            reason = "resource_limit"
        elif any(len(a) < 4 or not np.isfinite(a).all() for a in arrays):
            reason = "short_or_nonfinite_no_compaction"
        elif any(len(t.x) != len(t.values) or not np.isfinite(t.x).all() or not np.all(np.diff(t.x) > 0) or not np.allclose(np.diff(t.x), np.diff(t.x)[0], rtol=1e-6, atol=1e-12) for t in group):
            reason = "irregular_sampling"
        elif not all(t.shared_basis for t in group) and len(group) > 1 and s.analysis in {"cluster", "pca"}:
            reason = "cluster_labels_or_pca_basis_not_shared"
        if reason:
            omissions.append([identity, reason]); continue
        used_series += len(group)
        used_pairs += pair_count
        if s.circular:
            arrays = [(a + 180) % 360 - 180 for a in arrays]
        states = np.unique(np.concatenate(arrays)) if s.discrete else None
        if states is None:
            lo, hi = min(a.min() for a in arrays), max(a.max() for a in arrays)
            edges = np.linspace(-180, 180, limits.get("histogram_bins", 32)+1) if s.circular else np.linspace(lo-0.5 if lo == hi else lo, hi+0.5 if lo == hi else hi, limits.get("histogram_bins", 32)+1)
        def probability(a):
            counts = np.array([(a == state).sum() for state in states]) if states is not None else np.histogram(a, edges)[0]
            return counts / counts.sum()
        reference_arrays = arrays
        if opt["reference_mode"] == "tail":
            if s.xunit != "time (ns)" or any(t.x[-1]-t.x[0] < opt["tail_ns"] for t in group):
                omissions.append([identity, "tail_requires_ns_and_full_requested_duration"]); continue
            reference_arrays = [a[t.x > t.x[-1]-opt["tail_ns"]] for t, a in zip(group, arrays)]
        ref_probability = None
        ref_mean = ref_sd = None
        if opt["reference_mode"] in {"pooled", "tail"}:
            ref_probability = np.mean([probability(a) for a in reference_arrays], axis=0)
            if not s.discrete:
                if s.circular:
                    z = np.mean([np.mean(np.exp(1j*np.deg2rad(a))) for a in reference_arrays])
                    if abs(z) >= 1e-12:
                        ref_mean, ref_sd = float(np.rad2deg(np.angle(z))), float(np.rad2deg(np.sqrt(-2*np.log(min(1, abs(z))))))
                else:
                    ref_mean = float(np.mean([a.mean() for a in reference_arrays]))
                    ref_sd = float(np.sqrt(np.mean([np.mean((a-ref_mean)**2) for a in reference_arrays])))
        else:
            ref = refs.get(identity)
            if ref is None or ref["unit"] != s.unit or s.circular or s.discrete:
                omissions.append([identity, "reference_missing_unit_mismatch_or_non_linear"])
            else:
                ref_mean, ref_sd = ref["mean"], ref["sd"]
        reference_rows.append([identity, s.unit, opt["reference_mode"], ref_mean, ref_sd])
        if ref_probability is not None:
            histogram_rows.extend([identity, states[i] if states is not None else edges[i],
                                   None if states is not None else edges[i+1], p]
                                  for i, p in enumerate(ref_probability))
        duration = min(t.x[-1]-t.x[0] for t in group)
        # Common elapsed coordinates: unequal replicas never extend the comparison horizon.
        ends = [t.x[-1]-t.x[0] for t in group]
        checkpoints = np.unique(np.concatenate((
            np.linspace(0, duration, opt["checkpoints"]+1)[1:],
            np.linspace(0, max(ends), opt["checkpoints"]+1)[1:], ends)))
        for elapsed in checkpoints:
            prefixes = [a[t.x-t.x[0] <= elapsed+1e-12] for t, a in zip(group, arrays)]
            if min(map(len, prefixes)) < 4:
                continue
            for t, a in zip(group, prefixes):
                if elapsed > t.x[-1]-t.x[0] + 1e-12:
                    continue
                for scope, v in (("cumulative", a), ("window", a[-opt["window_frames"]:])):
                    mean, sd = (None, None) if s.discrete else moments(v, s.circular)
                    error = None if mean is None or ref_mean is None else mean-ref_mean
                    if error is not None and s.circular:
                        error = (error+180) % 360-180
                    diag = sampling_diagnostic(v, limits.get("max_lag", 10000))[0] if not s.discrete and not s.circular else {}
                    within.append([identity, t.replica, elapsed, s.xunit, scope, len(v), mean, sd, ref_mean, ref_sd, error,
                                   None if sd is None or ref_sd is None else sd-ref_sd,
                                   js_distance(probability(v), ref_probability) if ref_probability is not None else None,
                                   diag.get("ess"), diag.get("sem_stationary")])
                    if s.discrete:
                        populations.extend([identity, t.replica, elapsed, scope, state, p] for state, p in zip(states, probability(v)))
            if elapsed > duration + 1e-12:
                continue
            for scope in ("cumulative", "window"):
                samples = prefixes if scope == "cumulative" else [a[-opt["window_frames"]:] for a in prefixes]
                rhat = split_rhat(samples) if not s.circular and not s.discrete else None
                for i, j in itertools.combinations(range(len(group)), 2):
                    between.append([identity, group[i].replica, group[j].replica, elapsed, s.xunit, scope,
                                    len(samples[i]), len(samples[j]), js_distance(probability(samples[i]), probability(samples[j])), rhat])
    _write(folder / "within_replica.dat", ["metric", "replica", "elapsed", "coordinate_unit", "scope", "n", "mean", "sd", "reference_mean", "reference_sd", "mean_error", "sd_error", "reference_js_distance", "ess_stationary", "sem_stationary"], within)
    _write(folder / "between_replicas.dat", ["metric", "replica_a", "replica_b", "elapsed", "coordinate_unit", "scope", "n_a", "n_b", "js_distance", "all_replica_classical_split_rhat"], between)
    _write(folder / "state_populations.dat", ["metric", "replica", "elapsed", "scope", "state", "fraction"], populations)
    _write(folder / "reference_histograms.dat", ["metric", "bin_low_or_state", "bin_high", "probability"], histogram_rows)
    _write(folder / "references.dat", ["metric", "unit", "mode", "mean", "sd"], reference_rows)
    _write(folder / "omissions.dat", ["metric", "reason"], omissions)
    (folder / "README.md").write_text(
        "Convergence curves are descriptive, with no pass/fail threshold. See the package docs/CONVERGENCE.md.\n"
        "Cumulative prefixes and trailing windows retain original sampling; no automatic equilibration removal.\n"
        "Within-replica curves cover each replica; between-replica curves stop at the shortest elapsed duration.\n"
        "Pooled/tail references weight replicas equally, using population SD of the mixture, not SD of replica means.\n"
        "Tail intervals are (end - tail_ns, end]. Self-derived references are not independent gold standards.\n"
        "Circular means/SD use degree-based circular statistics; categories use populations only.\n"
        "JS distance uses shared fixed bins (saved in reference_histograms.dat), base-2 logs, range 0 to 1.\n"
        "Custom/file mean and SD do not define a distribution; their JS distance is unavailable.\n"
        "Classical split R-hat uses equal-length prefix halves from all replicas, not rank normalization.\n"
        "Constant halves have undefined R-hat; ESS and SEM assume stationarity, not independence of frames.\n"
        "Inspect omissions.dat for exclusions; missing values are not compacted.\n",
        encoding="utf-8")
    return {"enabled": True, "within_rows": len(within), "between_rows": len(between), "omissions": len(omissions)}


def plot_convergence(plt, output, save_figure, files, figure_analyses):
    """Bounded per-metric SVGs, also discoverable in the existing HTML gallery."""
    output = Path(output)
    folder = output / "reports" / "convergence"
    def read(name):
        with (folder / name).open(encoding="utf-8") as handle:
            return list(csv.DictReader(handle, delimiter="\t"))
    within = read("within_replica.dat")
    between = read("between_replicas.dat")
    units = {r["metric"]: r["unit"] for r in read("references.dat")}
    for metric in dict.fromkeys(r["metric"] for r in within):
        rows = [r for r in within if r["metric"] == metric]
        pairs = [r for r in between if r["metric"] == metric and r["scope"] == "cumulative"]
        fig, axes = plt.subplots(3, 1, figsize=(8, 9), sharex=True)
        for replica in dict.fromkeys(r["replica"] for r in rows):
            for scope, style in (("cumulative", "-"), ("window", "--")):
                curve = [r for r in rows if r["replica"] == replica and r["scope"] == scope]
                for ax, field in zip(axes[:2], ("mean_error", "reference_js_distance")):
                    valid = [r for r in curve if r[field] != ""]
                    if valid:
                        ax.plot([float(r["elapsed"]) for r in valid], [float(r[field]) for r in valid],
                                style, label=f"{replica} {scope}")
        for pair in dict.fromkeys((r["replica_a"], r["replica_b"]) for r in pairs):
            curve = [r for r in pairs if (r["replica_a"], r["replica_b"]) == pair]
            axes[2].plot([float(r["elapsed"]) for r in curve], [float(r["js_distance"]) for r in curve], label=" / ".join(pair))
        axes[0].set(title=metric, ylabel=f"Mean − reference ({units.get(metric, '')})")
        axes[1].set(ylabel="Reference JS distance", ylim=(0, 1))
        axes[2].set(ylabel="Between-replica JS distance", ylim=(0, 1), xlabel="Elapsed " + rows[0]["coordinate_unit"])
        for ax in axes:
            if ax.lines:
                ax.legend(fontsize=7)
            else:
                ax.text(.5, .5, "Unavailable for this metric/reference", ha="center", transform=ax.transAxes)
        token = metric_figure_stem(*metric.rsplit("/", 2))
        path = output / "figures" / "convergence" / f"{token}__convergence-over-time.svg"
        save_figure(plt, fig, path, files, output)
        figure_analyses[path.relative_to(output).as_posix()] = metric.split('/')[0]
