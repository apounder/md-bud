"""Exploratory sampling diagnostics and histogram-based replica comparisons.

This is independent NumPy code, not a PyMBAR/PENSA implementation or wrapper.
No equilibration cutoff is chosen and no frames are automatically removed.
"""
from __future__ import annotations

import csv
import hashlib
import itertools
import math
from pathlib import Path

import numpy as np


def sampling_diagnostic(values, max_lag=10000):
    """Biased FFT ACF with Geyer initial-positive, monotone paired sequence.

    Pairing starts at lags 0/1, then 2/3. g=-1+2*sum(pair sums),
    constrained to [1,N]; negative-correlation gains are not claimed. Stationary
    scalar sampling and adequate mixing are assumptions, not conclusions.
    """
    values = np.asarray(values, dtype=float)
    n = len(values)
    result = {"status": "too_short", "n": n, "g": None, "ess": None,
              "sem_stationary": None, "suggested_block_frames": None,
              "half_mean_difference": None, "lag_limit_reached": False}
    if n < 32:
        return result, np.array([])
    if not np.all(np.isfinite(values)):
        result["status"] = "missing_values_no_compaction"
        return result, np.array([])
    scale = float(np.max(np.abs(values)))
    if scale == 0:
        result["status"] = "constant"
        return result, np.array([])
    centered = values / scale - np.mean(values / scale)
    variance = float(np.dot(centered, centered) / n)
    if variance <= np.finfo(float).eps ** 2:
        result["status"] = "constant_or_numerically_constant"
        return result, np.array([])
    size = 1 << (2 * n - 1).bit_length()
    spectrum = np.fft.rfft(centered, n=size)
    lag = min(max_lag, n // 2)
    acf = np.fft.irfft(spectrum * spectrum.conjugate(), n=size)[:lag + 1]
    acf /= acf[0]
    pairs = []
    terminated = False
    for i in range(0, len(acf) - 1, 2):
        pair = float(acf[i] + acf[i + 1])
        if pair <= 0:
            terminated = True
            break
        pairs.append(min(pair, pairs[-1]) if pairs else pair)
    g = min(float(n), max(1.0, -1.0 + 2.0 * sum(pairs)))
    result.update(status="estimate_requires_stationarity", g=g, ess=n / g,
                  sem_stationary=float(np.std(values / scale, ddof=1) * scale * math.sqrt(g / n)),
                  suggested_block_frames=int(math.ceil(5 * g)),
                  half_mean_difference=float((np.mean(values[n // 2:] / scale) - np.mean(values[:n // 2] / scale)) * scale),
                  lag_limit_reached=not terminated)
    return result, acf


def js_distance(p, q):
    """Base-2 Jensen-Shannon DISTANCE (sqrt divergence), bounded by [0,1]."""
    p, q = np.asarray(p, float), np.asarray(q, float)
    p, q = p / p.sum(), q / q.sum()
    middle = (p + q) / 2
    keep_p, keep_q = p > 0, q > 0
    divergence = 0.5 * (np.sum(p[keep_p] * np.log2(p[keep_p] / middle[keep_p])) +
                        np.sum(q[keep_q] * np.log2(q[keep_q] / middle[keep_q])))
    return float(np.sqrt(np.clip(divergence, 0, 1)))


def _write(path, columns, rows):
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(columns)
        writer.writerows(rows)


def write_diagnostics(config, series_list, output):
    options = config.get("diagnostics", {})
    if not options.get("enabled", True):
        return {"enabled": False}
    directory = Path(output) / "reports" / "diagnostics"
    directory.mkdir(parents=True, exist_ok=True)
    max_series = options.get("max_series", 100)
    max_samples = options.get("max_samples", 1000000)
    max_lag = options.get("max_lag", 10000)
    bins = options.get("histogram_bins", 32)
    rows, curves = [], []
    processed = 0
    for series in series_list:
        if series.kind != "timeseries" or series.discrete or series.circular:
            continue
        prefix = [series.replica, series.analysis, series.section, series.name, series.unit]
        if processed >= max_series or len(series.values) > max_samples:
            rows.append([*prefix, "resource_limit", len(series.values), *([None] * 7)])
            continue
        processed += 1
        delta = np.diff(series.x)
        if len(delta) and (not np.all(np.isfinite(delta)) or delta[0] <= 0 or not np.allclose(delta, delta[0], rtol=1e-6, atol=1e-12)):
            rows.append([*prefix, "irregular_sampling", len(series.values), *([None] * 7)])
            continue
        result, acf = sampling_diagnostic(series.values, max_lag)
        rows.append([*prefix, result["status"], result["n"], result["g"], result["ess"],
                     result["sem_stationary"], result["suggested_block_frames"],
                     result["half_mean_difference"], result["lag_limit_reached"],
                     bool(result["g"] and config.get("stats", {}).get("block_size", 50) < result["suggested_block_frames"])])
        if len(acf):
            token = hashlib.sha256("|".join(prefix).encode()).hexdigest()[:16]
            path = directory / f"acf_{token}.dat"
            _write(path, ["lag_analyzed_frames", "lag_coordinate", "acf"],
                   [(i, i * delta[0], value) for i, value in enumerate(acf)])
            curves.append([*prefix[:4], str(path.relative_to(output)), series.xunit])
    _write(directory / "sampling.dat",
           ["replica", "analysis", "section", "series", "unit", "status", "n", "g_estimate",
            "effective_samples_estimate", "sem_under_stationarity", "suggested_block_frames",
            "second_minus_first_half_mean", "lag_limit_reached", "configured_block_shorter_than_suggestion"], rows)
    _write(directory / "acf_index.dat", ["replica", "analysis", "section", "series", "path", "lag_unit"], curves)
    groups = {}
    for series in series_list:
        if series.kind == "timeseries":
            groups.setdefault(series.key, []).append(series)
    comparisons, distributions, omitted = [], [], []
    for key, group in groups.items():
        if len(group) < 2:
            continue
        if len(comparisons) + len(group) * (len(group) - 1) // 2 > options.get("max_pairs", 10000):
            omitted.append([*key[:3], "pair_limit"])
            continue
        if len(distributions) >= max_series or sum(len(s.values) for s in group) > max_samples:
            omitted.append([*key[:3], "resource_limit"])
            continue
        arrays = [s.values[np.isfinite(s.values)] for s in group]
        if not all(len(a) for a in arrays):
            omitted.append([*key[:3], "empty_replica"])
            continue
        if group[0].circular:
            arrays = [(a + 180) % 360 - 180 for a in arrays]
            edges = np.linspace(-180, 180, bins + 1)
        elif group[0].discrete:
            states = np.unique(np.concatenate(arrays))
            edges = None
        else:
            low, high = min(float(a.min()) for a in arrays), max(float(a.max()) for a in arrays)
            if low == high:
                low, high = low - 0.5, high + 0.5
            edges = np.linspace(low, high, bins + 1)
        counts = [np.array([(a == state).sum() for state in states]) if edges is None else np.histogram(a, bins=edges)[0] for a in arrays]
        probabilities = [c / c.sum() for c in counts]
        token = hashlib.sha256("|".join(key).encode()).hexdigest()[:16]
        hist_path = directory / f"distribution_{token}.dat"
        _write(hist_path, ["bin_low_or_state", "bin_high", *[s.replica for s in group]],
               [(states[i] if edges is None else edges[i], None if edges is None else edges[i + 1],
                 *[p[i] for p in probabilities]) for i in range(len(probabilities[0]))])
        distributions.append(str(hist_path.relative_to(output)))
        for i, j in itertools.combinations(range(len(group)), 2):
            comparisons.append([*key[:3], group[i].replica, group[j].replica, len(arrays[i]), len(arrays[j]),
                                js_distance(probabilities[i], probabilities[j]), len(probabilities[i]),
                                "circular_fixed_bins" if group[0].circular else "categorical" if edges is None else "shared_linear_bins",
                                str(hist_path.relative_to(output))])
    _write(directory / "replica_distribution_distance.dat",
           ["analysis", "section", "series", "replica_a", "replica_b", "n_a", "n_b", "js_distance_base2",
            "bins_or_states", "binning", "distribution_file"], comparisons)
    _write(directory / "distribution_omissions.dat", ["analysis", "section", "series", "reason"], omitted)
    (directory / "README.md").write_text(
        "# Exploratory sampling diagnostics\n\n"
        "ACFs are computed separately for finite, regularly sampled, nonconstant linear observables with at least 32 samples. "
        "Missing values are not removed to manufacture contiguous sampling. Circular and categorical series are excluded from ESS. "
        "FFT autocovariances use a fixed N denominator. Initial-positive monotone pairs start at lags 0/1, 2/3, etc.; "
        "g = max(1, -1 + 2 sum(pairs)), capped at N. ESS=N/g and SEM=sqrt(sample variance*g/N) assume stationarity. "
        "Neither is proof of equilibration or a calibrated confidence interval. Truncation at the lag limit needs review. "
        "The suggested block size ceil(5*g) is a heuristic; no configuration or data is changed.\n\n"
        "Jensen-Shannon DISTANCE uses base-2 logs and is the square root of the divergence (0 to 1). "
        "Each replica histogram is normalized independently; common bins are retained alongside results. "
        "Circular bins cover [-180,180); categories use shared states. Values are descriptive, sensitive to binning and sample size, "
        "and are not significance tests or convergence thresholds. Identical sampled distributions do not establish adequate exploration.\n\n"
        "Methods: https://www.stat.umn.edu/geyer/mcmc/library/mcmc/html/initseq.html\n"
        "Context: https://pymbar.readthedocs.io/en/latest/timeseries.html and https://github.com/drorlab/pensa\n",
        encoding="utf-8")
    return {"enabled": True, "sampling_rows": len(rows), "acf_files": len(curves),
            "distribution_comparisons": len(comparisons), "distribution_omissions": len(omitted)}
