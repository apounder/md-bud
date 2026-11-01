"""Self-contained browser previews built from numerical results, without Matplotlib."""
from __future__ import annotations

import base64
import hashlib
import importlib.resources
import json
from pathlib import Path

import numpy as np

from .report_pages import GUIDANCE, _read
from .figure_style import PALETTES
from .labels import analysis_label, observable_label, readable, unit_label
from .distributions import shared_edges, histogram


def write_preview_data(config, series_list, output):
    options = config.get("reports", {})
    maximum = options.get("preview_series", 500)
    points = options.get("preview_points", 600)
    groups = {}
    for series in series_list:
        groups.setdefault(series.key, []).append(series)
    records = []
    omitted = 0
    for key, group in groups.items():
        # Keep all replicas of an observable together, or omit the whole group.
        if len(records) + len(group) > maximum:
            omitted += len(group)
            continue
        finite_arrays = [s.values[np.isfinite(s.values)] for s in group]
        categorical = group[0].discrete
        circular = group[0].circular
        if categorical:
            states = sorted({float(v) for a in finite_arrays for v in np.unique(a)})
            if len(states) > 100:
                omitted += len(group)
                continue
            bins = None
        else:
            bins = shared_edges(group)
        for series, finite in zip(group, finite_arrays):
            values = (series.values + 180) % 360 - 180 if circular else series.values
            samples = []
            # Ordered min/max envelopes preserve extrema; any missing bucket creates a gap.
            for indices in np.array_split(np.arange(len(values)), max(1, min(len(values), points // 4))):
                if not len(indices):
                    continue
                valid = np.isfinite(values[indices]) & np.isfinite(series.x[indices])
                if not np.all(valid):
                    samples.append([None, None])
                    continue
                positions = sorted({int(indices[0]), int(indices[-1]), int(indices[np.argmin(values[indices])]), int(indices[np.argmax(values[indices])])})
                samples.extend([[float(series.x[i]), float(values[i])] for i in positions])
            if categorical:
                counts = np.array([np.count_nonzero(finite == state) for state in states])
                labels = [f"State {state:g}" for state in states]
            else:
                counts = histogram(series, bins)["count"]
                labels = [f"{a:.4g}–{b:.4g}" for a, b in zip(bins[:-1], bins[1:])]
            fractions = (counts / len(finite)).tolist() if len(finite) else [0.0] * len(counts)
            records.append({"id": hashlib.sha256("|".join(key).encode()).hexdigest()[:20],
                            "analysis": series.analysis, "section": series.section, "name": series.name,
                            "replica": series.replica, "kind": series.kind, "unit": series.unit,
                            "analysis_label": analysis_label(series.analysis),
                            "observable_label": observable_label(series), "section_label": readable(series.section),
                            "unit_label": unit_label(series.unit),
                            "xunit": series.xunit, "categorical": categorical, "circular": circular,
                            "points": samples, "total": len(values), "finite": len(finite),
                            "distribution": fractions, "labels": labels,
                            "density": None if categorical else (np.asarray(fractions) / np.diff(bins)).tolist(),
                            "distribution_applicable": series.kind == "timeseries",
                            "edges": None if categorical else bins.tolist(),
                            "preview_reduced": len(values) > points // 4})
    payload = {"series": records, "omitted_series": omitted, "limits": {"series": maximum, "trace_points": points},
               "note": "Traces use ordered min/max envelopes; buckets with missing values remain gaps. Histograms, state fractions and statistics use all finite samples of each included observable. Independent replicas remain separate."}
    (Path(output) / "reports" / "preview_data.json").write_text(json.dumps(payload, allow_nan=False, separators=(",", ":")), encoding="utf-8")


def write_interactive_report(config, manifest, summary, output):
    output = Path(output)
    from .structure_view import structure_payload
    structure = structure_payload(config, manifest, output)
    vendor = importlib.resources.files("mdwb").joinpath("vendor")
    path = output / "reports" / "preview_data.json"
    preview = json.loads(path.read_text(encoding="utf-8")) if path.exists() else {"series": [], "omitted_series": 0, "note": "Regenerate this report to add numerical previews."}
    figures, used = [], 0
    for filename in summary.get("figures", []):
        relative = Path(filename)
        if relative.is_absolute() or ".." in relative.parts:
            continue
        source = output / relative
        if not source.is_file():
            continue
        size = source.stat().st_size
        embedded = None
        if used + size <= 20_000_000:
            embedded = "data:image/svg+xml;base64," + base64.b64encode(source.read_bytes()).decode("ascii")
            used += size
        category = ("Probability distributions" if "distributions" in relative.parts else
                    "Sampling diagnostics" if "diagnostics" in relative.parts else
                    "Replica overlays" if "overlays" in relative.parts else
                    "Replicas in sequence" if "concatenated" in relative.parts else
                    "Matrix summaries" if "matrices" in relative.parts else "Individual replicas")
        figures.append({"path": relative.as_posix(), "analysis": summary.get("figure_analyses", {}).get(filename, ""),
                        "title": summary.get("figure_titles", {}).get(filename, readable(relative.stem)),
                        "category": category, "image": embedded})
    payload = {"title": config.get("study", {}).get("title", "MD analysis"), "description": config.get("study", {}).get("description", ""),
               "structure": structure,
               "viewerScript": base64.b64encode(vendor.joinpath("3Dmol-2.5.3-min.js").read_bytes()).decode("ascii") if structure.get("pdb") else None,
               "viewerLicense": vendor.joinpath("3Dmol-LICENSE.txt").read_text(encoding="utf-8") if structure.get("pdb") else "",
               "palette": list(PALETTES[config.get("plots", {}).get("palette", "lagoon")]),
               "status": "synthetic" if manifest.get("synthetic") else manifest.get("status", "unknown"), "summary": summary, "preview": preview,
               "replicas": [{"name": r["name"], "status": r.get("status", "unknown"), "frames": r.get("frame_count")} for r in manifest.get("replicas", [])],
               "guidance": GUIDANCE, "figures": figures,
               "analysisLabels": {a: analysis_label(a) for a in {r.get("analysis", "") for r in preview["series"]} | set(summary.get("figure_analyses", {}).values()) | {r.get("analysis", "") for r in _read(output, "artifact_index.dat")}},
               "files": [{"path": name, "label": label} for name, label in (
                   ("reports/all_replicas.dat", "Combined data"), ("reports/replicate_statistics.dat", "Replica summary"),
                   ("reports/probability_distributions.dat", "Probability distributions: exact bins and replica counts"),
                   ("reports/artifact_index.dat", "Artifact index"), ("METHODS.md", "Methods"),
                   ("manifest.json", "Manifest"), ("synthetic_manifest.json", "Synthetic manifest")) if (output / name).is_file()],
               "tables": {name: _read(output, filename) for name, filename in {
                   "individual": "descriptive_statistics.dat", "replicate": "replicate_statistics.dat",
                   "circular": "circular_statistics.dat", "circularReplicate": "circular_replicate_statistics.dat",
                   "states": "categorical_fractions.dat", "sampling": "diagnostics/sampling.dat",
                   "findings": "review_findings.dat", "artifacts": "artifact_index.dat"}.items()}}
    if not summary.get("diagnostics", {}).get("enabled"):
        payload["tables"]["sampling"] = []
    template = importlib.resources.files("mdwb").joinpath("report_template.html").read_text(encoding="utf-8")
    encoded = json.dumps(payload, allow_nan=False, separators=(",", ":")).replace("<", "\\u003c").replace(">", "\\u003e").replace("&", "\\u0026")
    (output / "index.html").write_text(template.replace("__REPORT_DATA__", encoded), encoding="utf-8")
