"""Streaming descriptive AMBER production-log checks; no equilibrium verdict."""
from __future__ import annotations

import math
from pathlib import Path
import re

from . import campaign
from .project import inside

FIELDS = {"TEMP(K)": "temperature_K", "PRESS": "pressure_bar", "Etot": "total_energy_kcal_mol",
          "EPtot": "potential_energy_kcal_mol", "VOLUME": "volume_A3", "Density": "density_g_cm3"}
VALUE = re.compile(r"(TEMP\(K\)|PRESS|Etot|EPtot|VOLUME|Density)\s*=\s*(\S+)")
FAILURE = re.compile(r"(?:^\s*(?:ERROR|FATAL)\b|SHAKE cannot|vlimit exceeded|SANDER BOMB|PMEMD BOMB)", re.I)


def inspect_log(path):
    metrics, errors = {}, []
    samples = 0
    in_samples = False
    with Path(path).open(encoding="utf-8", errors="replace") as stream:
        for line_number, line in enumerate(stream, 1):
            if FAILURE.search(line) and len(errors) < 20:
                errors.append(f"line {line_number}: {line.strip()[:240]}")
            if re.search(r"A\s+V\s+E\s+R\s+A\s+G\s+E\s+S|R\s+M\s+S\s+F\s+L\s+U\s+C", line):
                in_samples = False
            if re.search(r"\bNSTEP\s*=\s*\d+", line):
                in_samples = True
                samples += 1
            if not in_samples:
                continue
            for field, raw in VALUE.findall(line):
                try:
                    value = float(raw.replace("D", "E").replace("d", "e"))
                    if not math.isfinite(value):
                        raise ValueError()
                except ValueError:
                    if len(errors) < 20:
                        errors.append(f"line {line_number}: invalid {field}={raw}")
                    continue
                key = FIELDS[field]
                item = metrics.setdefault(key, {"count": 0, "mean": 0.0, "min": value, "max": value, "last": value})
                item["count"] += 1
                item["mean"] += (value - item["mean"]) / item["count"]
                item["min"], item["max"], item["last"] = min(item["min"], value), max(item["max"], value), value
    return {"samples": samples, "metrics": metrics, "errors": errors}


def inspect_campaign(project, replica=None):
    root = Path(project).expanduser().resolve()
    manifest = campaign.load(root)
    rows, issues = [], []
    target = manifest["config"]["temperature"]
    for rep in campaign.selected(manifest, replica):
        for stage in rep["stages"]:
            if not re.fullmatch(r"md-\d+", stage["name"]):
                continue
            path = inside(root, rep["name"], stage["name"] + ".out")
            row = {"replica": rep["name"], "stage": stage["name"], "file": str(path)}
            if not path.is_file():
                row["status"] = "missing"
                issues.append(f"Missing log: {rep['name']}/{stage['name']}")
            else:
                row.update(inspect_log(path))
                row["status"] = "errors" if row["errors"] else "readable"
                issues.extend(f"{rep['name']}/{stage['name']}: {e}" for e in row["errors"])
                if row["samples"] == 0 or "temperature_K" not in row["metrics"]:
                    issues.append(f"No temperature samples: {rep['name']}/{stage['name']}")
                temperature = row["metrics"].get("temperature_K")
                if temperature and abs(temperature["mean"] - target) > 20:
                    issues.append(f"Mean temperature differs from target by >20 K: {rep['name']}/{stage['name']}")
            rows.append(row)
    if not rows:
        issues.append("No production stages in this campaign.")
    return {"target_temperature_K": target, "stages": rows, "issues": issues,
            "note": "Descriptive production-log checks only. The 20 K temperature flag is a heuristic; no energy-drift, equilibration, convergence or ensemble-validity test is performed. Missing quantities are not inferred."}
