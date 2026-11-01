"""Dependency-free validation for convergence reference configuration."""
import json
import math
from pathlib import Path

def validate(options, base, dt_ps, check_files):
    if type(options["enabled"]) is not bool:
        raise ValueError("convergence.enabled must be boolean")
    for key in ("checkpoints", "window_frames"):
        if type(options[key]) is not int or not 2 <= options[key] <= 10000:
            raise ValueError(f"convergence.{key} must be an integer from 2 to 10000")
    if options["reference_mode"] not in {"pooled", "tail", "custom", "file"}:
        raise ValueError("convergence.reference_mode: pooled, tail, custom, file")
    if type(options["tail_ns"]) not in (int, float) or not math.isfinite(options["tail_ns"]) or options["tail_ns"] <= 0:
        raise ValueError("convergence.tail_ns must be finite and positive")
    if options["enabled"] and options["reference_mode"] == "tail" and not dt_ps:
        raise ValueError("Last-ns references require frames.dt_ps (saved-frame interval)")
    if not isinstance(options["metrics"], list) or not all(isinstance(s, str) and s for s in options["metrics"]):
        raise ValueError("convergence.metrics must list analysis names or analysis/section/series IDs")
    if options["reference_path"] is not None:
        if not isinstance(options["reference_path"], str) or not options["reference_path"].strip():
            raise ValueError("convergence.reference_path must be a path")
        options["reference_path"] = str((Path(base) / Path(options["reference_path"]).expanduser()).resolve())
    if options["reference_mode"] == "file":
        if not options["reference_path"]:
            raise ValueError("File reference requires convergence.reference_path")
        if check_files:
            validate_references(json.loads(Path(options["reference_path"]).read_text()))
    validate_references(options["references"])
    if options["reference_mode"] == "custom" and not options["references"]:
        raise ValueError("Custom reference requires metric-specific convergence.references")


def validate_references(refs):
    if not isinstance(refs, dict):
        raise ValueError("References must map analysis/section/series IDs to mean, sd and unit")
    for key, ref in refs.items():
        if not isinstance(key, str) or len(key.rsplit('/', 2)) != 3 or not all(key.split('/')) or not isinstance(ref, dict) or set(ref) != {"mean", "sd", "unit"}:
            raise ValueError("Each reference needs analysis/section/series: {mean, sd, unit}")
        if not isinstance(ref["unit"], str) or not ref["unit"] or any(type(ref[k]) not in (int, float) or not math.isfinite(ref[k]) for k in ("mean", "sd")) or ref["sd"] < 0:
            raise ValueError("Reference mean/SD must be finite, SD nonnegative, with explicit units")


