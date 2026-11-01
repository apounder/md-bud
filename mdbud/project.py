"""Explicit campaign-to-analysis handoff, without filename discovery heuristics."""
from __future__ import annotations

import copy
import hashlib
import json
import math
import os
import re
from decimal import Decimal, ROUND_FLOOR
from pathlib import Path

from . import campaign
from mdwb.config import normalize_config


def inside(root, *parts):
    path = root.joinpath(*parts).resolve()
    if not path.is_relative_to(root):
        raise ValueError("Campaign paths must stay inside the campaign directory.")
    return path


def campaign_config(project, *, template=None, level=None, discard_ns=None, stride=None,
                    output=None, replica=None, check_files=True):
    root = Path(project).expanduser().resolve()
    manifest_path = root / "md-bud.json"
    manifest_bytes = manifest_path.read_bytes()
    manifest = campaign.load(root)
    settings = manifest["config"]
    # Inputs are deliberately identical to what generation copied and submission checks.
    for name, expected in manifest["inputs"].items():
        if not isinstance(name, str) or Path(name).name != name:
            raise ValueError("Invalid campaign input filename.")
        path = inside(root, "inputs", name)
        with path.open("rb") as stream:
            actual = hashlib.file_digest(stream, "sha256").hexdigest()
        if actual != expected:
            raise ValueError(f"Campaign input changed: {name}. Use the original inputs or generate a new campaign.")
    dt = settings.get("dt")
    ntwx = settings.get("ntwx")
    steps = settings.get("steps")
    if (type(dt) not in (int, float) or not math.isfinite(dt) or dt <= 0 or
            type(ntwx) is not int or ntwx <= 0 or type(steps) is not int or steps < ntwx):
        raise ValueError("Invalid campaign timestep, output interval or segment length.")
    if steps % ntwx:
        raise ValueError("Segment length is not divisible by ntwx. A uniform time axis would be incorrect; reprocess trajectories with explicit time metadata first.")
    replicas = []
    for rep in campaign.selected(manifest, replica):
        names = [stage["name"] for stage in rep["stages"] if re.fullmatch(r"md-\d+", stage["name"])]
        if not names:
            raise ValueError("This campaign contains no production trajectories (prep-only workflow).")
        indices = [int(name.split("-")[1]) for name in names]
        if indices != list(range(1, settings["segments"] + 1)):
            raise ValueError("Production segments must be consecutive, ordered and complete in the manifest.")
        paths = [inside(root, rep["name"], name + ".nc") for name in names]
        if check_files:
            for name, path in zip(names, paths):
                stage = next(s for s in rep["stages"] if s["name"] == name)
                if not campaign.complete(root / rep["name"], stage):
                    raise ValueError(f"Production stage {rep['name']}/{name} is not complete. Wait for completion or use --plan/--save-config to prepare analysis.")
                if path.stat().st_size == 0:
                    raise ValueError(f"Empty production trajectory: {path}")
        replicas.append({"name": rep["name"], "trajectories": [str(p) for p in paths]})
    config = copy.deepcopy(template or {})
    if not isinstance(config, dict):
        raise ValueError("Analysis template must be a JSON object.")
    config.update(topology=str(inside(root, "inputs", "topology.prmtop")), replicas=replicas,
                  output=str(Path(output).expanduser().resolve() if output else root / "analysis_results"))
    if level is not None or not template:
        config["analysis"] = {"level": level or "basic"}
    frames = config.setdefault("frames", {})
    # This is the source saved-frame interval, before analysis stride.
    frames["dt_ps"] = float(Decimal(str(dt)) * ntwx)
    if discard_ns is not None:
        cutoff = Decimal(str(discard_ns)) * 1000 / Decimal(str(frames["dt_ps"]))
        frames["start"] = int(cutoff.to_integral_value(rounding=ROUND_FLOOR)) + 1
    if stride is not None:
        frames["stride"] = stride
    if frames.get("start", 1) > (steps // ntwx) * settings["segments"]:
        raise ValueError("Discard/range selection removes every expected production frame.")
    config.setdefault("study", {}).setdefault("title", f"MD-bud: {root.name}")
    config["notes"] = {"campaign_manifest_sha256": hashlib.sha256(manifest_bytes).hexdigest(),
                       "campaign_expected_frames_per_segment": steps // ntwx,
                       "source": "MD-bud campaign; ordered production segments only",
                       "replica_caution": "Shared starting coordinates/velocities do not establish statistical independence."}
    return normalize_config(config, root, check_files=check_files)


def portable_config(config, destination):
    """Relative paths allow exporting and moving a whole campaign together."""
    data = copy.deepcopy(config)
    parent = Path(destination).resolve().parent
    def relative(value):
        try:
            return os.path.relpath(value, parent)
        except ValueError:  # Windows volumes cannot be expressed relatively.
            return str(value)
    for key in ("topology", "output", "reference"):
        if data.get(key):
            data[key] = relative(data[key])
    for rep in data["replicas"]:
        rep["trajectories"] = [relative(p) for p in rep["trajectories"]]
    if data.get("convergence", {}).get("reference_path"):
        data["convergence"]["reference_path"] = relative(data["convergence"]["reference_path"])
    return data


def write_new_json(path, data):
    path = Path(path).expanduser().resolve()
    path.parent.mkdir(parents=True, exist_ok=True)
    # Refuse existing configs/profiles, including symlink aliases.
    with path.open("x", encoding="utf-8", newline="\n") as stream:
        json.dump(data, stream, indent=2, allow_nan=False)
        stream.write("\n")
    return path
