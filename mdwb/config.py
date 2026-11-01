"""Portable JSON configuration and input discovery, with no third-party imports."""
from __future__ import annotations

import copy
import json
import math
import os
import re
from pathlib import Path

from .topology import inspect_topology
from .selections import normalize_mask, normalize_config_masks

TOPOLOGY_SUFFIXES = {".prmtop", ".parm7", ".top", ".psf", ".pdb", ".mol2", ".gro"}
TRAJECTORY_SUFFIXES = {".nc", ".netcdf", ".mdcrd", ".crd", ".dcd", ".xtc", ".trr", ".binpos", ".h5", ".lammpstrj"}
DEFAULTS = {
    "version": 1, "cpptraj": "cpptraj", "output": "analysis_results",
    "frames": {"start": 1, "stop": -1, "stride": 1, "dt_ps": None},
    "imaging": {"mode": "auto", "anchor": None},
    "reference": None, "fit_mask": "representative",
    "sections": [], "interactions": [], "monitors": [],
    "solvation": [],
    "stacking": {"rings": [], "targets": [], "distance_cutoff": 5.0, "angle_cutoff": 30.0,
                 "center": "geometry", "exclude_same_residue": True, "max_pairs": 2000},
    "analysis": {"level": "standard", "enabled": None, "additional": []},
    "progress": {"enabled": True},
    "nucleic": {"resrange": None, "resmap": {}},
    "advanced": {"matrix_mask": "representative", "pca_modes": 3, "cluster_count": 5,
                 "cluster_sieve": 10, "cluster_metric": "rms", "cluster_features": [],
                 "cluster_algorithm": "kmeans", "pairwise_rmsd": False, "max_frames": 20000,
                 "max_matrix_atoms": 2000, "max_memory_gb": 4.0},
    "plots": {"dpi": 600, "alpha": 0.85, "enabled": True, "font": "Arial",
              "font_size": 11, "width_mm": 180, "palette": "lagoon", "svg_text": "editable"},
    "stats": {"block_size": 50, "metric_correlations": False},
    "reports": {"max_matrix_cells": 1000000, "preview_series": 500, "preview_points": 600},
    "study": {"title": "MD analysis", "description": ""},
    "structure_view": {"enabled": True, "max_atoms": 12000, "max_bytes": 1500000},
    "convergence": {"enabled": True, "checkpoints": 20, "window_frames": 50,
                    "reference_mode": "pooled", "tail_ns": 10.0, "reference_path": None,
                    "references": {}, "metrics": []},
    "diagnostics": {"enabled": True, "max_series": 100, "max_samples": 1000000,
                    "max_lag": 10000, "histogram_bins": 32, "max_pairs": 10000},
}


def natural_key(value):
    return [(0, int(part)) if part.isdigit() else (1, part.lower()) for part in re.split(r"(\d+)", str(value))]


def discover(root):
    root = Path(root).expanduser().resolve()
    found = {"topologies": [], "trajectories": []}
    for path in root.rglob("*"):
        if not path.is_file() or any(part in {".git", ".venv", "__pycache__"} for part in path.relative_to(root).parts):
            continue
        suffix = path.suffix.lower()
        if suffix == ".gz":
            suffix = Path(path.stem).suffix.lower()
        if suffix in TOPOLOGY_SUFFIXES:
            found["topologies"].append(path)
        elif suffix in TRAJECTORY_SUFFIXES:
            found["trajectories"].append(path)
    for paths in found.values():
        paths.sort(key=natural_key)
    return found


def _merge(default, incoming):
    result = copy.deepcopy(default)
    for key, value in incoming.items():
        result[key] = _merge(result[key], value) if isinstance(value, dict) and isinstance(result.get(key), dict) else copy.deepcopy(value)
    return result


def safe_name(value):
    if not isinstance(value, str) or not re.fullmatch(r"[A-Za-z][A-Za-z0-9_-]{0,63}", value):
        raise ValueError(f"Use a name starting with a letter and containing only letters, numbers, _ or -: {value!r}")
    return value


def safe_text(value, label="value"):
    if not isinstance(value, str) or not value.strip() or any(c in value for c in '\r\n\x00";$#'):
        raise ValueError(f"Invalid {label}: blank text or CPPTRAJ control characters (newline, double quote, ;, $, #).")
    return value


def quote(value):
    return '"' + safe_text(str(value), "CPPTRAJ argument") + '"'


def load_config(path, check_files=True):
    path = Path(path).expanduser().resolve()
    return normalize_config(json.loads(path.read_text(encoding="utf-8-sig")), path.parent, check_files)


def normalize_config(data, base=None, check_files=True):
    base = Path(base or ".").resolve()
    if not isinstance(data, dict):
        raise ValueError("Configuration must be a JSON object.")
    unknown = set(data) - (set(DEFAULTS) | {"topology", "replicas", "notes"})
    if unknown:
        raise ValueError(f"Unknown configuration fields: {', '.join(sorted(unknown))}")
    for section, defaults in DEFAULTS.items():
        if isinstance(defaults, dict) and section in data:
            if not isinstance(data[section], dict):
                raise ValueError(f"{section} must be a JSON object.")
            unknown = set(data[section]) - set(defaults)
            if unknown:
                raise ValueError(f"Unknown {section} fields: {', '.join(sorted(unknown))}")
    config = _merge(DEFAULTS, data)
    if type(config["structure_view"]["enabled"]) is not bool:
        raise ValueError("structure_view.enabled must be true or false.")
    for key in ("max_atoms", "max_bytes"):
        if type(config["structure_view"][key]) is not int or config["structure_view"][key] < 1:
            raise ValueError(f"structure_view.{key} must be a positive integer.")
    for key in ("title", "description"):
        if not isinstance(config["study"][key], str):
            raise ValueError(f"study.{key} must be text.")
    if type(config["diagnostics"]["enabled"]) is not bool:
        raise ValueError("diagnostics.enabled must be true or false.")
    for key in ("max_series", "max_samples", "max_lag", "histogram_bins", "max_pairs"):
        value = config["diagnostics"][key]
        if type(value) is not int or value < 2:
            raise ValueError(f"diagnostics.{key} must be an integer >= 2.")
    if type(config["version"]) is not int or config["version"] != 1:
        raise ValueError("Unsupported configuration version; expected 1.")
    safe_text(config["cpptraj"], "cpptraj executable")
    def resolve(value):
        path = Path(safe_text(value, "path")).expanduser()
        return str((base / path).resolve())
    if not config.get("topology"):
        raise ValueError("Specify one matching topology for this comparison.")
    config["topology"] = resolve(config["topology"])
    config["output"] = resolve(config["output"])
    if config["reference"]:
        config["reference"] = resolve(config["reference"])
    replicas = config.get("replicas", [])
    if not isinstance(replicas, list) or not replicas:
        raise ValueError("At least one replica with an ordered trajectories list is required.")
    names, paths, path_ids = set(), set(), set()
    for rep in replicas:
        if not isinstance(rep, dict):
            raise ValueError("Each replica must be an object with name and trajectories.")
        if set(rep) - {"name", "trajectories"}:
            raise ValueError("Replica fields are name and trajectories. Use one common topology and frame interval per comparison.")
        name = safe_name(rep.get("name"))
        if name.casefold() in names:
            raise ValueError(f"Duplicate replica name: {name}")
        names.add(name.casefold())
        if not isinstance(rep.get("trajectories"), list) or not rep["trajectories"]:
            raise ValueError(f"Replica {name} has no trajectories.")
        rep["trajectories"] = [resolve(p) for p in rep["trajectories"]]
        for p in rep["trajectories"]:
            # Windows case aliases and resolved relative/symlink paths must not
            # manufacture independent samples by loading the same file twice.
            path_id = os.path.normcase(p)
            if path_id in path_ids:
                raise ValueError(f"Trajectory used more than once (would duplicate samples): {p}")
            path_ids.add(path_id)
            paths.add(p)
    if config["imaging"]["mode"] not in {"auto", "on", "off"}:
        raise ValueError("imaging.mode must be auto, on, or off.")
    if config["analysis"]["level"] not in {"basic", "standard", "advanced"}:
        raise ValueError("analysis.level must be basic, standard, or advanced.")
    from .analysis import enabled_analyses
    if not isinstance(config["analysis"]["additional"], list) or not all(isinstance(x, str) for x in config["analysis"]["additional"]):
        raise ValueError("analysis.additional must be a list of analysis names.")
    if config["analysis"]["enabled"] is not None and config["analysis"]["additional"]:
        raise ValueError("Use either analysis.enabled (exact list) or a preset with analysis.additional.")
    enabled = enabled_analyses(config)
    from .stacking import validate_stacking
    validate_stacking(config["stacking"], "stacking" in enabled)
    if not isinstance(config["solvation"], list):
        raise ValueError("solvation must be a list of named rdf or watershell specifications.")
    solvation_names = set()
    for item in config["solvation"]:
        if not isinstance(item, dict) or set(item) - {"name", "type", "mask1", "mask2", "spacing", "maximum", "lower", "upper"}:
            raise ValueError("Solvation fields: name, type, mask1, mask2, spacing, maximum, lower, upper.")
        name = safe_name(item.get("name"))
        if name.casefold() in solvation_names:
            raise ValueError("Solvation names must be unique.")
        solvation_names.add(name.casefold())
        kind = item.get("type")
        if kind not in {"rdf", "watershell"} or kind not in enabled:
            raise ValueError("Enable each solvation type (rdf or watershell) explicitly in analysis.additional.")
        for key in ("mask1", "mask2"):
            safe_text(item.get(key), "solvation mask")
        keys = ("spacing", "maximum") if kind == "rdf" else ("lower", "upper")
        if any(k in item for k in ({"spacing", "maximum", "lower", "upper"} - set(keys))):
            raise ValueError("RDF uses spacing/maximum; watershell uses lower/upper.")
        for key in keys:
            value = item.get(key)
            if type(value) not in (float, int) or not math.isfinite(value) or value <= 0:
                raise ValueError(f"solvation.{key} must be positive and finite.")
        if item[keys[0]] >= item[keys[1]]:
            raise ValueError("Solvation spacing/lower cutoff must be below maximum/upper cutoff.")
        if kind == "rdf" and item["maximum"] / item["spacing"] > 10000:
            raise ValueError("RDF requests over 10,000 bins; increase spacing or reduce maximum.")
    for kind in ("rdf", "watershell"):
        if kind in enabled and not any(s["type"] == kind for s in config["solvation"]):
            raise ValueError(f"Configure at least one {kind} solvation selection.")
    plots = config["plots"]
    if not isinstance(plots["font"], str) or not plots["font"].strip():
        raise ValueError("plots.font must be a nonempty font family name.")
    for key, low, high in (("font_size", 7, 24), ("width_mm", 80, 300)):
        if type(plots[key]) not in (int, float) or not math.isfinite(plots[key]) or not low <= plots[key] <= high:
            raise ValueError(f"plots.{key} must be between {low} and {high}.")
    if plots["palette"] not in {"lagoon", "mineral", "colorblind"} or plots["svg_text"] not in {"editable", "paths"}:
        raise ValueError("Choose palette lagoon/mineral/colorblind and svg_text editable/paths.")
    if type(config["progress"]["enabled"]) is not bool:
        raise ValueError("progress.enabled must be true or false.")
    if not isinstance(config["monitors"], list):
        raise ValueError("monitors must be a list.")
    monitor_names = set()
    for monitor in config["monitors"]:
        if not isinstance(monitor, dict) or set(monitor) - {"name", "type", "masks", "center", "image", "range360"}:
            raise ValueError("Monitor fields: name, type, masks, center, image (distance), range360 (dihedral).")
        name = safe_name(monitor.get("name"))
        if name.casefold() in monitor_names:
            raise ValueError(f"Duplicate monitor name: {name}")
        monitor_names.add(name.casefold())
        kind = monitor.get("type")
        if kind not in {"distance", "angle", "dihedral"}:
            raise ValueError("Monitor type must be distance, angle, or dihedral.")
        masks = monitor.get("masks")
        required = {"distance": 2, "angle": 3, "dihedral": 4}[kind]
        if not isinstance(masks, list) or len(masks) != required:
            raise ValueError(f"{name}: {kind} requires an ordered list of {required} masks.")
        for mask in masks:
            safe_text(mask, "monitor mask")
        if monitor.get("center", "geometry") not in {"mass", "geometry"}:
            raise ValueError("Monitor center must be mass or geometry.")
        for option, applicable in (("image", "distance"), ("range360", "dihedral")):
            if option in monitor and (kind != applicable or type(monitor[option]) is not bool):
                raise ValueError(f"{option} is a boolean option for {applicable} monitors only.")
        if kind not in enabled:
            raise ValueError(f"Monitor {name} requires {kind} in analysis.additional or analysis.enabled.")
    for kind in ("angle", "dihedral"):
        if kind in enabled and not any(m["type"] == kind for m in config["monitors"]):
            raise ValueError(f"Add at least one {kind} monitor when enabling {kind}.")
    frames = config["frames"]
    for key in ("start", "stride"):
        if type(frames[key]) is not int or frames[key] < 1:
            raise ValueError(f"frames.{key} must be a positive integer.")
    if type(frames["stop"]) is not int or (frames["stop"] != -1 and frames["stop"] < frames["start"]):
        raise ValueError("frames.stop must be -1 or >= frames.start.")
    if frames["dt_ps"] is not None and (type(frames["dt_ps"]) not in (int, float) or not math.isfinite(frames["dt_ps"]) or frames["dt_ps"] <= 0):
        raise ValueError("frames.dt_ps must be positive or null (saved-frame interval, not the integration timestep).")
    from .convergence_config import validate as validate_convergence
    validate_convergence(config["convergence"], base, frames["dt_ps"], check_files)
    advanced = config["advanced"]
    if advanced["cluster_metric"] not in {"rms", "dme", "srmsd", "features"}:
        raise ValueError("cluster_metric must be rms, dme, srmsd, or features")
    if advanced["cluster_algorithm"] not in {"kmeans", "hieragglo"}:
        raise ValueError("cluster_algorithm must be kmeans or hieragglo")
    if not isinstance(advanced["cluster_features"], list):
        raise ValueError("cluster_features must list monitor/weight objects")
    seen_features = set()
    for feature in advanced["cluster_features"]:
        if not isinstance(feature, dict) or set(feature) != {"monitor", "weight"}:
            raise ValueError("Each cluster feature requires monitor and weight")
        name, weight = feature["monitor"], feature["weight"]
        if not isinstance(name, str) or name not in {m["name"] for m in config["monitors"]} or name in seen_features:
            raise ValueError("Cluster features must name distinct configured monitors")
        seen_features.add(name)
        if type(weight) not in (int, float) or not math.isfinite(weight) or weight <= 0:
            raise ValueError("Cluster feature weights must be finite and positive")
    if advanced["cluster_metric"] == "features" and not advanced["cluster_features"]:
        raise ValueError("Feature clustering requires at least one geometry monitor")
    for key in ("pca_modes", "cluster_count", "cluster_sieve", "max_frames", "max_matrix_atoms"):
        if type(config["advanced"][key]) is not int or config["advanced"][key] < 1:
            raise ValueError(f"advanced.{key} must be a positive integer.")
    if type(config["advanced"]["max_memory_gb"]) not in (int, float) or not 0 < config["advanced"]["max_memory_gb"] < 100000:
        raise ValueError("advanced.max_memory_gb must be positive and finite.")
    if type(config["advanced"]["pairwise_rmsd"]) is not bool:
        raise ValueError("advanced.pairwise_rmsd must be true or false.")
    if type(config["stats"]["block_size"]) is not int or config["stats"]["block_size"] < 2:
        raise ValueError("stats.block_size must be an integer >= 2.")
    if type(config["plots"]["alpha"]) not in (int, float) or not 0 < config["plots"]["alpha"] <= 1:
        raise ValueError("plots.alpha must be in (0, 1].")
    if type(config["plots"]["dpi"]) is not int or config["plots"]["dpi"] < 1:
        raise ValueError("plots.dpi must be a positive integer.")
    if type(config["plots"]["enabled"]) is not bool or type(config["stats"]["metric_correlations"]) is not bool:
        raise ValueError("plots.enabled and stats.metric_correlations must be true or false.")
    if type(config["reports"]["max_matrix_cells"]) is not int or config["reports"]["max_matrix_cells"] < 1:
        raise ValueError("reports.max_matrix_cells must be a positive integer.")
    for key, minimum in (("preview_series", 1), ("preview_points", 8)):
        if type(config["reports"][key]) is not int or config["reports"][key] < minimum:
            raise ValueError(f"reports.{key} must be an integer >= {minimum}.")
    for collection_name, mask_keys, extra_keys in (("sections", ("mask",), set()), ("interactions", ("mask1", "mask2"), {"contact_cutoff"})):
        collection = config[collection_name]
        if not isinstance(collection, list):
            raise ValueError(f"{collection_name} must be a list.")
        collection_names = set()
        for entry in collection:
            if not isinstance(entry, dict):
                raise ValueError(f"Each {collection_name} entry must be an object.")
            unknown = set(entry) - {"name", *mask_keys, *extra_keys}
            if unknown:
                raise ValueError(f"Unknown {collection_name} fields: {', '.join(sorted(unknown))}")
            name = safe_name(entry.get("name"))
            if name.casefold() in collection_names:
                raise ValueError(f"Duplicate selection/interaction name: {name}")
            collection_names.add(name.casefold())
            for key in mask_keys:
                safe_text(entry.get(key), key)
            if "contact_cutoff" in entry:
                cutoff = entry["contact_cutoff"]
                if type(cutoff) not in (int, float) or not math.isfinite(cutoff) or cutoff <= 0:
                    raise ValueError("interaction.contact_cutoff must be a finite positive number in Angstroms.")
    nucleic = config["nucleic"]
    if nucleic["resrange"] is not None and (not isinstance(nucleic["resrange"], str) or not re.fullmatch(r"[1-9]\d*(?:-[1-9]\d*)?(?:,[1-9]\d*(?:-[1-9]\d*)?)*", nucleic["resrange"])):
        raise ValueError("nucleic.resrange must be a positive topology residue range, e.g. 1-12,15-26.")
    if not isinstance(nucleic["resmap"], dict):
        raise ValueError("nucleic.resmap must map residue names to A, C, G, T or U.")
    for residue, base_name in nucleic["resmap"].items():
        if not isinstance(residue, str) or not re.fullmatch(r"[A-Za-z0-9_+\-]+", residue) or not isinstance(base_name, str) or base_name not in {"A", "C", "G", "T", "U"}:
            raise ValueError("nucleic.resmap must map residue names to A, C, G, T or U.")
    for value in (config["fit_mask"], config["advanced"]["matrix_mask"]):
        safe_text(value, "mask")
    if config["imaging"]["anchor"]:
        safe_text(config["imaging"]["anchor"], "imaging anchor")
    if check_files:
        for path in [config["topology"], *paths, *([config["reference"]] if config["reference"] else [])]:
            if not Path(path).is_file():
                raise ValueError(f"Input file not found: {path}")
    if check_files and "stacking" in enabled:
        from .rings import read_bond_graph
        atoms, _, _ = read_bond_graph(config["topology"])
        atom_residues = {a["index"]: a["residue"] for a in atoms}
        for ring in config["stacking"]["rings"]:
            if any(atom_residues.get(a) != ring["residue"] for a in ring["atoms"]):
                raise ValueError(f"Pi system {ring['name']}: atom indices do not belong to topology residue {ring['residue']}")
    output = Path(config["output"])
    if config["convergence"]["reference_path"] and Path(config["convergence"]["reference_path"]).is_relative_to(output):
        raise ValueError("Output directory must not contain convergence reference input data")
    for p in [config["topology"], *paths, *([config["reference"]] if config["reference"] else [])]:
        if Path(p).is_relative_to(output):
            raise ValueError("Output directory must not contain input data. Choose a separate results directory.")
    normalize_config_masks(config)
    return config


def resolve_selections(config):
    normalize_config_masks(config)
    info = inspect_topology(config["topology"])
    selections = dict(info["selections"])
    builtin = set(selections)
    def resolve(mask):
        return normalize_mask(selections.get(mask, mask))
    for section in config["sections"]:
        if section["name"] in builtin and section["mask"] != section["name"] and resolve(section["mask"]) != selections[section["name"]]:
            raise ValueError(f"Section {section['name']} would overwrite a built-in selection; choose another name.")
        selections[section["name"]] = resolve(section["mask"])
    if not config["sections"]:
        suggested = [name for name in ("complex", "backbone", "nucleic_backbone") if name in selections]
        if not suggested:
            raise ValueError("No protein/nucleic sections detected. Add explicit sections and fit_mask to the JSON config.")
        config["sections"] = [{"name": name, "mask": selections[name]} for name in suggested]
    else:
        config["sections"] = [{**s, "mask": resolve(s["mask"])} for s in config["sections"]]
    config["fit_mask"] = resolve(config["fit_mask"])
    config["advanced"]["matrix_mask"] = resolve(config["advanced"]["matrix_mask"])
    if config["advanced"]["matrix_mask"] == "representative" and "representative" not in selections:
        config["advanced"]["matrix_mask"] = config["fit_mask"]
        info["warnings"].append("No automatic representative atoms: matrix selection uses the explicitly supplied fit mask.")
    for pair in config["interactions"]:
        pair["mask1"], pair["mask2"] = resolve(pair["mask1"]), resolve(pair["mask2"])
    for monitor in config.get("monitors", []):
        monitor["masks"] = [resolve(mask) for mask in monitor["masks"]]
    for item in config.get("solvation", []):
        item["mask1"], item["mask2"] = resolve(item["mask1"]), resolve(item["mask2"])
        for mask in (item["mask1"], item["mask2"]):
            if not any(c in mask for c in ":@*^"):
                raise ValueError(f"Unknown solvation selection {mask!r}.")
    if config["imaging"]["anchor"]:
        config["imaging"]["anchor"] = resolve(config["imaging"]["anchor"])
    if not config["interactions"] and "protein" in selections and "nucleic" in selections:
        config["interactions"] = [{"name": "protein_nucleic", "mask1": selections["protein"], "mask2": selections["nucleic"]}]
    # An unresolved symbolic name otherwise produces a misleading CPPTRAJ error.
    for mask in [config["fit_mask"], config["advanced"]["matrix_mask"], *([config["imaging"]["anchor"]] if config["imaging"]["anchor"] else []), *[s["mask"] for s in config["sections"]], *[p[k] for p in config["interactions"] for k in ("mask1", "mask2")], *[mask for m in config.get("monitors", []) for mask in m["masks"]]]:
        if not any(c in mask for c in ":@*^"):
            raise ValueError(f"Unknown selection {mask!r}; use a detected selection name or explicit CPPTRAJ mask.")
    if config["imaging"]["mode"] == "auto":
        if info["periodic"] is None:
            raise ValueError("Periodicity is unknown for this topology. Set imaging.mode to on (box in trajectory) or off.")
        config["imaging"]["mode"] = "on" if info["periodic"] else "off"
    if any(s["type"] == "rdf" for s in config.get("solvation", [])) and config["imaging"]["mode"] != "on":
        raise ValueError("Volume-normalized RDF requires periodic trajectories; set imaging.mode=on only when a valid box is present.")
    return selections, info


def save_config(config, path):
    Path(path).write_text(json.dumps(config, indent=2) + "\n", encoding="utf-8")
