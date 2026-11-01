"""Environment checks, provenance, methods notes and portable result archives."""
from __future__ import annotations

import copy
import importlib.metadata
import json
import platform
import re
import shutil
import sys
from pathlib import Path, PurePosixPath, PureWindowsPath
from zipfile import ZIP_DEFLATED, ZipFile

from . import __version__
from .analysis import SOURCES, build_batches, build_pooled_batches, enabled_analyses
from .config import resolve_selections


def environment_info(cpptraj="cpptraj"):
    packages = {}
    for name in ("numpy", "matplotlib"):
        try:
            packages[name] = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            packages[name] = None
    binary = shutil.which(cpptraj)
    if binary is None and Path(cpptraj).is_file():
        binary = str(Path(cpptraj).resolve())
    minimum = {"numpy": (1, 24), "matplotlib": (3, 7)}
    compatible = {}
    for name, required in minimum.items():
        match = re.match(r"^(\d+)\.(\d+)", packages[name] or "")
        compatible[name] = bool(match and tuple(map(int, match.groups())) >= required)
    return {"workbench": __version__, "python": platform.python_version(),
            "python_executable": sys.executable, "platform": platform.platform(),
            "packages": packages, "meets_minimum_version": compatible, "cpptraj_executable": binary,
            "python_supported": sys.version_info >= (3, 10),
            "note": "Package presence/version only; no native engine or numerical calculation was executed."}


def inspect_config(config):
    cfg = copy.deepcopy(config)
    selections, info = resolve_selections(cfg)
    batches = build_batches(cfg, selections)
    pooled = build_pooled_batches(cfg, selections)
    implemented = [b["name"] for b in batches + pooled]
    return {"status": "configuration_valid_engine_not_checked", "topology": cfg["topology"],
            "replicas": [{"name": r["name"], "segments": len(r["trajectories"])} for r in cfg["replicas"]],
            "natom": info["natom"], "fit_mask": cfg["fit_mask"], "imaging": cfg["imaging"],
            "sections": cfg["sections"], "monitors": cfg.get("monitors", []), "solvation": cfg.get("solvation", []),
            "requested": sorted(enabled_analyses(cfg)), "applicable_batches": implemented,
            "warnings": info["warnings"],
            "notes": ["File existence and Python-side configuration checked; atom correspondence and native mask readability still need CPPTRAJ/user review.",
                      "Inapplicable specialist analyses are omitted from applicable_batches."]}


def write_methods(config, manifest, output):
    output = Path(output)
    environment = manifest.setdefault("environment", environment_info(config.get("cpptraj", "cpptraj")))
    environment["cpptraj_version"] = manifest.get("cpptraj_version")
    (output / "environment.json").write_text(json.dumps(environment, indent=2) + "\n", encoding="utf-8")
    frames = config.get("frames", {})
    title = config.get("study", {}).get("title", "MD analysis")
    lines = ["# Methods record: " + title, "", "Generated record for review before use in a manuscript. Execution status: " + manifest.get("status", "unknown") + ".", "",
             f"Analysis was orchestrated with MD-bud {manifest.get('workbench_version', __version__)}.",
             "The exact CPPTRAJ version is retained in cpptraj_version.log; Python/package/platform versions are in environment.json.", "",
             f"Frame selection: start={frames.get('start', 1)}, stop={frames.get('stop', -1)}, stride={frames.get('stride', 1)} across joined segments within each replica.",
             f"Saved source-frame interval (ps): {frames.get('dt_ps')}; null means no inferred physical time.",
             f"Global fit mask: `{config.get('fit_mask', '')}`. Imaging: `{config.get('imaging', {})}`.",
             "All comparisons use the common prepared reference. Section global RMSD retains relative motion; local RMSD fits each section without changing subsequent coordinates.", "", "## Actual batch outcomes", ""]
    for rep in manifest.get("replicas", []):
        lines.append(f"- {rep['name']}: {rep.get('frame_count', 'unknown')} selected frames; " +
                     ", ".join(f"{b['name']}={b['status']}" for b in rep.get("batches", [])))
    notes = config.get("notes", {})
    if isinstance(notes, dict) and notes.get("campaign_manifest_sha256"):
        lines += ["", "## Campaign provenance", "",
                  "Source campaign manifest SHA-256: `" + notes["campaign_manifest_sha256"] + "`.",
                  f"Expected production frames per segment: {notes.get('campaign_expected_frames_per_segment', 'unspecified')}.",
                  "Only manifest-listed production segments were selected; replica boundaries were retained.",
                  notes.get("replica_caution", "")]
    for job in manifest.get("pooled", []):
        lines.append(f"- Shared {job['name']}: {job['status']}")
    lines += ["", "## Selections and parameters", "", "```json", json.dumps({key: config.get(key) for key in ("sections", "interactions", "monitors", "nucleic", "advanced", "diagnostics", "convergence", "stacking", "stats")}, indent=2), "```", "",
              "## Interpretation", "", "Independent replicas remain the sampling units. Replica summaries weight available replica means equally; pooled PCA/clustering bases weight input frames. Appended plots do not imply temporal continuity between replicas. Failed analyses must not be described as successfully completed.", "",
              "Sampling diagnostics, when enabled, assume stationarity and do not select an equilibration cutoff. Histogram Jensen–Shannon distances describe observed distribution differences; they are not hypothesis tests or convergence certificates.", "",
              "## References", "", "- Roe & Cheatham (2013), PTRAJ and CPPTRAJ: https://doi.org/10.1021/ct400341p",
              "- Diagnostic estimation: Geyer (1992), https://doi.org/10.1214/ss/1177011137",
              "- Ensemble-comparison context: PENSA, https://github.com/drorlab/pensa",
              "- Workbench software citation: record this package version, real authors/maintainers and your eventual repository/archive identifier; no DOI has been assigned by this tool.", "",
              "Do not cite MDAnalysis, MDTraj, ProLIF, PyMBAR or PENSA as executed dependencies: this workbench uses CPPTRAJ and independent NumPy reporting code. Consult their papers when discussing the reviewed methods or using those tools separately.", "", "Command references:"]
    for key in sorted(enabled_analyses(config)):
        if key in SOURCES:
            lines.append(f"- {key}: {SOURCES[key]}")
    if config.get("solvation"):
        lines += ["", "## Solvation definitions", "", "These batches preserve the imaged coordinate/box orientation (RMS check uses nomod). RDF uses average-volume normalization; users verify the maximum radius against the periodic cell. Shell counts are geometric criteria, not residence times.", "", "```json", json.dumps(config["solvation"], indent=2), "```"]
    (output / "METHODS.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def _relative_artifact(value, old_output):
    # Interpret original path syntax even when a result moves between operating systems.
    cls = PureWindowsPath if PureWindowsPath(str(old_output)).is_absolute() else PurePosixPath
    path, old = cls(str(value)), cls(str(old_output))
    if path.is_absolute():
        try:
            path = path.relative_to(old)
        except ValueError:
            raise ValueError("Artifact lies outside the original results directory.")
    if path.anchor or ".." in path.parts:
        raise ValueError("Artifact path escapes the results directory.")
    return Path(*path.parts)


def rebase_manifest(manifest, old_output, new_output):
    """Allow reporting after moving a results directory, without touching inputs."""
    result = copy.deepcopy(manifest)
    new = Path(new_output).resolve()
    for entry in result.get("replicas", []) + result.get("pooled", []):
        for artifact in entry.get("artifacts", []):
            artifact["path"] = str(new / _relative_artifact(artifact["path"], old_output))
    return result


def bundle_results(directory, destination):
    """Bundle generated results and relative artifact paths; never follow symlinks."""
    directory, destination = Path(directory).resolve(), Path(destination).resolve()
    if destination.exists():
        raise ValueError("Bundle destination already exists; choose a new ZIP filename.")
    if destination.is_relative_to(directory):
        raise ValueError("Put the archive outside the results directory.")
    config = json.loads((directory / "resolved_config.json").read_text(encoding="utf-8"))
    manifest = json.loads((directory / "manifest.json").read_text(encoding="utf-8"))
    original_root = config["output"]
    for entry in manifest.get("replicas", []) + manifest.get("pooled", []):
        for artifact in entry.get("artifacts", []):
            path = _relative_artifact(artifact["path"], original_root)
            artifact["path"] = path.as_posix()
    config["output"] = "."
    input_paths = {Path(p).resolve() for p in [config["topology"], *([config["reference"]] if config.get("reference") else []),
                                            *[p for r in config["replicas"] for p in r["trajectories"]]]}
    files = []
    for path in directory.rglob("*"):
        if path.is_symlink():
            raise ValueError(f"Symlinks are not bundled: {path}")
        if path.is_file():
            if not path.resolve().is_relative_to(directory) or path.resolve() in input_paths:
                raise ValueError("Results directory includes an original input or an external file.")
            files.append(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    with ZipFile(destination, "x", compression=ZIP_DEFLATED) as archive:
        for path in files:
            relative = path.relative_to(directory).as_posix()
            if relative == "BUNDLE_README.txt":
                continue
            if relative == "resolved_config.json":
                archive.writestr(relative, json.dumps(config, indent=2) + "\n")
            elif relative == "manifest.json":
                archive.writestr(relative, json.dumps(manifest, indent=2) + "\n")
            else:
                archive.write(path, relative)
        archive.writestr("BUNDLE_README.txt", "Open index.html to browse results offline. Run mdworkbench report <extracted-directory> to regenerate reports.\nOriginal trajectories/topology are not included. This is a results archive, not a self-contained simulation.\nOriginal input paths and metadata remain in provenance and logs; review them before sharing.\n")
        file_count = len(archive.infolist())
    return {"archive": str(destination), "files": file_count, "original_trajectories_included": False}
