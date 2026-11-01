"""Explicitly synthetic fixtures for exercising reports without running cpptraj."""
from pathlib import Path
import json

import numpy as np

from .reporting import report_results


def create_demo(output_dir: str | Path, render: bool = True) -> dict:
    """Create three unequal synthetic replicas and run the standard reporter.

    No topology or coordinate simulation is generated; these are artificial
    analysis observables for demonstrating formatting and statistical behavior.
    """
    output_dir = Path(output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    marker = output_dir / "SYNTHETIC_DEMO_README.md"
    if any(output_dir.iterdir()) and not marker.is_file():
        raise ValueError("Demo output directory must be empty, or an existing marked synthetic demo directory")
    marker.write_text("# SYNTHETIC DEMONSTRATION — NOT MOLECULAR DYNAMICS RESULTS\n\nAll observables in this directory are generated from simple mathematical functions and seeded random numbers. No molecular system was simulated, and no cpptraj command was executed. Use these files to inspect output formats, plotting, unequal replica lengths, circular-angle handling, and matrix summaries.\n", encoding="utf-8")
    rng = np.random.default_rng(80211)
    config = {"frames": {"start": 1, "stop": -1, "stride": 1, "dt_ps": 20}, "plots": {"enabled": render, "dpi": 300, "alpha": 0.8}, "stats": {"block_size": 30, "metric_correlations": True}}
    config["study"] = {"title": "Synthetic report preview", "description": "DESIGN DEMONSTRATION — artificial observables, not molecular-dynamics results."}
    manifest = {"synthetic": True, "description": "SYNTHETIC mathematical observables; not molecular dynamics results", "replicas": []}
    for index, length in enumerate((300, 240, 360), 1):
        replica_name = f"SYNTHETIC_replica_{index}"
        raw_dir = output_dir / "raw" / replica_name
        raw_dir.mkdir(parents=True, exist_ok=True)
        frame = np.arange(1, length + 1)
        smooth_noise = np.convolve(rng.normal(0, 0.14, length + 10), np.ones(11) / 11, mode="valid")
        rmsd = 1.1 + 0.1 * index + 0.5 * (1 - np.exp(-frame / 45)) + 0.15 * np.sin(frame / 22 + index) + smooth_noise
        rg = 18 + 0.04 * index + 0.2 * np.cos(frame / 31) + rng.normal(0, 0.02, length)
        artifacts = []

        def write_table(filename, values, header, kind, analysis, section, unit, metadata=None):
            path = raw_dir / filename
            np.savetxt(path, values, fmt="%.8g", header=header, comments="#")
            artifacts.append({"path": str(path), "kind": kind, "analysis": f"SYNTHETIC {analysis}", "section": section, "unit": unit, "metadata": metadata or {}})

        write_table("rmsd.dat", np.column_stack([frame, rmsd]), "Frame RMSD", "timeseries", "RMSD", "protein backbone", "angstrom")
        write_table("rg.dat", np.column_stack([frame, rg]), "Frame Rg", "timeseries", "radius of gyration", "protein", "angstrom")
        residues = np.arange(1, 25)
        rmsf = 0.35 + 0.8 * np.exp(-((residues - 12) / 3) ** 2) + rng.uniform(0, 0.12, len(residues))
        write_table("rmsf.dat", np.column_stack([residues, rmsf]), "Res RMSF", "profile", "RMSF", "protein backbone", "angstrom")
        phi = (175 + 12 * np.sin(frame / 16 + index) + rng.normal(0, 2, length) + 180) % 360 - 180
        write_table("torsion.dat", np.column_stack([frame, phi]), "Frame Phi_12", "timeseries", "torsion", "protein residue 12", "degree", {"circular": True})
        state = ((frame // (35 + index * 5)) % 3).astype(int)
        write_table("states.dat", np.column_stack([frame, state]), "Frame State_12", "timeseries", "categorical states", "synthetic residue 12", "state", {"categorical": True})
        latent = rng.normal(size=(length, 2))
        loadings = np.column_stack([np.cos(residues / 4), np.sin(residues / 4)])
        signals = latent @ loadings.T + rng.normal(0, 0.35, (length, len(residues)))
        matrix = np.corrcoef(signals, rowvar=False)
        write_table("correlation.dat", matrix, "SYNTHETIC dense correlation matrix", "matrix", "correlation", "protein C-alpha", "correlation", {"matrix_format": "dense", "matrix_type": "correlation", "x_labels": residues.tolist(), "y_labels": residues.tolist(), "x_unit": "residue", "y_unit": "residue"})
        manifest["replicas"].append({"name": replica_name, "status": "synthetic", "frame_count": length, "artifacts": artifacts})
    (output_dir / "synthetic_config.json").write_text(json.dumps(config, indent=2) + "\n", encoding="utf-8")
    (output_dir / "synthetic_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    result = report_results(config, manifest, output_dir)
    result["synthetic"] = True
    result["description"] = "SYNTHETIC demonstration; no molecular dynamics or cpptraj execution"
    (output_dir / "reports" / "report_summary.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    return result
