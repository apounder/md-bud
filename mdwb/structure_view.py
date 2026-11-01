"""Bounded reference-structure embedding; no trajectory frames enter the browser."""
from pathlib import Path


def structure_payload(config, manifest, output):
    options = config.get("structure_view", {})
    if not options.get("enabled", True):
        return {"status": "disabled", "message": "Structure preview was disabled in the configuration."}
    path = Path(output) / "reference" / "viewer.pdb"
    record = manifest.get("structure_view", {})
    if record.get("status") not in {None, "complete", "synthetic", "example"}:
        return {"status": record.get("status"), "message": record.get("message", "Structure preview unavailable.")}
    if not path.is_file():
        return {"status": "unavailable", "message": "No reference PDB was exported for this run; numerical results remain available."}
    if path.stat().st_size > options.get("max_bytes", 1500000):
        return {"status": "size_limit", "message": "Reference structure exceeds the configured embedding size limit."}
    lines, atoms = [], 0
    with path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("ENDMDL"):
                break
            if line.startswith(("ATOM  ", "HETATM")):
                atoms += 1
                if atoms > options.get("max_atoms", 12000):
                    return {"status": "atom_limit", "message": "Reference structure exceeds the configured atom limit."}
                lines.append(line.rstrip())
            elif line.startswith("TER"):
                lines.append("TER")
    if not atoms:
        return {"status": "empty", "message": "No readable atoms found in the reference export."}
    message = "One common reference frame, selected heavy solute atoms; displayed residue numbers may change after stripping."
    if record.get("status") == "synthetic":
        message = "Synthetic geometry for design demonstration."
    elif record.get("status") == "example":
        message = record.get("message", "Example structure, unrelated to the numerical observables.")
    return {"status": "complete", "pdb": "\n".join(lines) + "\nEND\n", "atoms": atoms,
            "message": message,
            "selection": record.get("selection", "exported reference")}
