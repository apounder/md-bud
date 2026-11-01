"""Conservative Amber/PDB topology inspection; CPPTRAJ remains the final authority."""
from __future__ import annotations

import gzip
import re
from pathlib import Path

PROTEIN = set("ALA ARG ASN ASP ASH CYS CYM CYX GLN GLU GLH GLY HIS HID HIE HIP ILE LEU LYS LYN MET PHE PRO SER THR TRP TYR VAL HYP MSE".split())
DNA = set("DA DC DG DT DI DU ADE CYT GUA THY".split())
RNA = set("A C G U I RA RC RG RU RI URA".split())
WATER = set("WAT HOH SOL TIP3 TIP4 TIP5 OPC SPC SPCE".split())
IONS = set("NA NA+ CL CL- K K+ MG MG2+ CA CA2+ ZN ZN2+ CS LI RB F BR IOD SOD CLA POT CAL".split())


def residue_range(indices):
    """Compact a sorted set of one-based topology residue indices."""
    vals = sorted(set(indices))
    chunks = []
    if not vals:
        return ""
    start = end = vals[0]
    for n in vals[1:]:
        if n == end + 1:
            end = n
        else:
            chunks.append(str(start) if start == end else f"{start}-{end}")
            start = end = n
    chunks.append(str(start) if start == end else f"{start}-{end}")
    return ",".join(chunks)


def _amber_flags(path, extra_flags=()):
    result = {}
    key = None
    width = 0
    keep = {"POINTERS", "RESIDUE_LABEL", "RESIDUE_POINTER", "ATOM_NAME"} | set(extra_flags)
    opener = gzip.open if Path(path).suffix.lower() == ".gz" else open
    with opener(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("%FLAG "):
                key = line.split()[1]
                width = 0
                if key in keep:
                    result[key] = []
            elif line.startswith("%FORMAT"):
                m = re.search(r"\d+[aAiIeEfF](\d+)", line)
                width = int(m.group(1)) if m else 0
            elif key in keep and width and not line.startswith("%"):
                raw = line.rstrip("\r\n")
                result[key].extend(raw[i:i + width].strip() for i in range(0, len(raw), width) if raw[i:i + width].strip())
    return result


def inspect_topology(path):
    path = Path(path)
    opener = gzip.open if path.suffix.lower() == ".gz" else open
    format_path = Path(path.stem) if path.suffix.lower() == ".gz" else path
    with opener(path, "rb") as handle:
        head = handle.read(256)
    residues = []
    natom = None
    periodic = None
    if b"%VERSION" in head or b"%FLAG" in head:
        flags = _amber_flags(path)
        labels = flags.get("RESIDUE_LABEL", [])
        pointers = [int(n) for n in flags.get("RESIDUE_POINTER", [])]
        atoms = flags.get("ATOM_NAME", [])
        general = [int(n) for n in flags.get("POINTERS", [])]
        if not labels or len(labels) != len(pointers) or not general:
            raise ValueError("Incomplete Amber topology: cannot read residue labels/pointers.")
        natom = general[0]
        if natom < 1 or len(atoms) != natom or pointers[0] != 1 or pointers != sorted(set(pointers)) or pointers[-1] > natom:
            raise ValueError("Inconsistent Amber topology: atom count or residue pointers do not match atom names.")
        if len(general) > 11 and general[11] != len(labels):
            raise ValueError("Inconsistent Amber topology: NRES does not match residue labels.")
        if len(general) > 27:
            periodic = bool(general[27])
        for i, (label, first) in enumerate(zip(labels, pointers)):
            end = pointers[i + 1] - 1 if i + 1 < len(pointers) else natom
            residues.append({"index": i + 1, "name": label, "atoms": atoms[first - 1:end]})
        fmt = "amber"
    elif format_path.suffix.lower() in {".pdb", ".ent"}:
        last_key = None
        with opener(path, "rt", encoding="utf-8", errors="replace") as handle:
            for line in handle:
                if line.startswith("ENDMDL"):
                    break
                if line.startswith("TER"):
                    last_key = None
                if line.startswith(("ATOM  ", "HETATM")):
                    key = (line[21:22], line[22:27], line[17:20])
                    if key != last_key:
                        residues.append({"index": len(residues) + 1, "name": line[17:20].strip(), "atoms": []})
                        last_key = key
                    residues[-1]["atoms"].append(line[12:16].strip())
        natom = sum(len(r["atoms"]) for r in residues)
        fmt = "pdb"
    else:
        return {"format": "cpptraj", "natom": None, "periodic": None, "residues": [], "selections": {}, "warnings": ["Automatic residue detection is available for Amber parm7/prmtop and PDB. Supply explicit masks for this topology format; CPPTRAJ will validate it."]}
    groups = {key: [] for key in ("protein", "dna", "rna", "water", "ions", "other")}
    inferred = []
    for residue in residues:
        label = residue["name"].upper()
        names = {a.upper().replace("*", "'") for a in residue["atoms"]}
        amino_label = label[1:] if len(label) == 4 and label[0] in "NC" else label
        na_label = label[:-1] if label.endswith(("3", "5")) else label
        if amino_label in PROTEIN:
            kind = "protein"
        elif na_label in DNA:
            kind = "dna"
        elif na_label in RNA:
            kind = "rna"
        elif label in WATER:
            kind = "water"
        elif label in IONS and len(names) == 1:
            kind = "ions"
        elif {"N", "CA", "C", "O"} <= names:
            kind = "protein"
            inferred.append(residue["index"])
        elif {"C1'", "C2'", "C3'", "C4'", "O4'"} <= names and ("N1" in names or "N9" in names):
            kind = "rna" if "O2'" in names else "dna"
            inferred.append(residue["index"])
        else:
            kind = "other"
        residue["kind"] = kind
        groups[kind].append(residue["index"])
    groups["nucleic"] = groups["dna"] + groups["rna"]
    groups["complex"] = groups["protein"] + groups["nucleic"]
    groups["solute"] = groups["complex"] + groups["other"]
    selections = {key: ":" + residue_range(vals) for key, vals in groups.items() if vals}
    if "protein" in selections:
        selections["backbone"] = f"({selections['protein']})&@N,CA,C,O"
        selections["ca"] = f"({selections['protein']})&@CA"
    if "nucleic" in selections:
        selections["nucleic_backbone"] = f"({selections['nucleic']})&@P,O5',C5',C4',C3',O3'"
        selections["bases"] = f"({selections['nucleic']})&!@P,OP1,OP2,OP3,O1P,O2P,O3P,O5',C5',C4',O4',C3',O3',C2',O2',C1'&!@/H"
    # A literal '*' in old sugar atom names is a CPPTRAJ mask wildcard.
    # Exact topology atom indices avoid broadening the selection accidentally.
    legacy_nucleic = any("*" in atom for residue in residues if residue["index"] in groups["nucleic"] for atom in residue["atoms"])
    legacy_representative = None
    if legacy_nucleic:
        backbone_names = {"P", "O5'", "C5'", "C4'", "C3'", "O3'"}
        excluded_base_names = backbone_names | {"OP1", "OP2", "OP3", "O1P", "O2P", "O3P", "O4'", "C2'", "O2'", "C1'"}
        atom_groups = {"backbone": [], "bases": [], "representative": []}
        atom_index = 0
        for residue in residues:
            for atom in residue["atoms"]:
                atom_index += 1
                if residue["index"] not in groups["nucleic"]:
                    continue
                canonical = atom.upper().replace("*", "'")
                if canonical in backbone_names:
                    atom_groups["backbone"].append(atom_index)
                if canonical not in excluded_base_names:
                    atom_groups["bases"].append(atom_index)
                if canonical == "C4'":
                    atom_groups["representative"].append(atom_index)
        selections["nucleic_backbone"] = "@" + (residue_range(atom_groups["backbone"]) or "0")
        selections["bases"] = "(@" + (residue_range(atom_groups["bases"]) or "0") + ")&!@/H"
        legacy_representative = "@" + (residue_range(atom_groups["representative"]) or "0")
    representatives = []
    if "ca" in selections:
        representatives.append(f"({selections['ca']})")
    if "nucleic" in selections:
        representatives.append(f"({legacy_representative})" if legacy_representative else f"(({selections['nucleic']})&@C4')")
    if representatives:
        selections["representative"] = "|".join(representatives)
    if "solute" in selections:
        selections["heavy"] = f"({selections['solute']})&!@/H"
    warnings = []
    if legacy_nucleic:
        warnings.append("Legacy nucleic '*' atom names: automatic sugar/base masks use topology atom numbers to preserve exact atom identity.")
    if inferred:
        warnings.append("Modified residues inferred from backbone atom names; review: " + residue_range(inferred))
    if groups["other"]:
        warnings.append("Unclassified residues (ligands/cofactors/modifications): " + residue_range(groups["other"]))
    if fmt == "pdb":
        warnings.append("PDB topology bonding/atom types are incomplete; check CPPTRAJ's interpretation before chemical analyses.")
    return {"format": fmt, "natom": natom, "periodic": periodic, "residues": residues, "selections": selections, "warnings": warnings}
