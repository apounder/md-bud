"""Presentation labels; native identifiers and filenames remain unchanged."""
import re

ANALYSIS_LABELS = {
    "stacking": "Pi-stacking occupancy", "stacking_any": "Pi-stacking: any partner",
    "stacking_distance": "Pi-system centroid distance", "stacking_angle": "Pi-system plane angle",
    "rg": "Radius of gyration", "rog": "Radius of gyration",
    "radius of gyration": "Radius of gyration",
    "rmsd": "Root-mean-square deviation (RMSD)",
    "rmsd_global": "RMSD · global alignment", "rmsd_local": "RMSD · local alignment",
    "rmsf": "Residue fluctuations (RMSF)", "sasa": "Solvent-accessible surface area",
    "hbond": "Hydrogen bonds", "contacts": "Residue contacts",
    "distance": "Distance", "angle": "Angle", "dihedral": "Dihedral angle",
    "torsion": "Torsion angle", "torsions": "Torsion angles",
    "dssp": "Protein secondary structure", "nastruct": "Nucleic-acid structure",
    "rdf": "Radial distribution function", "watershell": "Solvent-shell occupancy",
    "pca": "Principal component analysis", "cluster": "Conformational clusters",
    "dccm": "Dynamic cross-correlation", "correlation": "Correlation matrix",
    "pairwise_rmsd": "Pairwise RMSD", "categorical states": "State populations",
    "monitor_distance": "Distance monitor", "monitor_angle": "Angle monitor",
    "monitor_dihedral": "Dihedral-angle monitor",
}
UNITS = {"angstrom": "Å", "angstrom^2": "Å²", "degree": "°", "ps": "ps", "ns": "ns"}


def readable(value):
    text = re.sub(r"_+", " ", str(value)).strip()
    return text[:1].upper() + text[1:]


def analysis_label(value):
    raw = str(value)
    synthetic = raw.lower().startswith("synthetic ")
    key = raw[10:] if synthetic else raw
    label = ANALYSIS_LABELS.get(key.lower(), readable(key))
    if key.lower() not in ANALYSIS_LABELS and "_" in key:
        base, suffix = key.lower().split("_", 1)
        if base in ANALYSIS_LABELS:
            label = f"{ANALYSIS_LABELS[base]} · {readable(suffix)}"
    return label + (" [synthetic]" if synthetic else "")


def observable_label(series):
    name = series.name
    if name.startswith("A_RMS_"):
        return "RMSD"
    if name.startswith("A_MON_"):
        return readable(series.section)
    # Expand familiar CPPTRAJ scalar headings while keeping monitor identities.
    aliases = {"rg": "Radius of gyration", "rog": "Radius of gyration",
               "rmsd": "RMSD", "rmsf": "RMSF", "sasa": "Surface area"}
    automatic = re.match(r"^A_(RG|RMSD|RMSF|SASA)_", name, flags=re.IGNORECASE)
    if automatic:
        return aliases[automatic[1].lower()]
    return aliases.get(name.lower(), readable(name))


def series_title(series, context=""):
    heading = analysis_label(series.analysis)
    observable = observable_label(series)
    detail = readable(series.section)
    if observable.lower() not in heading.lower() and observable.lower() != detail.lower():
        detail += f" · {observable}"
    return f"{heading}\n{detail}" + (f"\n{context}" if context else "")


def unit_label(unit):
    return UNITS.get(unit, unit)


def coordinate_label(unit):
    return f"Time ({unit})" if unit in {"ps", "ns", "fs"} else readable(unit)
