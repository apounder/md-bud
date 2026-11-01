"""Build reviewable CPPTRAJ analysis batches without running external programs.

The runner loads a topology, trajectories and the common reference [COMMON],
images coordinates if configured, and performs the global RMS fit first. It
then appends ``commands``, ``run``, and ``post_commands`` in that order.

PCA and clustering belong to a *single pooled job*. Independent eigenvectors
or independently numbered clusters cannot be overlaid across replicas.
"""
from __future__ import annotations

import hashlib
import math
import re
from typing import Any


SOURCES = {
    "stacking": "https://amberhub.chpc.utah.edu/vector/",
    "rdf": "https://amberhub.chpc.utah.edu/radial-rdf/",
    "watershell": "https://amberhub.chpc.utah.edu/watershell/",
    "rmsd": "https://amberhub.chpc.utah.edu/rmsd/",
    "rg": "https://amberhub.chpc.utah.edu/radgyr-rog/",
    "rmsf": "https://amberhub.chpc.utah.edu/atomicfluct-rmsf/",
    "sasa": "https://amberhub.chpc.utah.edu/surf/",
    "hbond": "https://amberhub.chpc.utah.edu/hbond/",
    "dssp": "https://amberhub.chpc.utah.edu/secstruct/",
    "nastruct": "https://amberhub.chpc.utah.edu/nastruct/",
    "contacts": "https://amberhub.chpc.utah.edu/nativecontacts/",
    "distance": "https://amberhub.chpc.utah.edu/distance/",
    "angle": "https://amberhub.chpc.utah.edu/angle/",
    "dihedral": "https://amberhub.chpc.utah.edu/dihedral/",
    "torsions": "https://amberhub.chpc.utah.edu/multidihedral/",
    "dccm": "https://amberhub.chpc.utah.edu/matrix-2/",
    "pca": "https://amberhub.chpc.utah.edu/introduction-to-principal-component-analysis/",
    "cluster": "https://amberhub.chpc.utah.edu/cluster/",
    "pairwise_rmsd": "https://amberhub.chpc.utah.edu/rms2d/",
}

BASIC = frozenset({"rmsd", "rg", "rmsf"})
STANDARD = BASIC | {"hbond", "dssp", "distance", "contacts"}
ADVANCED = STANDARD | {"sasa", "nastruct", "torsions", "dccm", "pca", "cluster"}
KNOWN_ANALYSES = ADVANCED | {"stacking", "pairwise_rmsd", "angle", "dihedral", "rdf", "watershell"}
LEVELS = {"basic": BASIC, "standard": STANDARD, "advanced": ADVANCED}


def safe_name(value: str) -> str:
    """Stable file-safe name; preserve distinct labels after normalization."""
    value = str(value)
    token = re.sub(r"[^A-Za-z0-9_]+", "_", value).strip("_") or "selection"
    if token != value or len(token) > 72:
        digest = hashlib.sha256(value.encode("utf-8")).hexdigest()[:8]
        token = token[:64] + "_" + digest
    return token


def quote_mask(mask: str) -> str:
    """Quote a CPPTRAJ mask, retaining nucleic-acid atom-name apostrophes."""
    if not isinstance(mask, str) or not mask.strip():
        raise ValueError("An atom selection must be a nonempty string.")
    if any(ord(char) < 32 for char in mask) or any(char in mask for char in '\\"#$'):
        raise ValueError("Atom masks cannot contain control characters, double quotes, \\, #, or $.")
    return '"' + mask.strip() + '"'


def resolve_mask(value: str, selections: dict[str, str]) -> str:
    """Resolve a named selection, or accept an explicit CPPTRAJ expression."""
    from .selections import normalize_mask
    value = normalize_mask(value)
    mask = normalize_mask(selections.get(value, value))
    quote_mask(mask)
    if not any(char in mask for char in ":@!*^/"):
        raise ValueError(f"Unknown selection {value!r}; use a named selection or CPPTRAJ mask.")
    return mask


def enabled_analyses(config: dict[str, Any]) -> set[str]:
    analysis = config.get("analysis", {})
    level = analysis.get("level", "standard")
    if level not in LEVELS:
        raise ValueError(f"Unknown analysis level {level!r}.")
    explicit = analysis.get("enabled")
    if explicit is not None:
        if not isinstance(explicit, list) or not all(isinstance(x, str) for x in explicit):
            raise ValueError("analysis.enabled must be a list of analysis names.")
        names = set(explicit)
        unknown = names - KNOWN_ANALYSES
        if unknown:
            raise ValueError("Unknown analyses: " + ", ".join(sorted(unknown)))
        return names
    names = set(LEVELS[level])
    if config.get("advanced", {}).get("pairwise_rmsd", False):
        names.add("pairwise_rmsd")
    additional = analysis.get("additional", [])
    if not isinstance(additional, list) or not all(isinstance(x, str) for x in additional):
        raise ValueError("analysis.additional must be a list of analysis names.")
    unknown = set(additional) - KNOWN_ANALYSES
    if unknown:
        raise ValueError("Unknown additional analyses: " + ", ".join(sorted(unknown)))
    names.update(additional)
    return names


def _positive(value: Any, name: str, integer: bool = False) -> float | int:
    number = float(value)
    if isinstance(value, bool) or not math.isfinite(number) or number <= 0:
        raise ValueError(f"{name} must be a finite positive number.")
    if integer:
        if not number.is_integer():
            raise ValueError(f"{name} must be an integer.")
        return int(number)
    return number


def _artifact(analysis: str, section: str, kind: str, unit: str,
              *, path: str | None = None, **metadata: Any) -> dict[str, Any]:
    return {"path": path or f"{analysis}__{safe_name(section)}.dat", "kind": kind,
            "analysis": analysis, "section": section, "unit": unit,
            "metadata": metadata}


def _batch(name: str, description: str, commands: list[str], artifacts: list[dict],
           post_commands: list[str] | None = None, **extra: Any) -> dict[str, Any]:
    return {"name": name, "description": description, "commands": commands,
            "artifacts": artifacts, "post_commands": post_commands or [], **extra}


def _solute(selections: dict[str, str]) -> str | None:
    for key in ("solute", "biopolymer"):
        if selections.get(key):
            return selections[key]
    parts = [selections[key] for key in ("protein", "nucleic") if selections.get(key)]
    return "|".join(f"({part})" for part in parts) or None


def _matrix_mask(config: dict[str, Any], selections: dict[str, str]) -> str:
    return resolve_mask(config.get("advanced", {}).get("matrix_mask", "representative"), selections)


def build_batches(config: dict, selections: dict[str, str]) -> list[dict]:
    """Return per-replica jobs; expensive pooled analyses are a separate API.

    RMSF is the mass-weighted residue average of atomic RMSFs in the common
    alignment frame. It is not the fluctuation of a residue's center of mass.
    SASA is an LCPO contribution to the configured solute's exposed area.
    """
    enabled = enabled_analyses(config)
    sections = [(str(item["name"]), resolve_mask(item["mask"], selections))
                for item in config.get("sections", [])]
    if len({name for name, _ in sections}) != len(sections):
        raise ValueError("Section names must be unique.")
    batches: list[dict] = []
    for item in config.get("solvation", []):
        if item["type"] not in enabled:
            continue
        kind, name = item["type"], item["name"]
        mask1, mask2 = resolve_mask(item["mask1"], selections), resolve_mask(item["mask2"], selections)
        if kind == "rdf":
            art = _artifact("rdf", name, "profile", "g(r)", x_column=0, x_unit="angstrom", normalization="average box volume", selections=[mask1, mask2])
            command = f"radial out {art['path']} {item['spacing']:g} {item['maximum']:g} {quote_mask(mask2)} {quote_mask(mask1)} volume A_RDF_{safe_name(name)}"
            description = "Pair RDF with volume normalization; user must keep maximum below half the shortest periodic cell height across frames."
        else:
            art = _artifact("watershell", name, "timeseries", "solvent count", x_column=0, selections=[mask1, mask2], cutoffs=[item["lower"], item["upper"]], definition="Cumulative inner/outer counts; select complete one-residue solvent molecules for CPU/CUDA consistency")
            noimage = " noimage" if config.get("imaging", {}).get("mode") == "off" else ""
            command = f"watershell {quote_mask(mask1)} out {art['path']} lower {item['lower']:g} upper {item['upper']:g} {quote_mask(mask2)}{noimage} A_WS_{safe_name(name)}"
            description = "Distance-based solvent shell counts using explicit solvent selection; not residence times."
        batches.append(_batch(kind + "_" + name, description, [command], [art], preserve_box=True))
    commands, artifacts = [], []
    for section, raw_mask in sections:
        token, mask = safe_name(section), quote_mask(raw_mask)
        if "rmsd" in enabled:
            for mode, flag in (("global", "nofit"), ("local", "nomod")):
                art = _artifact(f"rmsd_{mode}", section, "timeseries", "angstrom",
                                alignment=mode, reference="COMMON")
                commands.append(f"rms A_RMS_{mode}_{token} {mask} ref [COMMON] {flag} out {art['path']}")
                artifacts.append(art)
        if "rg" in enabled:
            art = _artifact("rg", section, "timeseries", "angstrom", mass_weighted=True)
            commands.append(f"radgyr A_RG_{token} {mask} mass nomax out {art['path']}")
            artifacts.append(art)
        if "rmsf" in enabled:
            art = _artifact("rmsf", section, "profile", "angstrom", x_label="Topology residue number",
                            definition="Mass-weighted residue average of atomic RMSFs; globally aligned")
            commands.append(f"atomicfluct A_RMSF_{token} {mask} byres out {art['path']}")
            artifacts.append(art)
    if commands:
        batches.append(_batch("structure", "Section RMSD, mass-weighted radius of gyration, and residue RMSF.", commands, artifacts))

    solute = _solute(selections)
    if "sasa" in enabled and sections:
        commands, artifacts = [], []
        # Explicit solutemask keeps a section's environment consistent across jobs.
        context = solute or "|".join(f"({mask})" for _, mask in sections)
        for section, raw_mask in sections:
            art = _artifact("sasa", section, "timeseries", "angstrom^2", method="LCPO",
                            context_mask=context, definition="Section contribution to exposed solute surface")
            commands.append(f"surf A_SASA_{safe_name(section)} {quote_mask(raw_mask)} solutemask {quote_mask(context)} out {art['path']}")
            artifacts.append(art)
        batches.append(_batch("sasa", "LCPO surface contributions; requires suitable topology atom types and bonds.", commands, artifacts))

    if "hbond" in enabled and solute:
        art = _artifact("hbond", "solute", "timeseries", "count")
        avg = _artifact("hbond_occupancy", "solute", "table", "fraction", format="hbond_occupancy")
        commands = [f"hbond A_HB {quote_mask(solute)} dist 3.0 angle 135 out {art['path']} avgout {avg['path']}"]
        artifacts = [art, avg]
        # Do not request uuresmatrix here: CPPTRAJ sizes it using the highest
        # non-solvent residue index in the topology, even outside this mask.
        # With interleaved solvent this can allocate an unexpectedly huge matrix.
        batches.append(_batch("hbond", "Geometric solute H-bonds (3.0 A, 135 degrees); inspect donor/acceptor chemistry.", commands, artifacts))

    if "dssp" in enabled and selections.get("protein"):
        art = _artifact("dssp", "protein", "timeseries", "state", categorical=True, palette="dssp")
        fraction = _artifact("dssp_fraction", "protein", "profile", "fraction")
        total = _artifact("dssp_content", "protein", "timeseries", "fraction")
        commands = [f"secstruct A_DSSP {quote_mask(selections['protein'])} out {art['path']} sumout {fraction['path']} totalout {total['path']}"]
        batches.append(_batch("dssp", "Protein secondary-structure assignment and per-residue fractions.", commands, [art, fraction, total]))

    if "nastruct" in enabled and selections.get("nucleic"):
        # NAStruct takes a residue RANGE, not an atom mask. Omitting resrange
        # searches recognized nucleic residues in the topology; sections apply
        # to the section observables above, not this whole-system chemistry.
        na = config.get("nucleic", {})
        command = "nastruct A_NA ref [COMMON] calcnohb sscalc groovecalc 3dna naout nastruct.dat noframespaces"
        if na.get("resrange"):
            residue_range = str(na["resrange"])
            if not re.fullmatch(r"[0-9,\-]+", residue_range):
                raise ValueError("nucleic.resrange must contain topology residue numbers, commas and hyphens.")
            command += f" resrange {residue_range}"
        for residue, base in na.get("resmap", {}).items():
            if (not re.fullmatch(r"[A-Za-z0-9_+\-]+", str(residue))
                    or not isinstance(base, str) or len(base) != 1 or base not in "ACGTU"):
                raise ValueError("nucleic.resmap requires residue names mapped to A, C, G, T or U.")
            command += f" resmap {residue}:{base}"
        artifacts = [_artifact(f"nastruct_{kind.lower()}", "nucleic", "table", "mixed",
                               path=f"{kind}.nastruct.dat", format="nastruct", optional=True,
                               reference="COMMON", parameter_group=kind, calcnohb=True,
                               note="Reference-defined pair geometry is retained after H-bond loss; assess BP and HB together")
                     for kind in ("BP", "BPstep", "Helix", "SS")]
        pucker = _artifact("pucker", "nucleic", "timeseries", "degree", circular=True, optional=True)
        artifacts.append(pucker)
        batches.append(_batch("nastruct", "Nucleic base pairs, helical parameters, grooves, strand geometry and sugar pucker.",
                              [command], artifacts,
                              [f"writedata {pucker['path']} A_NA[pucker]"]))

    for interface in config.get("interactions", []):
        name = str(interface["name"])
        token = safe_name(name)
        mask1 = resolve_mask(interface["mask1"], selections)
        mask2 = resolve_mask(interface["mask2"], selections)
        commands, artifacts = [], []
        if "distance" in enabled:
            art = _artifact("distance", name, "timeseries", "angstrom", definition="Center-of-mass distance, imaged when box present")
            commands.append(f"distance A_DIST_{token} {quote_mask(mask1)} {quote_mask(mask2)} out {art['path']}")
            artifacts.append(art)
        if "contacts" in enabled:
            cutoff = _positive(interface.get("contact_cutoff", 4.5), "contact_cutoff")
            heavy1, heavy2 = f"({mask1})&!@/H", f"({mask2})&!@/H"
            art = _artifact("contacts", name, "timeseries", "count", reference="COMMON", cutoff_angstrom=cutoff)
            avg = _artifact("contact_pairs", name, "table", "fraction", format="contact_pairs", reference="COMMON")
            maps = f"contacts_map__{token}.dat"
            command = (f"nativecontacts name A_CONTACT_{token} {quote_mask(heavy1)} {quote_mask(heavy2)} "
                       f"ref [COMMON] distance {cutoff:g} out {art['path']} writecontacts {avg['path']} byresidue map mapout {maps}")
            commands.append(command)
            artifacts.extend([art, avg])
            for category in ("native", "nonnative"):
                artifacts.append(_artifact(f"contacts_{category}_map", name, "matrix", "mean atom-pair contacts",
                                            path=f"{category}.{maps}", reference="COMMON",
                                            definition="Sum of atom-pair contact fractions within each residue pair; may exceed one"))
        if "hbond" in enabled:
            # Search both directions explicitly; a single donor/acceptor mask
            # pair would miss bonds donated by the opposite component.
            for direction, donor, acceptor in (("1to2", mask1, mask2), ("2to1", mask2, mask1)):
                art = _artifact(f"hbond_{direction}", name, "timeseries", "count")
                avg = _artifact(f"hbond_{direction}_occupancy", name, "table", "fraction", format="hbond_occupancy")
                # Restrict heavy donors/acceptors to F/O/N to preserve the
                # command's automatic chemical heuristic with explicit masks.
                commands.append(f"hbond A_HB_{token}_{direction} donormask {quote_mask(f'({donor})&@/N,O,F')} "
                                f"acceptormask {quote_mask(f'({acceptor})&@/N,O,F')} dist 3.0 angle 135 out {art['path']} avgout {avg['path']}")
                artifacts.extend([art, avg])
        if commands:
            extra = ({"contact_masks": [mask1, mask2], "resource_class": "contact_map"}
                     if "contacts" in enabled else {})
            batches.append(_batch(f"interface_{token}", f"Interactions between the two masks for {name}.", commands, artifacts, **extra))

    commands, artifacts = [], []
    for monitor in config.get("monitors", []):
        kind = monitor["type"]
        if kind not in enabled:
            continue
        name = monitor["name"]
        masks = [quote_mask(resolve_mask(m, selections)) for m in monitor["masks"]]
        center = monitor.get("center", "mass" if kind == "distance" else "geometry")
        flags = []
        if kind == "distance":
            if center == "geometry":
                flags.append("geom")
            if not monitor.get("image", True):
                flags.append("noimage")
        elif center == "mass":
            flags.append("mass")
        if kind == "dihedral" and monitor.get("range360", False):
            flags.append("range360")
        art = _artifact(f"monitor_{kind}", name, "timeseries",
                        "angstrom" if kind == "distance" else "degree",
                        circular=kind == "dihedral", center=center,
                        masks=monitor["masks"], definition=f"Ordered {kind} monitor: {name}")
        commands.append(f"{kind} A_MON_{safe_name(name)} {' '.join(masks)} {' '.join(flags)} out {art['path']}")
        artifacts.append(art)
    if commands:
        batches.append(_batch("geometry_monitors", "User-defined distances, angles and signed dihedral angles.", commands, artifacts))

    if "torsions" in enabled:
        for section, types in (("protein", "phi psi omega chip"), ("nucleic", "alpha beta gamma delta epsilon zeta chin")):
            if selections.get(section):
                art = _artifact("torsions", section, "timeseries", "degree", circular=True)
                command = f"multidihedral A_TORSION_{section} {types} resrange {quote_mask(selections[section])} out {art['path']}"
                batches.append(_batch(f"torsions_{section}", "Backbone and glycosidic/side-chain torsions, with circular statistics required.", [command], [art]))

    if "dccm" in enabled:
        mask = _matrix_mask(config, selections)
        art = _artifact("dccm", "representative", "matrix", "correlation", coordinate_mask=mask,
                        definition="Correlation of representative-atom displacement vectors after global alignment",
                        limits=[-1, 1])
        batches.append(_batch("dccm", "Dynamic cross-correlation of the selected representative atoms.",
                              [f"matrix correl {quote_mask(mask)} name A_DCCM out {art['path']}"], [art],
                              coordinate_mask=mask, resource_class="atom_matrix"))
    if "stacking" in enabled:
        from .stacking import build_stacking_batches
        batches.extend(build_stacking_batches(config))
    return batches


def build_pooled_batches(config: dict, selections: dict[str, str]) -> list[dict]:
    """Build pooled PCA/clustering/pairwise jobs with a common coordinate basis.

    The runner must preserve replica boundaries, enforce its frame/atom caps,
    and split timeseries back into replicas after execution. Pooled covariance
    is frame-weighted: longer replicas contribute more to the estimated basis.
    Cluster representatives contain only the selected atoms (RAM reduction).
    """
    enabled = enabled_analyses(config)
    if not enabled & {"pca", "cluster", "pairwise_rmsd"}:
        return []
    options = config.get("advanced", {})
    mask = _matrix_mask(config, selections)
    strip = f"strip {quote_mask(f'!({mask})')}"
    commands = [strip, "createcrd A_POOLED"]
    batches = []
    if "pca" in enabled:
        nmodes = _positive(options.get("pca_modes", 3), "pca_modes", integer=True)
        projection = _artifact("pca", "representative", "timeseries", "angstrom", pooled=True,
                               coordinate_mask=mask, basis="pooled_frame_weighted_common_alignment")
        covariance = _artifact("pca_covariance", "representative", "matrix", "angstrom^2", pooled=True,
                               coordinate_mask=mask, axis="selected Cartesian coordinate")
        modes = _artifact("pca_modes", "representative", "table", "mixed", pooled=True,
                          format="eigenvectors", note="Requested modes only; do not normalize these to claim total explained variance")
        post = ["crdaction A_POOLED matrix covar @* name A_COVAR",
                f"writedata {covariance['path']} A_COVAR",
                f"runanalysis diagmatrix A_COVAR out {modes['path']} vecs {nmodes} name A_MODES",
                f"crdaction A_POOLED projection A_PC modes A_MODES beg 1 end {nmodes} @*",
                f"writedata {projection['path']} A_PC"]
        batches.append(_batch("pca_pooled", "One pooled PCA basis for all replicas, using globally aligned coordinates.",
                              commands[:], [projection, covariance, modes], post,
                              coordinate_mask=mask, resource_class="pooled_pca", pooled=True))
    if "cluster" in enabled:
        count = _positive(options.get("cluster_count", 5), "cluster_count", integer=True)
        sieve = _positive(options.get("cluster_sieve", 10), "cluster_sieve", integer=True)
        assignment = _artifact("cluster", "representative", "timeseries", "cluster", pooled=True, categorical=True)
        summary = _artifact("cluster_population", "representative", "table", "mixed", pooled=True, format="cluster_summary")
        metric = options.get("cluster_metric", "rms")
        algorithm = options.get("cluster_algorithm", "kmeans")
        algorithm_args = (f"kmeans clusters {count} randompoint kseed 2026 maxit 200"
                          if algorithm == "kmeans" else f"hieragglo clusters {count} averagelinkage")
        cluster_commands = commands[:]
        metric_args = f"{metric} @*"
        if metric == "features":
            # Measure on original topology before stripping: masks retain atom IDs.
            geometry = next(b for b in build_batches(config, selections) if b["name"] == "geometry_monitors")
            names = [f["monitor"] for f in options["cluster_features"]]
            feature_commands = [c for c in geometry["commands"] if c.split()[1] in {f"A_MON_{n}" for n in names}]
            cluster_commands = feature_commands + cluster_commands
            datasets = ','.join(f"A_MON_{n}" for n in names)
            weights = ','.join(str(f["weight"]) for f in options["cluster_features"])
            metric_args = f"data {datasets} euclid wgt {weights}"
        command = (f"runanalysis cluster A_CLUSTER crdset A_POOLED {algorithm_args} "
                   f"{metric_args} sieve {sieve} random sieveseed 2026 "
                   f"out {assignment['path']} summary {summary['path']} "
                   "info cluster_info.dat singlerepout cluster_representatives.pdb singlerepfmt pdb")
        batches.append(_batch("cluster_pooled", f"Pooled {algorithm} clustering using {metric}; shared labels across replicas.",
                              cluster_commands, [assignment, summary], [command],
                              coordinate_mask=mask, resource_class="pooled_cluster", pooled=True))
    if "pairwise_rmsd" in enabled:
        art = _artifact("pairwise_rmsd", "representative", "matrix", "angstrom", pooled=True,
                        axis="Pooled sampled frame", coordinate_mask=mask)
        batches.append(_batch("pairwise_pooled", "Optional quadratic frame-by-frame RMSD over the pooled selected coordinates.",
                              commands[:], [art], [f"runanalysis rms2d A_PAIRWISE crdset A_POOLED @* out {art['path']}"],
                              coordinate_mask=mask, resource_class="pooled_pairwise", pooled=True))
    return batches
