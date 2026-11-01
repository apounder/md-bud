"""Interactive wizard and reproducible command-line entry points."""
from __future__ import annotations

import argparse
import copy
import json
import math
import tempfile
import shutil
import re
import shlex
import sys
from pathlib import Path

from . import __version__
from .analysis import KNOWN_ANALYSES, LEVELS, enabled_analyses
from .config import DEFAULTS, discover, load_config, natural_key, normalize_config, safe_name, save_config
from .topology import inspect_topology
from .selections import normalize_mask


def ask(prompt, default=None, convert=str, valid=None, *, help_text=None):
    from .wizard_help import explanation
    print("  " + (help_text or explanation(prompt)))
    while True:
        suffix = f" [{default}]" if default is not None else ""
        raw = input(prompt + suffix + ": ").strip()
        raw = str(default) if not raw and default is not None else raw
        try:
            value = convert(raw)
            if valid and not valid(value):
                raise ValueError("outside the allowed range")
            return value
        except (ValueError, TypeError) as exc:
            print(f"Please enter a valid value ({exc}).")


def yes(prompt, default=True):
    return ask(prompt + " (y/n)", "y" if default else "n", str.lower, lambda x: x in {"y", "n"}) == "y"


def _indices(raw, maximum):
    if raw.lower() == "all":
        return list(range(maximum))
    indices = []
    for part in raw.replace(" ", "").split(","):
        m = re.fullmatch(r"(\d+)(?:-(\d+))?", part)
        if not m:
            raise ValueError("use all, 1,3,5 or 1-4")
        start, stop = int(m[1]), int(m[2] or m[1])
        if start < 1 or stop > maximum or stop < start:
            raise ValueError("selection index outside the displayed list")
        indices.extend(range(start - 1, stop))
    if len(set(indices)) != len(indices):
        raise ValueError("do not repeat an index")
    return indices


def numbered(paths):
    for i, path in enumerate(paths, 1):
        print(f"  {i:3}. {path}")


def _group_files(paths, mode):
    if mode == "files":
        return [{"name": f"rep{i:02}", "trajectories": [str(path)]} for i, path in enumerate(paths, 1)]
    if mode == "single":
        return [{"name": "rep01", "trajectories": [str(path) for path in paths]}]
    groups = {}
    for path in paths:
        groups.setdefault(path.parent, []).append(path)
    return [{"name": f"rep{i:02}", "trajectories": [str(path) for path in chunks]}
            for i, (_, chunks) in enumerate(groups.items(), 1)]


def prompt_mask(prompt, default=None):
    return ask(prompt, default, normalize_mask)


def _quick_geometry(raw, count):
    from .config import safe_text
    values = [v.strip() for v in raw.split(";")]
    if len(values) != count:
        raise ValueError(f"Enter exactly {count} ordered selections separated by semicolons.")
    result = []
    for value in values:
        match = re.fullmatch(r"([1-9]\d*):([A-Za-z0-9'+_-]+)", value)
        value = f":{match[1]}@{match[2]}" if match else normalize_mask(value)
        result.append(safe_text(value, "geometry selection"))
    return result


def _geometry_atom(info, position):
    residues = {r["index"]: r for r in (info or {}).get("residues", [])}
    residue = ask(f"Point {position}: residue number", None, int,
                  lambda n: n > 0 and (not residues or n in residues))
    atoms = residues[residue]["atoms"] if residue in residues else []
    if atoms:
        print(f"  {residues[residue]['name']} {residue}: " + ", ".join(dict.fromkeys(atoms)))
    atom = ask(f"Point {position}: atom name", None, str,
               lambda name: bool(re.fullmatch(r"[A-Za-z0-9'+_-]+", name)) and (not atoms or name in atoms))
    return f":{residue}@{atom}"


def configure_monitors(config, kind, info=None):
    count = {"distance": 2, "angle": 3, "dihedral": 4}[kind]
    print(f"\n{kind.capitalize()} monitoring: provide {count} ordered atom/group masks.")
    print("Single atoms: :25@CA or @120. Groups: backbone, :1-20, or a named section.")
    if kind == "angle":
        print("Mask 2 is the vertex: mask1--mask2--mask3 (0 to 180 degrees).")
    elif kind == "dihedral":
        print("Masks 2--3 define the central axis; mask order controls the signed torsion.")
    if kind != "distance":
        print("Angles use the prepared coordinates; keep all groups in a consistent periodic image.")
    while True:
        name = ask("Monitor name", f"{kind}_{len(config['monitors']) + 1}", safe_name)
        if name.casefold() in {m["name"].casefold() for m in config["monitors"]}:
            print("That monitor name is already used.")
            continue
        print("  1. Guided residue + atom names\n  2. Named regions / groups\n  3. Quick ordered selections (e.g. 25:CA; 41:CA)")
        mode = ask("Geometry entry method", 1, int, lambda n: n in {1, 2, 3})
        if mode == 1:
            masks = [_geometry_atom(info, i) for i in range(1, count + 1)]
        elif mode == 2:
            groups = list(dict.fromkeys([s["name"] for s in config["sections"]] + list((info or {}).get("selections", {}))))
            numbered(groups)
            masks = []
            for i in range(1, count + 1):
                index = ask(f"Point {i}: choose group number (0 = explicit mask)", 0, int, lambda n: 0 <= n <= len(groups))
                masks.append(groups[index - 1] if index else prompt_mask(f"Point {i}: residue/group mask"))
        else:
            masks = ask("Ordered selections, separated by semicolons", None, lambda raw: _quick_geometry(raw, count))
        monitor = {"name": name, "type": kind, "masks": masks}
        monitor["center"] = "geometry" if mode == 1 else ask("For groups use mass or geometry centers", "mass" if kind == "distance" else "geometry", str.lower, lambda s: s in {"mass", "geometry"})
        if kind == "distance":
            monitor["image"] = yes("Use minimum-image distance when a periodic box exists?", config["imaging"]["mode"] != "off")
        if kind == "dihedral":
            monitor["range360"] = yes("Write raw angles in 0..360 instead of -180..180 degrees?", False)
        print(f"  Review {kind}: " + " → ".join(masks))
        if not yes("Keep this monitor?", True):
            continue
        config["monitors"].append(monitor)
        if not yes(f"Add another {kind} monitor?", False):
            break


def additional_analyses(config, info):
    menu = [
        ("distance", "Distance monitors: two atoms or group centers"),
        ("angle", "Angle monitors: three atoms/groups; second is the vertex"),
        ("dihedral", "Dihedral monitors: four ordered atoms/groups"),
        ("sasa", "SASA: LCPO surface area for your sections"),
        ("nastruct", "NAStruct: nucleic pair, step, helix, groove and pucker parameters"),
        ("hbond", "Hydrogen bonds within solute and configured interfaces"),
        ("contacts", "Interface native/nonnative contacts and residue maps"),
        ("dssp", "Protein secondary structure"),
        ("torsions", "Automatic protein/nucleic backbone and chi torsions"),
        ("dccm", "Dynamic cross-correlation matrix"),
        ("pca", "Shared pooled principal component analysis"),
        ("cluster", "Shared pooled RMSD clustering"),
        ("pairwise_rmsd", "Pooled frame-by-frame RMSD matrix; quadratic cost"),
        ("rmsd", "Global/local section RMSD"),
        ("rg", "Section radius of gyration"),
        ("rmsf", "Residue RMSF"),
        ("rdf", "Radial distribution: pair structure versus distance (periodic systems)"),
        ("watershell", "Solvent shells: local solvent counts around a selected region"),
        ("stacking", "Pi-stacking: discover and name rings from residue numbers"),
    ]
    if not yes("Would you like any additional analyses or custom monitors?", False):
        return
    preset = LEVELS[config["analysis"]["level"]]
    for i, (key, description) in enumerate(menu, 1):
        tag = " (in preset; configure monitors)" if key == "distance" and key in preset else " (already in preset)" if key in preset else ""
        print(f"  {i:2}. {description}{tag}")
    chosen = ask("Enter numbers, e.g. 1,2,5 or 1-3; 0 cancels", "0",
                 lambda s: [] if s == "0" else _indices(s, len(menu)))
    selected = [menu[i][0] for i in chosen]
    config["analysis"]["additional"] = sorted(set(selected) - set(preset))
    if "stacking" in selected:
        configure_stacking(config)
    for kind in ("distance", "angle", "dihedral"):
        if kind in selected:
            configure_monitors(config, kind, info)
    for kind in ("rdf", "watershell"):
        if kind in selected:
            if kind == "watershell":
                print("  Select complete solvent molecules with one residue per molecule, such as :WAT, for consistent CPU/CUDA behavior.")
            while True:
                item = {"type": kind, "name": ask("Solvation analysis name", f"{kind}_{len(config.get('solvation', [])) + 1}", safe_name),
                        "mask1": prompt_mask("Reference solute/region mask"),
                        "mask2": prompt_mask("Surrounding solvent/particle mask")}
                if kind == "rdf":
                    item["spacing"] = ask("RDF bin spacing in angstroms", 0.1, float, lambda n: n > 0,
                                          help_text="Smaller bins resolve finer structure but require more sampling per bin.")
                    item["maximum"] = ask("RDF maximum radius in angstroms", 10.0, float, lambda n: n > item["spacing"],
                                          help_text="Keep this below half the shortest periodic cell height over all analyzed frames; this geometric limit is your responsibility to verify.")
                    print("  RDF uses average-volume normalization and requires a valid periodic box; a heterogeneous solute environment need not approach g(r)=1.")
                else:
                    item["lower"] = ask("Inner solvent-shell cutoff in angstroms", 3.4, float, lambda n: n > 0,
                                        help_text="This distance defines the inner solvent shell around the selected solute atoms.")
                    item["upper"] = ask("Outer solvent-shell cutoff in angstroms", 5.0, float, lambda n: n > item["lower"],
                                        help_text="The outer count includes the inner count; subtract them for the intervening shell population.")
                config.setdefault("solvation", []).append(item)
                if not yes("Add another solvation selection?", False):
                    break
    if "sasa" in selected:
        print("SASA will use the sections you chose above, with the detected solute as surface context.")
    if "nastruct" in selected and "nucleic" not in info["selections"]:
        print("No nucleic residues were detected. NAStruct needs a section named nucleic with the appropriate mask.")
        if not any(s["name"] == "nucleic" for s in config["sections"]):
            config["sections"].append({"name": "nucleic", "mask": prompt_mask("Nucleic residue mask")})
    if "dssp" in selected and "protein" not in info["selections"]:
        if not any(s["name"] == "protein" for s in config["sections"]):
            config["sections"].append({"name": "protein", "mask": prompt_mask("Protein residue mask for DSSP")})
    if "contacts" in selected and not config["interactions"] and not {"protein", "nucleic"} <= set(info["selections"]):
        print("Contact analysis needs two regions. Define an interface for this calculation.")
        config["interactions"].append({"name": ask("Interface name", "interface1", safe_name),
                                       "mask1": prompt_mask("First interface mask"),
                                       "mask2": prompt_mask("Second interface mask"),
                                       "contact_cutoff": ask("Contact cutoff in angstroms", 4.5, float, lambda x: x > 0)})


def configure_stacking(config):
    from .rings import discover_pi_systems, read_bond_graph
    from .stacking import stacking_pairs
    atoms, _, _ = read_bond_graph(config["topology"])
    available = {a["residue"] for a in atoms}
    print("\nPi-system discovery uses topology bonds, not atom-name templates. Residues use topology numbering.")
    def residue_input(text):
        values = {i+1 for i in _indices(text.lstrip(':'), max(available))}
        if not values <= available:
            raise ValueError("Choose residues present in the topology")
        return values
    targets = ask("Target residue numbers (e.g. 24,30-32)", None, residue_input)
    candidates, notes = discover_pi_systems(config["topology"])
    target_candidates = [r for r in candidates if r['residue'] in targets]
    for note in notes:
        if any(note.startswith(f"Residue {r} ") for r in targets):
            print("  " + note)
    if not target_candidates:
        raise ValueError("No candidate pi systems in the selected residues. Check bonded topology/atom types; saturated rings are excluded.")
    def confirm(candidates, title):
        print("\n" + title)
        for i, ring in enumerate(candidates, 1):
            print(f"  {i}. Residue {ring['residue']} {ring['residue_name']}: {ring['kind']}; "
                  f"atoms {','.join(ring['atom_names'])}; topology indices {ring['atoms']}; {ring['evidence']}")
        print("Connectivity alone does not prove aromaticity. Confirm the chemistry of uncertain candidates; guanidinium groups are not rings.")
        # Prefer whole fused systems; users may select constituent rings instead.
        default = ','.join(str(i+1) for i,r in enumerate(candidates)
                           if not any(set(r['atoms']) < set(other['atoms']) for other in candidates))
        indices = ask("Confirm pi systems of interest (numbers/all; 0 = none)", default,
                      lambda value: [] if value == '0' else _indices(value,len(candidates)))
        return [dict(candidates[i]) for i in indices]
    selected_targets = confirm(target_candidates, "Candidate target rings / pi systems")
    if not selected_targets:
        raise ValueError("No target pi systems confirmed")
    named = set()
    for ring in selected_targets:
        ring['name'] = ask(f"Name the pi system in residue {ring['residue']}", ring['name'], safe_name,
                           lambda name: name.casefold() not in named)
        named.add(ring['name'].casefold())
    config["stacking"]["exclude_same_residue"] = not yes("Include contacts between distinct rings in the same residue?", False)
    partner_residues = ask("Partner residue scope: all, or residue numbers/ranges", "all",
                           lambda value: available if value == 'all' else residue_input(value))
    target_atom_sets = {tuple(r['atoms']) for r in selected_targets}
    partners = [r for r in candidates if r['residue'] in partner_residues and tuple(r['atoms']) not in target_atom_sets
                and (not config['stacking']['exclude_same_residue'] or r['residue'] not in targets)]
    selected_partners = confirm(partners, "Candidate partner pi systems (tested over every selected frame)") if partners else []
    for ring in selected_partners:
        # Automatic names are shown above; custom naming is optional for many partners.
        while ring['name'].casefold() in named:
            ring['name'] += '_partner'
        named.add(ring['name'].casefold())
    if selected_partners and yes("Customize partner ring names?", False):
        for ring in selected_partners:
            named.remove(ring['name'].casefold())
            ring['name'] = ask(f"Partner name for residue {ring['residue']}", ring['name'], safe_name,
                               lambda name: name.casefold() not in named)
            named.add(ring['name'].casefold())
    options = config['stacking']
    options['rings'] = selected_targets + selected_partners
    options['targets'] = [r['name'] for r in selected_targets]
    options['distance_cutoff'] = ask("Stacking centroid-distance cutoff (angstrom)", 5.0, float,
                                     lambda x: math.isfinite(x) and 0 < x <= 100)
    options['angle_cutoff'] = ask("Parallel/antiparallel plane-angle tolerance (degrees)", 30.0, float,
                                  lambda x: math.isfinite(x) and 0 < x <= 90)
    options['center'] = ask("Ring center: geometry/mass (donor script used mass)", "geometry", str.lower,
                            lambda x: x in {'geometry','mass'})
    pairs = stacking_pairs(options)
    if not pairs:
        raise ValueError("No nonoverlapping inter-residue pi-system pairs were confirmed")
    print(f"Stacking: {len(pairs)} unique pairs; per-target and per-residue any-partner occupancies, plus equal-replica mean and SD.")


def configure_frame_range(config):
    """Length inspection and trimming precede all analysis-selection questions."""
    from .runner import executable, read_replica_lengths, show_replica_lengths, length_warning, RunError
    replicas = None
    try:
        binary = executable(config)
    except RunError:
        print("Trajectory lengths are unverified: CPPTRAJ is unavailable. The run will check counts before analysis.")
    else:
        print("Loading trajectory length metadata with CPPTRAJ...")
        directory = tempfile.mkdtemp(prefix="mdwb-lengths-")
        try:
            replicas = read_replica_lengths(binary, config, directory)
        except BaseException:
            print(f"Trajectory length-check logs retained in {directory}")
            raise
        else:
            shutil.rmtree(directory)
        show_replica_lengths(replicas)
    print("Choose the part to analyze now; settings apply across joined segments in every replica.")
    dt = ask("Saved-frame interval in ps (0 = unknown; use frame axes)", 0.0, float, lambda x: math.isfinite(x) and x >= 0)
    config["frames"]["dt_ps"] = dt or None
    if replicas and dt:
        show_replica_lengths(replicas, dt)
    shortest = min(sum(r["lengths"]) for r in replicas) if replicas else None
    endpoint = -1
    if replicas and length_warning(replicas) and yes("Use the shortest replica's endpoint for all replicas?", False):
        endpoint = shortest
    unit = ask("Select analysis range in frames/ns", "frames", str.lower, lambda x: x in {"frames", "ns"}) if dt else "frames"
    if unit == "ns":
        interval = dt / 1000
        print("Times are relative to saved frame 1 = 0 ns. Start/end bounds are inclusive; only saved frames inside them are kept.")
        start_ns = ask("Start time in ns (discard the beginning)", 0.0, float,
                       lambda x: math.isfinite(x) and x >= 0 and (shortest is None or x <= (shortest-1)*interval))
        start = math.ceil(start_ns/interval - 1e-10) + 1
        end_ns = ask("End time in ns (-1 = each replica's end)", -1 if endpoint == -1 else (endpoint-1)*interval, float,
                     lambda x: math.isfinite(x) and (x == -1 or (x >= start_ns and math.floor(x/interval + 1e-10)+1 >= start and (shortest is None or math.floor(x/interval + 1e-10)+1 <= shortest))))
        stop = -1 if end_ns == -1 else math.floor(end_ns/interval + 1e-10) + 1
    else:
        start = ask("First frame to analyze (discard equilibration before this)", 1, int,
                    lambda x: x > 0 and (shortest is None or x <= shortest))
        stop = ask("Last frame (-1 = each replica's end)", endpoint, int, lambda x: x == -1 or (x >= start and (shortest is None or x <= shortest)))
    stride = ask("Global analysis stride (1 = analyze every saved frame)", 1, int, lambda x: x > 0,
                 help_text="Applies to every analysis. Keep 1 for full sampling; clustering has its own separate sieve setting.")
    config["frames"].update(start=start, stop=stop, stride=stride)
    if replicas:
        from .runner import select_windows
        for replica in replicas:
            count = sum(w["count"] for w in select_windows(replica["lengths"], config["frames"]) if w)
            print(f"  Analyze {replica['name']}: {count} selected frames")
    return replicas


def configure_resource_limits(config, replicas=None):
    enabled = enabled_analyses(config)
    if not enabled & {"dccm", "pca", "cluster", "pairwise_rmsd"}:
        return
    options = config["advanced"]
    print("\nAdvanced analysis resource limits (guards, not memory allocations).")
    if enabled & {"pca", "cluster", "pairwise_rmsd"}:
        from .runner import select_windows

        def pooled_count():
            frames = config["frames"]
            return sum(w["count"] for r in replicas
                       for w in select_windows(r["lengths"], frames) if w)

        total = pooled_count() if replicas else None
        print("The frame limit counts selected frames across ALL replicas, after trimming and stride; cluster sieve does not reduce this count.")
        if total is not None:
            print(f"Pooled selection: {total:,} frames across {len(replicas)} replicas; current limit {options['max_frames']:,}.")
            if total > options["max_frames"]:
                print("The suggested cap below covers all selected frames. Analysis sampling stays unchanged.")
        else:
            print("Pooled frame count is unverified; the run will check it against this limit.")
        options["max_frames"] = ask("Maximum pooled frames", max(options["max_frames"], total or 1), int,
            lambda x: x > 0 and (total is None or x >= total),
            help_text="Must cover the selected total when known. Raising this cap still leaves the memory guard in force.")
    if "cluster" in enabled:
        print(f"Clustering alone uses sieve {options['cluster_sieve']}: fit a random subset, then assign remaining frames. Other analyses keep their selected frames.")
    options["max_matrix_atoms"] = ask("Maximum selected atoms for advanced analyses", options["max_matrix_atoms"], int,
        lambda x: x > 0, help_text="The run checks the advanced coordinate mask against this atom limit.")
    options["max_memory_gb"] = ask("Maximum estimated memory in GiB", options["max_memory_gb"], float,
        lambda x: math.isfinite(x) and 0 < x < 100000,
        help_text="Choose a budget within the RAM available to your job. The run checks a conservative estimate; this does not allocate RAM.")


def configure_convergence(config):
    options = config["convergence"]
    print("\nConvergence over time: cumulative/window means, distributions, ESS and replica agreement.")
    options["enabled"] = yes("Analyze convergence within and between replicas?", True)
    if not options["enabled"]:
        return
    if yes("Set a gold-standard coordinate structure for RMSD?", False):
        config["reference"] = str(Path(ask("Reference structure path")).expanduser().resolve())
    print("Pooled = equal-replica mean and total fluctuation SD. Tail = each replica's last X ns.")
    print("Custom/file targets use exact analysis/section/series IDs and units from reports/all_replicas.dat.")
    options["reference_mode"] = ask("Statistical reference: pooled/tail/custom/file", "pooled", str.lower,
                                    lambda x: x in {"pooled", "tail", "custom", "file"})
    if options["reference_mode"] == "tail":
        if not config["frames"]["dt_ps"]:
            config["frames"]["dt_ps"] = ask("Saved-frame interval in ps (required for ns)", None, float, lambda x: x > 0)
        options["tail_ns"] = ask("Last X ns per replica for reference average and SD", 10.0, float, lambda x: x > 0)
    elif options["reference_mode"] == "file":
        options["reference_path"] = str(Path(ask("Reference statistics JSON path")).expanduser().resolve())
    elif options["reference_mode"] == "custom":
        while True:
            key = ask("Metric ID: analysis/section/series", None, str, lambda x: len(x.rsplit('/', 2)) == 3 and all(x.split('/')))
            options["references"][key] = {"mean": ask("Reference mean", None, float),
                                           "sd": ask("Reference standard deviation", None, float, lambda x: x >= 0),
                                           "unit": ask("Exact metric unit (e.g. angstrom)")}
            if not yes("Add another reference metric?", False):
                break
    metrics = ask("Convergence metrics: all, or comma-separated analysis names / exact metric IDs", "all")
    options["metrics"] = [] if metrics == "all" else [m.strip() for m in metrics.split(',') if m.strip()]
    options["window_frames"] = ask("Trailing window in analyzed frames", 50, int, lambda x: 2 <= x <= 10000)
    options["checkpoints"] = ask("Number of convergence checkpoints", 20, int, lambda x: 2 <= x <= 10000)


def configure_clustering(config, info):
    if "cluster" not in enabled_analyses(config):
        return
    options = config["advanced"]
    print("\nShared clustering: RMSD, distance-RMSD (DME), symmetry RMSD, or geometry features.")
    print("Geometry features can describe bonds, metal coordination, donor/acceptor distances, angles and torsions for QM/MM.")
    options["cluster_metric"] = ask("Cluster metric: rms/dme/srmsd/features", "rms", str.lower,
                                     lambda x: x in {"rms", "dme", "srmsd", "features"})
    if options["cluster_metric"] == "features":
        while yes("Add a bond/distance, angle or dihedral feature monitor?", not bool(config["monitors"])):
            kind = ask("Geometry: bond/distance/angle/dihedral", "distance", str.lower,
                       lambda x: x in {"bond", "distance", "angle", "dihedral"})
            kind = "distance" if kind == "bond" else kind
            configure_monitors(config, kind, info)
            config["analysis"]["additional"] = sorted(set(config["analysis"]["additional"]) | {kind})
        names = [m["name"] for m in config["monitors"]]
        numbered(names)
        chosen = ask("Choose feature monitors", "all", lambda x: _indices(x, len(names)))
        print("Choose weights explicitly: raw angstrom and degree scales differ. These are CPPTRAJ metric weights.")
        options["cluster_features"] = [{"monitor": names[i], "weight": ask(f"Weight for {names[i]}", 1.0, float, lambda x: x > 0)} for i in chosen]
    options["matrix_mask"] = prompt_mask("Cluster representative/coordinate mask (e.g. active-site residues)", options["matrix_mask"])
    options["cluster_algorithm"] = ask("Cluster algorithm: kmeans/hieragglo", "kmeans", str.lower,
                                        lambda x: x in {"kmeans", "hieragglo"})
    options["cluster_count"] = ask("Cluster count (exploratory; tune for your system)", 5, int, lambda x: x > 0)
    options["cluster_sieve"] = ask("Clustering sieve (fit a random 1/N subset, assign remaining frames)", 10, int, lambda x: x > 0)


def wizard(root, config_path, quick=None):
    root = Path(root).expanduser().resolve()
    target = Path(config_path).expanduser().resolve()
    existing = [target] if target.is_file() else []
    if not existing and str(config_path) == "analysis_config.json":
        # Only inspect nearby configuration files, not JSON in result trees.
        for path in sorted(set(Path.cwd().glob("*.json")) | set(root.glob("*.json")), key=natural_key):
            try:
                data = json.loads(path.read_text(encoding="utf-8-sig"))
            except (OSError, ValueError):
                continue
            if isinstance(data, dict) and "topology" in data and "replicas" in data:
                existing.append(path)
    if existing:
        if len(existing) > 1:
            print("\nSaved analysis configurations:")
            numbered(existing)
            index = ask("Choose saved configuration", 1, int, lambda n: 1 <= n <= len(existing))
            target = existing[index - 1]
        else:
            target = existing[0]
        print(f"\nExisting configuration found: {target}")
        print("To create another configuration, use wizard --config with a new filename.")
        return offer_run(target)
    from .console import panel
    panel("MD-BUD · SETUP", [("Basic path", "Detected biomolecular regions, RMSD/RMSF/Rg, optional extras"),
          ("Guided path", "Custom regions, fitting, frame selection and specialist parameters"),
          ("Data safety", "Inputs are read only; calculations write to a separate results folder")])
    if quick is None:
        quick = ask("Setup path: basic/guided", "basic", str.lower, lambda value: value in {"basic", "guided"},
                    help_text="Basic keeps common defaults and offers extras; guided exposes all selection and analysis choices.") == "basic"
    print("\nMD-bud — interactive setup")
    print("One configuration compares replicas of the SAME system with one matching atom order/topology.")
    print("Files in one replica are sequential chunks; independent replicas remain separate.")
    candidates = discover(root)
    print(f"\nSearching {root}")
    print("Topology candidates:")
    numbered(candidates["topologies"])
    if candidates["topologies"]:
        index = ask("Choose topology number (0 for an explicit path)", 1, int, lambda x: 0 <= x <= len(candidates["topologies"]))
        topology = candidates["topologies"][index - 1] if index else Path(ask("Topology path")).expanduser().resolve()
    else:
        topology = Path(ask("No topology discovered; enter its path")).expanduser().resolve()
    if not topology.is_file():
        raise ValueError(f"Topology not found: {topology}")
    info = inspect_topology(topology)
    print(f"\nTopology inspection: {info['natom'] or 'unknown'} atoms; {len(info['residues']) or 'unknown'} residues")
    for warning in info["warnings"]:
        print("  Review: " + warning)
    print("\nTrajectory candidates (natural filename order; review chronology yourself):")
    numbered(candidates["trajectories"])
    if candidates["trajectories"]:
        selected = ask("Choose trajectories; e.g. all or 1,3-6 (order retained)", "all", lambda s: _indices(s, len(candidates["trajectories"])))
        paths = [candidates["trajectories"][i] for i in selected]
    else:
        print("Enter one trajectory path at a time; blank finishes. Other CPPTRAJ-supported formats can be entered here.")
        paths = []
        while True:
            value = input("Trajectory path: ").strip()
            if not value:
                break
            paths.append(Path(value).expanduser().resolve())
        if not paths:
            raise ValueError("No trajectories selected.")
    print("\nGrouping choices: files = one replica/file; folders = one replica/parent folder; single = all chunks of one run; custom = explicit groups.")
    recommended = "folders" if len({p.parent for p in paths}) > 1 else ("single" if len(paths) == 1 else "custom")
    mode = ask("How should files be grouped?", recommended, str.lower, lambda s: s in {"files", "folders", "single", "custom"})
    if mode == "custom":
        numbered(paths)
        replicas, used = [], set()
        print("Assign each selected file exactly once. Indices within a group set its chronological order.")
        while len(used) < len(paths):
            name = ask("Replica name", f"rep{len(replicas) + 1:02}", safe_name)
            indices = ask("Trajectory indices for this replica", None, lambda s: _indices(s, len(paths)))
            if set(indices) & used or name in {r["name"] for r in replicas}:
                print("That name or file is already assigned; try again.")
                continue
            used.update(indices)
            replicas.append({"name": name, "trajectories": [str(paths[i]) for i in indices]})
    else:
        replicas = _group_files(paths, mode)
    print("\nReview the proposed grouping and file order:")
    for replica in replicas:
        print(f"  {replica['name']}")
        for path in replica["trajectories"]:
            print(f"    {path}")
    if not yes("Is this grouping and chronological order correct?"):
        raise ValueError("Setup stopped before saving. Restart and choose custom grouping to change the order.")
    config = copy.deepcopy(DEFAULTS)
    config["study"]["title"] = ask("Study title for reports", root.name or "MD analysis")
    config.update(topology=str(topology), replicas=replicas)
    length_replicas = configure_frame_range(config)
    if quick:
        names = [name for name in ("complex", "backbone", "nucleic_backbone") if name in info["selections"]]
        if not names or "representative" not in info["selections"]:
            raise ValueError("Automatic regions are unavailable; restart with the guided setup path and explicit masks.")
        config["sections"] = [{"name": name, "mask": name} for name in names]
        config["analysis"]["level"] = "basic"
        config["imaging"]["mode"] = "auto" if info["periodic"] is not None else ask("Imaging: on/off", "off", str.lower, lambda s: s in {"on", "off"}, help_text="Use on only when the trajectory contains a valid periodic box.")
        panel("BASIC ANALYSIS PLAN", [("Sections", ", ".join(names)), ("Analyses", "RMSD, RMSF, radius of gyration"),
              ("Alignment", "Detected representative atoms; first selected frame of replica 1"), ("Imaging", config["imaging"]["mode"])])
        additional_analyses(config, info)
        configure_convergence(config)
        configure_clustering(config, info)
        configure_resource_limits(config, length_replicas)
        return finish_wizard(config, config_path, root)
    print("\nAvailable selections use CPPTRAJ topology residue numbers (starting at 1):")
    for name, mask in info["selections"].items():
        print(f"  {name:20} {mask}")
    print("Examples: :1-120, :5,8,12-20, :1-100@CA, (:1-100)&!@/H")
    print("Bare residue ranges such as 10-50 are accepted here and converted to :10-50.")
    suggested = [n for n in ("complex", "backbone", "nucleic_backbone") if n in info["selections"]]
    config["sections"] = [{"name": name, "mask": name} for name in suggested]
    print("Suggested sections: " + (", ".join(suggested) or "none; add at least one"))
    if suggested and not yes("Keep the suggested sections?"):
        config["sections"] = []
    while yes("Add a custom section?", not bool(config["sections"])):
        name = ask("Section name", None, safe_name)
        mask = prompt_mask("Named selection or CPPTRAJ mask")
        if name in {s["name"] for s in config["sections"]}:
            print("A section with this name already exists.")
            continue
        config["sections"].append({"name": name, "mask": mask})
    if not config["sections"]:
        raise ValueError("Select at least one section.")
    default_fit = "representative" if "representative" in info["selections"] else config["sections"][0]["mask"]
    print("\nThe global fit removes translation/rotation. For relative domain motion, choose the stable reference domain/core.")
    config["fit_mask"] = prompt_mask("Global alignment mask", default_fit)
    reference = ask("Reference structure/trajectory path, or auto for first selected frame of replica 1", "auto")
    config["reference"] = None if reference.lower() == "auto" else str(Path(reference).expanduser().resolve())
    default_image = "auto" if info["periodic"] is not None else "off"
    print("Imaging requires periodic box data. Multi-chain complexes may need a carefully chosen anchor.")
    config["imaging"]["mode"] = ask("Imaging: auto/on/off", default_image, str.lower, lambda s: s in {"auto", "on", "off"})
    if config["imaging"]["mode"] != "off":
        anchor = prompt_mask("Imaging anchor mask, or default for CPPTRAJ's first molecule", "default")
        config["imaging"]["anchor"] = None if anchor == "default" else anchor
    print("\nAnalysis levels:")
    for level, names in LEVELS.items():
        print(f"  {level}: {', '.join(sorted(names))}")
    print("Inapplicable protein/DNA-specific actions are omitted. Advanced PCA/clustering can require substantial RAM.")
    config["analysis"]["level"] = ask("Analysis level", "standard", str.lower, lambda x: x in LEVELS)
    additional_analyses(config, info)
    configure_convergence(config)
    configure_clustering(config, info)
    enabled = enabled_analyses(config)
    if "nastruct" in enabled:
        residues = ask("NAStruct residue range (e.g. 121-144), or all", "all")
        config["nucleic"]["resrange"] = None if residues.lower() == "all" else residues.lstrip(":")
        print("Modified bases can be mapped to canonical bases if the required atoms exist.")
        while yes("Add a modified nucleic residue mapping?", False):
            label = ask("Modified residue name")
            config["nucleic"]["resmap"][label] = ask("Canonical base A/C/G/T/U", None, str.upper, lambda s: s in {"A", "C", "G", "T", "U"})
    if "protein" in info["selections"] and "nucleic" in info["selections"]:
        print("Protein–nucleic interface distances, contacts and H-bonds will be included automatically when enabled.")
    while enabled & {"contacts", "hbond", "distance"} and yes("Add an explicit interface for enabled contact/H-bond/distance analyses?", False):
        config["interactions"].append({"name": ask("Interface name", None, safe_name),
                                       "mask1": prompt_mask("First named selection or mask"), "mask2": prompt_mask("Second named selection or mask")})
        if "contacts" in enabled:
            config["interactions"][-1]["contact_cutoff"] = ask("Contact cutoff in angstroms", 4.5, float, lambda x: x > 0)
    if enabled & {"dccm", "pca", "cluster", "pairwise_rmsd"}:
        config["advanced"]["matrix_mask"] = prompt_mask("Advanced coordinate mask (CA + nucleic C4' recommended)", config["advanced"]["matrix_mask"] if "cluster" in enabled else default_fit)
    if "pca" in enabled:
        config["advanced"]["pca_modes"] = ask("Number of principal components", 3, int, lambda x: x > 0)
    configure_resource_limits(config, length_replicas)
    config["stats"]["block_size"] = ask("Block size in analyzed frames for block-mean diagnostics", 50, int, lambda x: x >= 2)
    return finish_wizard(config, config_path, root)


def finish_wizard(config, config_path, root):
    from .console import panel
    panel("PUBLICATION FIGURES", [("Font", "Arial preferred; local fallback recorded if unavailable"),
          ("Layout", "180 mm canvas, 11 pt labels, editable SVG text"),
          ("Palette", "Lagoon: coordinated teal, copper, blue and violet"), ("Raster detail", "600 DPI for embedded raster layers")])
    if yes("Customize figure style?", False):
        config["plots"]["palette"] = ask("Figure palette: lagoon/mineral/colorblind", "lagoon", str.lower, lambda s: s in {"lagoon", "mineral", "colorblind"}, help_text="Choose a coordinated palette; use colorblind for the established categorical alternative.")
        config["plots"]["font"] = ask("Preferred font family", "Arial", help_text="The font must already exist on the rendering machine; any fallback is recorded.")
        config["plots"]["width_mm"] = ask("Figure canvas width in mm", 180, float, lambda n: 80 <= n <= 300, help_text="Use about 90 mm for a single-column canvas or 180 mm for a double-column canvas, then check your journal's dimensions.")
        config["plots"]["font_size"] = ask("Axis label font size in points", 11, float, lambda n: 7 <= n <= 24, help_text="Choose text for its intended final printed size rather than enlarging a small figure afterward.")
        config["plots"]["svg_text"] = ask("SVG text: editable/paths", "editable", str.lower, lambda s: s in {"editable", "paths"}, help_text="Editable text can be revised in a graphics editor; paths preserve rendered glyph shapes across machines.")
    config["cpptraj"] = ask("CPPTRAJ executable or path", "cpptraj")
    config["output"] = str(Path(ask("Results directory", str(root / "analysis_results"))).expanduser().resolve())
    target = Path(config_path).expanduser().resolve()
    config = normalize_config(config, Path.cwd())
    if target.exists():
        raise ValueError(f"Configuration already exists: {target}. Use --config with a new filename.")
    target.parent.mkdir(parents=True, exist_ok=True)
    save_config(config, target)
    print(f"\nConfiguration saved: {target}")
    return offer_run(target)


def offer_run(target):
    entry = "mdworkbench.pyz" if sys.argv[0].endswith(".pyz") else "mdworkbench.py"
    prefix = ["python", entry]
    if Path(sys.argv[0]).name in {"md-bud", "md-bud.exe", "md-bud.pyz"}:
        prefix = ["md-bud"] if sys.argv[0].endswith(".exe") else ["python3", sys.argv[0]]
    elif Path(sys.argv[0]).name == "__main__.py" and Path(sys.argv[0]).parent.name == "mdbud":
        prefix = ["python", "-m", "mdbud"]
    command = shlex.join([*prefix, "run", str(target)])
    print("Preview: " + shlex.join([*prefix, "plan", str(target), "--output", str(target.parent / "analysis_plan")]))
    print("Run saved configuration: " + command)
    if yes("Run CPPTRAJ now?", True):
        config = load_config(target)
        from .runner import run_workflow
        manifest = run_workflow(config)
        from .console import result_summary
        if "report" in manifest:
            result_summary(manifest["report"], config["output"])
        print(f"Run status: {manifest['status']}. Results: {config['output']}")
        return 0 if manifest["status"] == "complete" else 2
    print("Configuration kept. Run it later with:\n  " + command)
    return 0


def main(argv=None):
    parser = argparse.ArgumentParser(description="Replica-aware CPPTRAJ analysis with an interactive wizard and publication SVG reports.")
    parser.add_argument("--version", action="version", version=__version__)
    sub = parser.add_subparsers(dest="command")
    setup = sub.add_parser("wizard", help="Discover files and answer guided setup questions (default)")
    setup.add_argument("--root", default=".")
    setup.add_argument("--config", default="analysis_config.json")
    setup.add_argument("--quick", action="store_true", help="Use the short basic setup path with optional add-ons")
    scan = sub.add_parser("discover", help="List candidate topology and trajectory files")
    scan.add_argument("root", nargs="?", default=".")
    doctor = sub.add_parser("doctor", help="Check environment availability without installing or running calculations")
    doctor.add_argument("--cpptraj", default="cpptraj")
    doctor.add_argument("--no-plots", action="store_true", help="Do not require Matplotlib")
    validate = sub.add_parser("validate", help="Check files, selections and applicable batches without running CPPTRAJ")
    validate.add_argument("config")
    catalog = sub.add_parser("catalog", help="List analyses, presets and official command documentation")
    catalog.add_argument("analysis", nargs="?", choices=sorted(KNOWN_ANALYSES))
    bundle = sub.add_parser("bundle", help="Create a portable ZIP of generated results, excluding original simulation inputs")
    bundle.add_argument("output", help="Existing results directory")
    bundle.add_argument("--destination", required=True, help="New ZIP filename outside the results directory")
    for name in ("plan", "run"):
        action = sub.add_parser(name, help="Write an offline review plan" if name == "plan" else "Validate inputs, execute CPPTRAJ, and report")
        action.add_argument("config")
        action.add_argument("--output", help="New output directory (must be empty)")
        if name == "run":
            action.add_argument("--cpptraj", help="Override CPPTRAJ executable")
            action.add_argument("--no-plots", action="store_true", help="Write data/statistics without SVG rendering")
            action.add_argument("--no-progress", action="store_true", help="Disable terminal progress animation")
    report = sub.add_parser("report", help="Regenerate tables and plots from existing outputs, without rerunning CPPTRAJ")
    report.add_argument("output")
    report.add_argument("--no-plots", action="store_true")
    demo = sub.add_parser("demo", help="Generate clearly labeled synthetic data to exercise reporting without CPPTRAJ")
    demo.add_argument("--output", default="synthetic_demo")
    demo.add_argument("--no-plots", action="store_true")
    args = parser.parse_args(argv)
    try:
        if args.command is None:
            if not sys.stdin.isatty():
                parser.print_help()
                return 0
            return wizard(".", "analysis_config.json")
        if args.command == "wizard":
            return wizard(args.root, args.config, quick=True if args.quick else None)
        if args.command == "discover":
            print(json.dumps({k: [str(p) for p in v] for k, v in discover(args.root).items()}, indent=2))
            return 0
        if args.command == "doctor":
            from .community import environment_info
            info = environment_info(args.cpptraj)
            missing = []
            if not info["python_supported"]:
                missing.append("Python >=3.10")
            for name in (["numpy"] if args.no_plots else ["numpy", "matplotlib"]):
                if not info["meets_minimum_version"][name]:
                    missing.append(name + " (missing or below the documented minimum version)")
            if not info["cpptraj_executable"]:
                missing.append("CPPTRAJ executable")
            info["missing"] = missing
            info["next_step"] = "Use an existing environment with the missing prerequisites." if missing else "Validate your configuration; package importability and CPPTRAJ compatibility are checked during use."
            print(json.dumps(info, indent=2))
            return 2 if missing else 0
        if args.command == "validate":
            from .community import inspect_config
            print(json.dumps(inspect_config(load_config(args.config)), indent=2))
            return 0
        if args.command == "catalog":
            from .analysis import SOURCES
            descriptions = {
                "rdf": "Optional volume-normalized radial distribution for two disjoint atom selections",
                "watershell": "Optional solvent shell counts around a selected region",
                "rmsd": "Structural deviation from the common reference, with global and local fits",
                "rg": "Radius of gyration: overall compactness of each section",
                "rmsf": "Residue fluctuations after alignment",
                "hbond": "Hydrogen bonds in solute and configured interfaces",
                "dssp": "Protein secondary structure and state fractions",
                "sasa": "LCPO solvent-accessible surface area per section",
                "nastruct": "Nucleic base-pair, step, helix, groove and pucker geometry",
                "distance": "Distance between two ordered atom/group masks",
                "angle": "Angle between three masks, with mask 2 as the vertex",
                "dihedral": "Signed torsion of four ordered masks; circular statistics",
                "contacts": "Interface native/nonnative contacts and residue maps",
                "torsions": "Protein backbone and nucleic backbone torsion series",
                "dccm": "Aligned-coordinate dynamic cross-correlation matrices",
                "pca": "Shared pooled principal components and replica projections",
                "cluster": "Shared pooled structural clustering and replica populations",
                "stacking": "Confirmed bond-graph pi systems, per-pair and per-residue replica occupancies",
                "pairwise_rmsd": "Pairwise frame RMSD matrix; potentially expensive",
            }
            for name in ([args.analysis] if args.analysis else sorted(KNOWN_ANALYSES)):
                presets = ", ".join(level for level, members in LEVELS.items() if name in members) or "optional add-on"
                print(f"{name}: {descriptions[name]}\n  Presets: {presets}\n  {SOURCES[name]}")
            print("Use the wizard's numbered add-on menu for selection and parameter prompts. Standard excludes SASA and NAStruct.")
            return 0
        if args.command == "bundle":
            from .community import bundle_results
            print(json.dumps(bundle_results(args.output, args.destination), indent=2))
            return 0
        if args.command == "demo":
            from .demo import create_demo
            summary = create_demo(Path(args.output).resolve(), render=not args.no_plots)
            print(json.dumps(summary, indent=2))
            return 0
        if args.command == "report":
            from .runner import replot
            from .console import result_summary
            result_summary(replot(args.output, make_plots=not args.no_plots), args.output)
            return 0
        config = load_config(args.config)
        if args.output:
            config["output"] = str(Path(args.output).expanduser().resolve())
            config = normalize_config(config)
        if args.command == "plan":
            from .runner import make_plan
            print(f"Plan saved to {make_plan(config)}")
            return 0
        from .runner import run_workflow
        if args.cpptraj:
            config["cpptraj"] = args.cpptraj
        if args.no_progress:
            config["progress"]["enabled"] = False
        manifest = run_workflow(config, make_plots=not args.no_plots)
        from .console import result_summary
        if "report" in manifest:
            result_summary(manifest["report"], config["output"])
        print(f"Run status: {manifest['status']}. Results: {config['output']}")
        return 0 if manifest["status"] == "complete" else 2
    except (ValueError, OSError, RuntimeError, ImportError) as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 2
    except (KeyboardInterrupt, EOFError):
        print("\nStopped. Existing results and logs are retained.", file=sys.stderr)
        return 130


if __name__ == "__main__":
    raise SystemExit(main())
