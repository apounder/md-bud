"""Bond-graph pi-system candidates without residue/atom-name templates.

Connectivity identifies cycles, not aromaticity. Atom types provide conservative
screening; users confirm candidate chemistry before enabling stacking analyses.
"""
from __future__ import annotations

import gzip
import math
from pathlib import Path
from .topology import _amber_flags, inspect_topology

# Force-field hybridization/type hints, never atom names or residue templates.
SATURATED = {'CT', 'CX', 'CI', 'c3', 'c4', 'c5', 'c6', 'n3', 'n4',
             'C.3', 'N.3', 'N.4'}
PI_TYPES = {'CA', 'CB', 'CC', 'CD', 'CK', 'CM', 'CN', 'CQ', 'CR', 'CV', 'CW',
            'C', 'N', 'NA', 'NB', 'NC', 'N*', 'ca', 'cp', 'cq', 'cc', 'cd', 'ce',
            'cf', 'c', 'c2', 'na', 'nb', 'nc', 'nd', 'ne', 'nf', 'n', 'n2', 'nh',
            'C.ar', 'N.ar', 'C.2', 'N.2', 'N.pl3', 'S.2'}
ELEMENTS = {1:'H', 6:'C', 7:'N', 8:'O', 9:'F', 15:'P', 16:'S', 17:'Cl', 35:'Br', 53:'I'}


def read_bond_graph(path):
    """Return ordered topology atoms and undirected bonds; no new dependencies."""
    path = Path(path)
    opener = gzip.open if path.suffix.lower() == '.gz' else open
    with opener(path, 'rt', encoding='utf-8', errors='replace') as handle:
        head = handle.read(256)
        handle.seek(0)
        lines = [] if '%FLAG' in head or '%VERSION' in head else handle.readlines()
    atoms, bonds = [], set()
    source = 'topology_bonds'
    if '%FLAG' in head or '%VERSION' in head:
        flags = _amber_flags(path, {'BONDS_INC_HYDROGEN', 'BONDS_WITHOUT_HYDROGEN',
                                     'AMBER_ATOM_TYPE', 'ATOMIC_NUMBER', 'MASS'})
        info = inspect_topology(path)
        types = flags.get('AMBER_ATOM_TYPE', [])
        numbers = flags.get('ATOMIC_NUMBER', [])
        masses = flags.get('MASS', [])
        for residue in info['residues']:
            for name in residue['atoms']:
                index = len(atoms)
                element = ELEMENTS.get(int(numbers[index]), '') if len(numbers) == info['natom'] else ''
                if not element and len(masses) == info['natom']:
                    mass = float(masses[index].replace('D','E'))
                    element = next((e for e, m in [('H',1.008), ('C',12.01), ('N',14.01), ('O',16), ('S',32.06)] if abs(mass-m) < .7), '')
                atoms.append(dict(index=index+1, name=name, residue=residue['index'],
                                  residue_name=residue['name'], element=element,
                                  atom_type=types[index] if len(types) == info['natom'] else ''))
        for key in ('BONDS_INC_HYDROGEN', 'BONDS_WITHOUT_HYDROGEN'):
            values = [int(v) for v in flags.get(key, [])]
            if len(values) % 3:
                raise ValueError(f'Malformed Amber {key} bond triples')
            for i in range(0, len(values), 3):
                a, b = values[i:i+2]
                if a < 0 or b < 0 or a % 3 or b % 3:
                    raise ValueError('Malformed Amber bond atom offset')
                bonds.add(tuple(sorted((a//3+1, b//3+1))))
    elif any(line.startswith('@<TRIPOS>ATOM') for line in lines):
        section, ids, residue_ids = '', {}, {}
        pending = []
        for line in lines:
            if line.startswith('@<TRIPOS>'):
                section = line.strip(); continue
            parts = line.split()
            if not parts or line.lstrip().startswith('#'):
                continue
            if section == '@<TRIPOS>ATOM':
                if len(parts) < 8:
                    raise ValueError('MOL2 ring detection requires atom residue identifiers')
                key = parts[6]
                residue_ids.setdefault(key, len(residue_ids)+1)
                if parts[0] in ids:
                    raise ValueError('MOL2 must contain one topology molecule with unique atom IDs')
                ids[parts[0]] = len(atoms)+1
                atoms.append(dict(index=len(atoms)+1, name=parts[1], residue=residue_ids[key],
                                  residue_name=parts[7], element=parts[5].split('.')[0], atom_type=parts[5]))
            elif section == '@<TRIPOS>BOND':
                pending.append((parts[1], parts[2]))
        for a,b in pending:
            if a not in ids or b not in ids:
                raise ValueError('MOL2 bond references an unknown atom')
            bonds.add(tuple(sorted((ids[a],ids[b]))))
    elif any(line.startswith(('ATOM  ', 'HETATM')) for line in lines):
        ids, last_key, residue = {}, None, 0
        conect, finished_model = [], False
        for line in lines:
            if line.startswith('ENDMDL'):
                finished_model = True
            if line.startswith('TER'):
                last_key = None
            if line.startswith(('ATOM  ', 'HETATM')) and not finished_model:
                key = line[21:27], line[17:20]
                if key != last_key:
                    residue += 1; last_key = key
                serial = int(line[6:11])
                if serial in ids or line[16:17] not in (' ', ''):
                    raise ValueError('Ring detection requires a PDB with unique atoms and resolved alternate locations')
                element = line[76:78].strip().title()
                if not element:
                    element = line[12:16].strip().lstrip('0123456789')[:1].upper()
                ids[serial] = len(atoms)+1
                atoms.append(dict(index=len(atoms)+1, name=line[12:16].strip(), residue=residue,
                                  residue_name=line[17:20].strip(), element=element, atom_type='',
                                  xyz=tuple(float(line[i:i+8]) for i in (30,38,46))))
            if line.startswith('CONECT'):
                values = [int(line[i:i+5]) for i in range(6,len(line.rstrip()),5) if line[i:i+5].strip()]
                conect.extend((values[0], value) for value in values[1:])
        for a,b in conect:
            if a in ids and b in ids:
                bonds.add(tuple(sorted((ids[a],ids[b]))))
        # Many biomolecular PDBs omit CONECT. Infer within-residue heavy bonds,
        # mark candidates as unverified, and require confirmation in the wizard.
        explicit_atoms = {a for pair in bonds for a in pair}
        radii = {'C':.76, 'N':.71, 'O':.66, 'S':1.05}
        by_residue = {}
        for atom in atoms:
            by_residue.setdefault(atom['residue'], []).append(atom)
        inferred = False
        for members in by_residue.values():
            heavy = [a for a in members if a['element'] in radii]
            if any(a['index'] in explicit_atoms for a in heavy):
                continue  # Do not silently repair a partially specified bond graph.
            if len(heavy) > 256:
                raise ValueError('PDB residue exceeds 256 heavy atoms; use a bonded Amber/MOL2 topology')
            for i,a in enumerate(heavy):
                for b in heavy[i+1:]:
                    distance = math.dist(a['xyz'], b['xyz'])
                    if .7 < distance <= radii[a['element']]+radii[b['element']]+.25:
                        bonds.add((a['index'],b['index'])); inferred = True
        if inferred:
            source = 'pdb_distance_inferred_bonds_unverified'
    else:
        raise ValueError('Automatic ring discovery supports bonded Amber, MOL2, or PDB topologies')
    if not atoms:
        raise ValueError('No topology atoms available for ring discovery')
    if any(a == b or a < 1 or b > len(atoms) for a,b in bonds):
        raise ValueError('Invalid topology bond indices')
    return atoms, bonds, source


def chordless_cycles(graph, minimum=5, maximum=7):
    """Small simple rings only; fused perimeter cycles with chords are excluded."""
    cycles = set()
    visits = 0
    def visit(path, start):
        nonlocal visits
        visits += 1
        if visits > 100000:
            raise ValueError('Ring search exceeds graph complexity limit')
        for atom in sorted(graph[path[-1]]):
            if atom == start:
                if len(path) >= minimum:
                    ring = frozenset(path)
                    if all(len(graph[n] & ring) == 2 for n in ring):
                        cycles.add(tuple(sorted(ring)))
                continue
            if len(path) < maximum and atom > start and atom not in path:
                visit(path+[atom], start)
    for atom in sorted(graph):
        visit([atom], atom)
    return sorted(cycles)


def discover_pi_systems(path, residues=None):
    atoms, bonds, source = read_bond_graph(path)
    by_id = {a['index']:a for a in atoms}
    groups = {}
    for atom in atoms:
        if residues is None or atom['residue'] in residues:
            groups.setdefault(atom['residue'], []).append(atom)
    adjacency = {a['index']:set() for a in atoms}
    for a,b in bonds:
        adjacency[a].add(b); adjacency[b].add(a)
    candidates, notes = [], []
    for residue, members in groups.items():
        eligible = {a['index'] for a in members if a['element'] in {'C','N','O','S'} or (not a['element'] and a['atom_type'] in PI_TYPES)}
        if len(eligible) > 256:
            raise ValueError(f'Residue {residue} exceeds the 256-heavy-atom ring search limit')
        graph = {a:adjacency[a] & eligible for a in eligible}
        rings = [ring for ring in chordless_cycles(graph) if not any(by_id[a]['atom_type'] in SATURATED for a in ring)]
        # Offer both constituent rings and each connected fused system.
        systems = [(ring, 'ring') for ring in rings]
        remaining = list(map(set, rings))
        while remaining:
            merged, count = remaining.pop(0), 1
            changed = True
            while changed:
                changed = False
                for ring in remaining[:]:
                    if len(merged & ring) >= 2:
                        merged |= ring; remaining.remove(ring); count += 1; changed = True
            if count > 1:
                systems.append((tuple(sorted(merged)), 'fused_system'))
        # Guanidinium is an acyclic conjugated pi group present in the supplied
        # donor script, not an aromatic ring. Connectivity identifies its candidate.
        for atom in sorted(graph):
            ns = [n for n in graph[atom] if by_id[n]['element'] == 'N']
            if by_id[atom]['element'] == 'C' and len(ns) == 3 and by_id[atom]['atom_type'] not in SATURATED:
                systems.append((tuple(sorted([atom,*ns])), 'guanidinium_candidate'))
        for index, (system, kind) in enumerate(systems, 1):
            evidence = 'pi_atom_types' if all(by_id[a]['atom_type'] in PI_TYPES for a in system) else 'connectivity_only_confirm_chemistry'
            if source.startswith('pdb_distance'):
                evidence = source
            candidates.append(dict(name=f'r{residue}_{kind}{index}', residue=residue,
                                   residue_name=members[0]['residue_name'], atoms=list(system),
                                   atom_names=[by_id[a]['name'] for a in system], kind=kind,
                                   evidence=evidence))
        if not systems:
            notes.append(f'Residue {residue} ({members[0]["residue_name"]}): no candidate 5–7 member pi rings or guanidinium groups')
    return candidates, notes
