"""Confirmed pi-system pairs, CPPTRAJ geometry, and replica occupancy summaries."""
from __future__ import annotations

import itertools
import math

DEFAULTS = {'rings': [], 'targets': [], 'distance_cutoff': 5.0, 'angle_cutoff': 30.0,
            'center': 'geometry', 'exclude_same_residue': True, 'max_pairs': 2000}


def stacking_pairs(options):
    targets = set(options['targets'])
    rings = options['rings']
    pairs = []
    for a,b in itertools.combinations(rings, 2):
        if a['name'] not in targets and b['name'] not in targets:
            continue
        if set(a['atoms']) & set(b['atoms']):
            continue
        if options.get('exclude_same_residue', True) and a['residue'] == b['residue']:
            continue
        pairs.append((a,b))
    if len(pairs) > options.get('max_pairs', 2000):
        raise ValueError(f'Stacking requests {len(pairs)} pairs; reduce selected rings or raise stacking.max_pairs explicitly')
    return pairs


def validate_stacking(options, enabled):
    from .config import safe_name
    names, atom_sets = set(), set()
    if not isinstance(options['rings'], list) or not isinstance(options['targets'], list):
        raise ValueError('stacking.rings and stacking.targets must be lists')
    for ring in options['rings']:
        if not isinstance(ring, dict) or set(ring) - {'name','residue','residue_name','atoms','atom_names','kind','evidence'}:
            raise ValueError('Stacking ring fields: name, residue, atoms, residue_name, atom_names, kind, evidence')
        name = safe_name(ring.get('name'))
        if name.casefold() in names:
            raise ValueError('Stacking ring names must be unique')
        names.add(name.casefold())
        if type(ring.get('residue')) is not int or ring['residue'] < 1:
            raise ValueError('Stacking residues are positive topology residue indices')
        atoms = ring.get('atoms')
        if not isinstance(atoms, list) or not 3 <= len(atoms) <= 128 or any(type(a) is not int or a < 1 for a in atoms) or len(set(atoms)) != len(atoms):
            raise ValueError('Each pi system requires 3–128 distinct positive topology atom indices')
        token = tuple(sorted(atoms))
        if token in atom_sets:
            raise ValueError('The same pi-system atom set cannot be configured twice')
        atom_sets.add(token)
    if any(not isinstance(n,str) or n not in {r['name'] for r in options['rings']} for n in options['targets']) or len(set(options['targets'])) != len(options['targets']):
        raise ValueError('stacking.targets must name distinct confirmed pi systems')
    for key, maximum in (('distance_cutoff',100), ('angle_cutoff',90)):
        value = options[key]
        if type(value) not in (float,int) or not math.isfinite(value) or not 0 < value <= maximum:
            raise ValueError(f'stacking.{key} must be finite, positive and <= {maximum}')
    if options['center'] not in {'geometry','mass'} or type(options['exclude_same_residue']) is not bool:
        raise ValueError('Stacking center must be geometry/mass; exclude_same_residue must be boolean')
    if type(options['max_pairs']) is not int or not 1 <= options['max_pairs'] <= 10000:
        raise ValueError('stacking.max_pairs must be an integer from 1 to 10000')
    if enabled and (not options['targets'] or not stacking_pairs(options)):
        raise ValueError('Stacking requires a target pi system and at least one nonoverlapping partner')


def pair_name(a,b):
    from .analysis import safe_name
    # Length-prefix encoding avoids collisions between names containing underscores.
    return safe_name(f"p{len(a['name'])}_{a['name']}_{len(b['name'])}_{b['name']}")


def build_stacking_batches(config):
    from .analysis import _artifact, _batch, quote_mask
    options = {**DEFAULTS, **config.get('stacking', {})}
    pairs = stacking_pairs(options)
    batches = []
    # Bound each job's retained vector/scalar datasets while amortizing trajectory I/O.
    for offset in range(0, len(pairs), 32):
        chunk = pairs[offset:offset+32]
        rings = {r['name']:r for pair in chunk for r in pair}
        masks = {name:'@'+','.join(map(str, r['atoms'])) for name,r in rings.items()}
        vectors = {name:f'A_PI_V{i}' for i,name in enumerate(rings)}
        commands = [f'vector {vectors[name]} corrplane {quote_mask(masks[name])}' for name in rings]
        post, artifacts = [], []
        for i,(a,b) in enumerate(chunk):
            token = pair_name(a,b)
            distance, angle = f'A_PI_D{i}', f'A_PI_A{i}'
            geom = 'geom' if options['center'] == 'geometry' else ''
            commands.append(f'distance {distance} {quote_mask(masks[a["name"]])} {quote_mask(masks[b["name"]])} {geom}')
            post.append(f'runanalysis vectormath vec1 {vectors[a["name"]]} vec2 {vectors[b["name"]]} dotangle name {angle}')
            artifact = _artifact('stacking', token, 'timeseries', 'fraction',
                                 format='stacking_geometry', ring_a=a['name'], ring_b=b['name'],
                                 ring_a_atoms=a['atoms'], ring_b_atoms=b['atoms'],
                                 distance_cutoff=options['distance_cutoff'], angle_cutoff=options['angle_cutoff'],
                                 center=options['center'], units=['angstrom','degree'])
            post.append(f'writedata {artifact["path"]} {distance} {angle}')
            artifacts.append(artifact)
        atoms = sorted({a for r in rings.values() for a in r['atoms']})
        batches.append(_batch(f'pi_stacking_{offset//32+1}',
                              'Confirmed pi systems: centroid distance and sign-independent plane angle, across all selected frames.',
                              commands, artifacts, post, preserve_box=True,
                              coordinate_mask='@'+','.join(map(str, atoms)), resource_class='stacking',
                              stacking_ring_count=len(rings), stacking_pair_count=len(chunk)))
    return batches


def geometry_series(replica, artifact, table, config, source):
    """Derive occupancy without dependence on CPPTRAJ calcstate state numbering."""
    import numpy as np
    from .reporting import Series, TableFormatError, _coordinate
    x, xunit, column = _coordinate(table, artifact, config, 'timeseries')
    if column != 0 or table.data.shape[1] != 3:
        raise TableFormatError('Stacking geometry requires Frame, distance, normal angle columns')
    distance, angle = table.data[:,1], table.data[:,2]
    valid = np.isfinite(distance) & np.isfinite(angle) & (distance >= 0) & (angle >= 0) & (angle <= 180)
    folded = np.minimum(angle, 180-angle)
    recorded = artifact['metadata']
    requested = config.get('stacking', {})
    if requested.get('center', recorded.get('center')) != recorded.get('center'):
        raise TableFormatError('Changing stacking center convention requires rerunning CPPTRAJ')
    if 'rings' in requested:
        rings = {r['name']:r for r in requested['rings']}
        for key in ('ring_a','ring_b'):
            ring = rings.get(recorded[key])
            if ring is None or (key+'_atoms' in recorded and set(ring['atoms']) != set(recorded[key+'_atoms'])):
                raise TableFormatError('Changing pi-system membership requires rerunning CPPTRAJ')
    metadata = {**recorded, **requested}
    values = np.full(len(distance), np.nan)
    values[valid] = ((distance[valid] <= metadata['distance_cutoff']) & (folded[valid] <= metadata['angle_cutoff'])).astype(float)
    kwargs = dict(replica=replica, section=artifact['section'], kind='timeseries', x=x, xunit=xunit, source=source)
    return [Series(analysis='stacking', name='occupancy', values=values, unit='fraction', **kwargs),
            Series(analysis='stacking_distance', name='centroid_distance', values=np.where(valid,distance,np.nan), unit='angstrom', **kwargs),
            Series(analysis='stacking_angle', name='plane_angle', values=np.where(valid,folded,np.nan), unit='degree', **kwargs)]


def write_stacking_report(config, series_list, manifest, output):
    """Equal-weight replica fractions; union occupancies never sum pair contacts."""
    from pathlib import Path
    import numpy as np
    from .diagnostics import _write
    from .reporting import Series
    from .analysis import enabled_analyses
    options = {**DEFAULTS, **config.get('stacking', {})}
    if 'stacking' not in enabled_analyses(config):
        return {'enabled':False}, []
    pairs = stacking_pairs(options)
    if not pairs:
        return {'enabled':False}, []
    folder = Path(output)/'reports'/'stacking'
    folder.mkdir(parents=True, exist_ok=True)
    replicas = [r['name'] for r in manifest.get('replicas', [])]
    lookup = {(s.replica,s.section):s for s in series_list if s.analysis == 'stacking' and s.name == 'occupancy'}
    # Include any input replicas omitted from the manifest as unavailable, not zero.
    replicas = list(dict.fromkeys([*replicas, *[r['name'] for r in config.get('replicas',[])]]))
    within, aggregate, unions, missing = [], [], [], []
    labels = {pair_name(a,b): f"{a['name']} / {b['name']}" for a,b in pairs}
    def summarize(scope, label, replica, s):
        if s is None:
            within.append([scope,label,replica,0,0,None,None,'missing'])
            missing.append([scope,label,replica,'missing_geometry'])
            return None
        finite = s.values[np.isfinite(s.values)]
        occupancy = float(finite.mean()) if len(finite) else None
        status = 'complete' if len(finite) == len(s.values) and len(finite) else 'incomplete'
        # Do not aggregate partial sampling as an unbiased whole-replica occupancy.
        within.append([scope,label,replica,len(s.values),len(finite),occupancy,
                       occupancy*100 if occupancy is not None else None,status])
        if status != 'complete':
            missing.append([scope,label,replica,'nonfinite_or_empty_geometry'])
            return None
        return occupancy
    def summary(scope,label,values):
        available = [v for v in values if v is not None]
        mean = float(np.mean(available)) if available else None
        sd = float(np.std(available,ddof=1)) if len(available)>1 else None
        aggregate.append([scope,label,len(replicas),len(available),mean,sd,
                          mean*100 if mean is not None else None, sd*100 if sd is not None else None,
                          'complete' if len(available)==len(replicas) else 'incomplete_replicas'])
    for a,b in pairs:
        token = pair_name(a,b)
        summary('pair',token,[summarize('pair',token,rep,lookup.get((rep,token))) for rep in replicas])
    # Per-target and per-target-residue occupancy: at least one selected contact.
    scopes = [('target',target,[pair_name(a,b) for a,b in pairs if target in (a['name'],b['name'])]) for target in options['targets']]
    residues = sorted({r['residue'] for r in options['rings'] if r['name'] in options['targets']})
    for residue in residues:
        names = {r['name'] for r in options['rings'] if r['residue']==residue and r['name'] in options['targets']}
        scopes.append(('residue',f'residue_{residue}',[pair_name(a,b) for a,b in pairs if a['name'] in names or b['name'] in names]))
    for scope,label,tokens in scopes:
        values = []
        for replica in replicas:
            members = [lookup.get((replica,token)) for token in tokens]
            union = None
            if members and all(s is not None for s in members) and all(np.array_equal(s.x,members[0].x) for s in members):
                valid = np.ones(len(members[0].values), dtype=bool)
                occupied = np.zeros(len(members[0].values), dtype=bool)
                for member in members:
                    valid &= np.isfinite(member.values)
                    occupied |= member.values == 1
                occupancy = np.where(valid,occupied.astype(float),np.nan)
                first = members[0]
                union = Series(replica,'stacking_any',f'{scope}_{label}','occupancy','timeseries',first.x,occupancy,'fraction',first.xunit,'derived: stacking geometry union')
                unions.append(union)
            values.append(summarize(scope,label,replica,union))
        summary(scope,label,values)
    _write(folder/'per_replica.dat',['scope','system','replica','n_frames','n_valid','occupancy','occupancy_percent','status','display_name'],
           [[*row,labels.get(row[1],row[1])] for row in within])
    _write(folder/'replica_summary.dat',['scope','system','n_replicas_requested','n_replicas_complete','mean_occupancy','sd_replica_occupancies','mean_percent','sd_percentage_points','status','display_name'],
           [[*row,labels.get(row[1],row[1])] for row in aggregate])
    _write(folder/'pairs.dat',['pair','ring_a','ring_b','residue_a','residue_b'],[[pair_name(a,b),a['name'],b['name'],a['residue'],b['residue']] for a,b in pairs])
    _write(folder/'rings.dat',['ring','residue','residue_name','kind','evidence','topology_atom_indices','atom_names'],
           [[r['name'],r['residue'],r.get('residue_name',''),r.get('kind','confirmed'),r.get('evidence','user_configured'),','.join(map(str,r['atoms'])),','.join(r.get('atom_names',[]))] for r in options['rings']])
    _write(folder/'omissions.dat',['scope','system','replica','reason'],missing)
    _write(folder/'criteria.dat',['distance_cutoff_angstrom','folded_plane_angle_cutoff_degree','center','distance_imaging','definition'],
           [[options['distance_cutoff'],options['angle_cutoff'],options['center'],'CPPTRAJ default minimum image when box exists','distance <= cutoff AND min(angle,180-angle) <= cutoff']])
    _write(folder/'any_contact_timeseries.dat',['replica','system','coordinate','coordinate_unit','occupancy'],
           ([s.replica,s.section,x,s.xunit,v] for s in unions for x,v in zip(s.x,s.values)))
    return {'enabled':True,'pairs':len(pairs),'replicas':len(replicas),'omissions':len(missing)}, unions
