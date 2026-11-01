"""Synthetic bond graphs and occupancies; these are not runnable MD topologies."""
import contextlib
import copy
import csv
import io
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import numpy as np
from mdwb.rings import chordless_cycles, discover_pi_systems, read_bond_graph
from mdwb.stacking import (DEFAULTS, stacking_pairs, pair_name, build_stacking_batches,
                           geometry_series, write_stacking_report, validate_stacking)
from mdwb.reporting import NumericTable, Series, report_results
from mdwb.config import normalize_config
from mdwb.cli import configure_stacking
from test_config_topology import amber_metadata


def ring(name, residue, first=1):
    return {'name':name, 'residue':residue, 'atoms':list(range(first,first+6))}


def options():
    return {**copy.deepcopy(DEFAULTS), 'rings':[ring('donor',1),ring('partner',2,7),ring('other',3,13)], 'targets':['donor']}


def amber_graph(path, residues, edges, types=None, numbers=None):
    text = amber_metadata(residues)
    n = sum(len(atoms) for _,atoms in residues)
    for flag, values, kind, width in (
        ('AMBER_ATOM_TYPE', types or ['ca']*n, 'a', 4),
        ('ATOMIC_NUMBER', numbers or [6]*n, 'I', 8),
        ('BONDS_WITHOUT_HYDROGEN',[x for a,b in edges for x in ((a-1)*3,(b-1)*3,1)],'I',8)):
        text += f'%FLAG {flag}\n%FORMAT(10{kind}{width})\n'
        for i in range(0,len(values),10):
            text += ''.join(f'{v:>{width}}' for v in values[i:i+10])+'\n'
    path.write_text(text)
    return path


def cycle(first=1, count=6):
    return [(first+i,first+(i+1)%count) for i in range(count)]


class RingDiscoveryTests(unittest.TestCase):
    def test_nonstandard_names_and_residue_numbers(self):
        with tempfile.TemporaryDirectory() as temp:
            path = amber_graph(Path(temp)/'lig.prmtop', [('X01',['Q7','Z1','K9','B8','L2','J0'])],cycle())
            found, _ = discover_pi_systems(path,{1})
        self.assertEqual(len(found),1)
        self.assertEqual(found[0]['atoms'],list(range(1,7)))
        self.assertIn('Q7',found[0]['atom_names'])
        self.assertEqual(found[0]['evidence'],'pi_atom_types')

    def test_saturated_sugar_ring_is_excluded(self):
        with tempfile.TemporaryDirectory() as temp:
            path = amber_graph(Path(temp)/'x.top',[('SUG',['X1','X2','X3','X4','X5'])],cycle(count=5),types=['CT']*5)
            found, notes = discover_pi_systems(path)
        self.assertFalse(found)
        self.assertTrue(notes)

    def test_lone_pair_oxygen_does_not_exclude_furan_candidate(self):
        with tempfile.TemporaryDirectory() as temp:
            path = amber_graph(Path(temp)/'furan.top',[('FUR',['X1','X2','X3','X4','X5'])],cycle(count=5),
                               types=['ca']*4+['os'],numbers=[6]*4+[8])
            found,_ = discover_pi_systems(path)
        self.assertEqual(len(found),1)

    def test_fused_rings_and_system_without_outer_perimeter_duplicates(self):
        # Two six-member rings share edge 5-6.
        edges = cycle()+[(5,7),(7,8),(8,9),(9,10),(10,6)]
        with tempfile.TemporaryDirectory() as temp:
            path = amber_graph(Path(temp)/'x.top',[('FUS',[f'X{i}' for i in range(10)])],edges)
            found,_ = discover_pi_systems(path)
        self.assertEqual([len(r['atoms']) for r in found],[6,6,10])
        self.assertEqual(found[-1]['kind'],'fused_system')

    def test_hydrogens_and_cross_residue_bonds_do_not_make_rings(self):
        with tempfile.TemporaryDirectory() as temp:
            path = amber_graph(Path(temp)/'x.top',[('A',['A1','A2','A3']),('B',['B1','B2','B3'])],cycle())
            found,_ = discover_pi_systems(path)
        self.assertFalse(found)

    def test_guanidinium_is_labeled_as_an_acyclic_candidate(self):
        with tempfile.TemporaryDirectory() as temp:
            path = amber_graph(Path(temp)/'x.top',[('UNK',['X','Y','Z','W'])],[(1,2),(1,3),(1,4)],types=['C','N2','N2','N2'],numbers=[6,7,7,7])
            found,_ = discover_pi_systems(path)
        self.assertEqual(found[0]['kind'],'guanidinium_candidate')

    def test_pdb_distance_inference_is_explicitly_unverified(self):
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp)/'ring.pdb'
            lines = []
            for i in range(6):
                angle = i*np.pi/3
                x,y = 1.4*np.cos(angle),1.4*np.sin(angle)
                lines.append(f'HETATM{i+1:5d}  X{i+1:<2} LIG A   1    {x:8.3f}{y:8.3f}{0.:8.3f}{1.:6.2f}{0.:6.2f}           C\n')
            path.write_text(''.join(lines))
            found,_ = discover_pi_systems(path)
        self.assertEqual(len(found),1)
        self.assertIn('inferred',found[0]['evidence'])

    def test_mol2_noncontiguous_serials_map_to_topology_atom_order(self):
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp)/'ring.mol2'
            text = '@<TRIPOS>MOLECULE\nLIG\n6 6 1\nSMALL\nNO_CHARGES\n@<TRIPOS>ATOM\n'
            text += ''.join(f'{10*(i+1)} Q{i} 0 0 0 C.ar 42 LIG\n' for i in range(6))
            text += '@<TRIPOS>BOND\n'+''.join(f'{i+1} {a*10} {b*10} ar\n' for i,(a,b) in enumerate(cycle()))
            path.write_text(text)
            found,_ = discover_pi_systems(path)
        self.assertEqual(found[0]['atoms'],list(range(1,7)))
        self.assertEqual(found[0]['residue'],1)

    def test_bad_bond_index_is_rejected(self):
        with tempfile.TemporaryDirectory() as temp:
            path = amber_graph(Path(temp)/'x.top',[('A',['X']*6)],[(1,50)])
            with self.assertRaisesRegex(ValueError,'bond indices'):
                read_bond_graph(path)


class StackingTests(unittest.TestCase):
    def config(self):
        return {'stacking':options(), 'analysis':{'enabled':['stacking']},
                'plots':{'enabled':False},'frames':{'dt_ps':1000,'start':1,'stride':1}}

    def read(self,root,name):
        with (Path(root)/'reports/stacking'/name).open() as handle:
            return list(csv.DictReader(handle,delimiter='\t'))

    def occupancy(self,replica,pair,values):
        return Series(replica,'stacking',pair,'occupancy','timeseries',np.arange(len(values),dtype=float),np.array(values,dtype=float),'fraction','time (ns)','synthetic')

    def test_all_pairs_unique_no_neighbor_prefilter_or_overlapping_self_pair(self):
        opt = options()
        opt['targets'] = ['donor','partner']
        self.assertEqual(len(stacking_pairs(opt)),3)
        opt['rings'].append({'name':'subring','residue':1,'atoms':[1,2,3,4,5]})
        self.assertFalse(any(set(a['atoms'])&set(b['atoms']) for a,b in stacking_pairs(opt)))
        opt['max_pairs'] = 1
        with self.assertRaises(ValueError):
            stacking_pairs(opt)

    def test_geometry_pipeline_preserves_box_and_defers_vector_math(self):
        batch = build_stacking_batches(self.config())[0]
        self.assertTrue(batch['preserve_box'])
        self.assertTrue(all('corrplane' in cmd for cmd in batch['commands'][:3]))
        self.assertEqual(sum(cmd.startswith('distance ') for cmd in batch['commands']),2)
        self.assertTrue(all('geom' in cmd for cmd in batch['commands'] if cmd.startswith('distance ')))
        self.assertTrue(batch['post_commands'][0].startswith('runanalysis vectormath'))
        self.assertTrue(batch['post_commands'][1].startswith('writedata'))
        self.assertNotIn('maskout','\n'.join(batch['commands']))

    def test_parallel_antiparallel_cutoff_and_invalid_geometry(self):
        artifact = build_stacking_batches(self.config())[0]['artifacts'][0]
        table = NumericTable(np.array([[1,5,30],[2,5,150],[3,5.01,0],[4,4,90],[5,4,np.nan],[6,-1,0]]),['Frame','D','A'],[])
        result = geometry_series('A',artifact,table,self.config(),'synthetic')
        np.testing.assert_allclose(result[0].values,[1,1,0,0,np.nan,np.nan],equal_nan=True)
        self.assertEqual(result[0].unit,'fraction')
        self.assertEqual(result[2].values[1],30)

    def test_geometry_cutoffs_can_be_changed_when_replotting(self):
        cfg = self.config()
        artifact = build_stacking_batches(cfg)[0]['artifacts'][0]
        cfg['stacking']['distance_cutoff'] = 4
        table = NumericTable(np.array([[1,4.5,10]]),['Frame','D','A'],[])
        self.assertEqual(geometry_series('A',artifact,table,cfg,'synthetic')[0].values[0],0)

    def test_replot_rejects_changed_centers_or_atom_membership(self):
        cfg = self.config()
        artifact = build_stacking_batches(cfg)[0]['artifacts'][0]
        table = NumericTable(np.array([[1,4,10]]),['Frame','D','A'],[])
        cfg['stacking']['center'] = 'mass'
        with self.assertRaisesRegex(ValueError,'rerunning CPPTRAJ'):
            geometry_series('A',artifact,table,cfg,'synthetic')
        cfg['stacking']['center'] = 'geometry'
        cfg['stacking']['rings'][0]['atoms'] = list(range(20,26))
        with self.assertRaisesRegex(ValueError,'rerunning CPPTRAJ'):
            geometry_series('A',artifact,table,cfg,'synthetic')

    def test_same_residue_disjoint_rings_can_be_enabled_explicitly(self):
        opt = options()
        opt['rings'] = [ring('donor',1),ring('partner',1,7)]
        self.assertFalse(stacking_pairs(opt))
        opt['exclude_same_residue'] = False
        self.assertEqual(len(stacking_pairs(opt)),1)

    def test_job_pair_batching_reuses_vectors_and_bounds_job_size(self):
        cfg = self.config()
        cfg['stacking']['rings'] = [ring('donor',1)] + [ring(f'partner{i}',i+2,7+i*6) for i in range(40)]
        batches = build_stacking_batches(cfg)
        self.assertEqual([len(b['artifacts']) for b in batches],[32,8])
        self.assertEqual([b['stacking_ring_count'] for b in batches],[33,9])
        self.assertEqual(len({a['section'] for b in batches for a in b['artifacts']}),40)

    def test_equal_replica_mean_sd_and_any_partner_union(self):
        cfg = self.config(); pairs = stacking_pairs(cfg['stacking'])
        keys = [pair_name(*p) for p in pairs]
        data = [self.occupancy('A',keys[0],[1,0,1,0]), self.occupancy('A',keys[1],[0,1,1,0]),
                self.occupancy('B',keys[0],[0]*8), self.occupancy('B',keys[1],[0]*8)]
        with tempfile.TemporaryDirectory() as temp:
            result, unions = write_stacking_report(cfg,data,{'replicas':[{'name':'A'},{'name':'B'}]},temp)
            rows = self.read(temp,'replica_summary.dat')
        pair = next(r for r in rows if r['system']==keys[0])
        self.assertAlmostEqual(float(pair['mean_occupancy']),.25)
        self.assertAlmostEqual(float(pair['sd_replica_occupancies']),np.std([.5,0],ddof=1))
        target = next(r for r in rows if r['scope']=='target')
        self.assertAlmostEqual(float(target['mean_occupancy']),.375)
        self.assertEqual(result['omissions'],0)
        a = next(s for s in unions if s.replica=='A' and s.section=='target_donor')
        np.testing.assert_array_equal(a.values,[1,1,1,0])

    def test_missing_replica_is_not_zero_and_one_replica_sd_is_blank(self):
        cfg = self.config(); keys = [pair_name(*p) for p in stacking_pairs(cfg['stacking'])]
        data = [self.occupancy('A',key,[1,1,0,0]) for key in keys]
        with tempfile.TemporaryDirectory() as temp:
            result,_ = write_stacking_report(cfg,data,{'replicas':[{'name':'A'},{'name':'B'}]},temp)
            rows = self.read(temp,'replica_summary.dat')
        self.assertTrue(result['omissions'])
        self.assertTrue(all(r['sd_replica_occupancies']=='' for r in rows))
        self.assertTrue(all(float(r['mean_occupancy'])==.5 for r in rows))
        self.assertTrue(all(r['status']=='incomplete_replicas' for r in rows))

    def test_nonfinite_geometry_excluded_from_replica_aggregate(self):
        cfg = self.config(); keys = [pair_name(*p) for p in stacking_pairs(cfg['stacking'])]
        data = [self.occupancy('A',key,[1,np.nan,0,0]) for key in keys]
        with tempfile.TemporaryDirectory() as temp:
            _,_ = write_stacking_report(cfg,data,{'replicas':[{'name':'A'}]},temp)
            rows = self.read(temp,'replica_summary.dat')
        self.assertTrue(all(r['mean_occupancy']=='' for r in rows))

    def test_configuration_validation_and_atom_residue_membership(self):
        with tempfile.TemporaryDirectory() as temp:
            path = amber_graph(Path(temp)/'x.top',[('A',[f'X{i}' for i in range(6)]),('B',[f'Y{i}' for i in range(6)])],cycle()+cycle(7))
            traj = Path(temp)/'a.nc'; traj.write_text('fixture')
            cfg = {'topology':str(path),'replicas':[{'name':'A','trajectories':[str(traj)]}],
                   'analysis':{'enabled':['stacking']}, 'stacking':{**options(),'rings':options()['rings'][:2]}}
            normalize_config(cfg,base=temp)
            cfg['stacking']['rings'][0]['residue']=3
            with self.assertRaisesRegex(ValueError,'do not belong'):
                normalize_config(cfg,base=temp)
        opt = options(); opt['angle_cutoff']=float('nan')
        with self.assertRaises(ValueError):
            validate_stacking(opt,True)

    def test_wizard_confirms_and_names_without_asking_for_atom_names(self):
        with tempfile.TemporaryDirectory() as temp:
            path = amber_graph(Path(temp)/'x.top',[('LIG',[f'Q{i}' for i in range(6)]),('MOD',[f'Z{i}' for i in range(6)])],cycle()+cycle(7))
            cfg = self.config(); cfg['topology']=str(path)
            answers = ['1','','donor_ring','n','all','','n','5','30','geometry']
            with patch('builtins.input',side_effect=answers) as prompt, contextlib.redirect_stdout(io.StringIO()):
                configure_stacking(cfg)
        self.assertEqual(cfg['stacking']['targets'],['donor_ring'])
        self.assertEqual(len(cfg['stacking']['rings']),2)
        self.assertFalse(any('atom name' in c.args[0].lower() for c in prompt.call_args_list))

    def test_report_end_to_end_geometry_occupancy_html_and_convergence(self):
        cfg = self.config()
        cfg['stacking']['rings']=cfg['stacking']['rings'][:2]
        artifacts = build_stacking_batches(cfg)[0]['artifacts']
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp); replicas=[]
            for rep,dist in [('A',4),('B',8)]:
                artifact=copy.deepcopy(artifacts[0]); path=root/f'{rep}.dat'
                path.write_text('#Frame D A\n'+''.join(f'{i+1} {dist} 170\n' for i in range(50)))
                artifact['path']=str(path)
                replicas.append({'name':rep,'status':'complete','frame_count':50,'artifacts':[artifact]})
            result=report_results(cfg,{'replicas':replicas},root)
            self.assertTrue(result['stacking']['enabled'])
            self.assertEqual(result['stacking']['omissions'],0)
            rows=self.read(root,'replica_summary.dat')
            self.assertAlmostEqual(float(rows[0]['mean_percent']),50)
            self.assertIn('Pi-stacking occupancies',(root/'report_pages/replicas.html').read_text())
            self.assertIn('stacking_any',(root/'reports/convergence/within_replica.dat').read_text())


if __name__ == '__main__':
    unittest.main()
