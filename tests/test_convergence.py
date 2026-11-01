"""Synthetic scientific counterexamples and reference/time-axis invariants."""
import csv
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch
import numpy as np
from mdwb.convergence import moments, split_rhat, write_convergence
from mdwb.config import normalize_config
from mdwb.reporting import Series
from mdwb.cli import configure_convergence, configure_clustering
from mdwb.analysis import build_pooled_batches
from test_analysis import SELECTIONS


def series(values, replica='A', **kwargs):
    values = np.asarray(values, float)
    defaults = dict(replica=replica, analysis='rmsd', section='protein', name='value',
                    kind='timeseries', x=np.arange(len(values), dtype=float), values=values,
                    unit='angstrom', xunit='time (ns)', source='synthetic')
    defaults.update(kwargs)
    return Series(**defaults)


class ConvergenceTests(unittest.TestCase):
    def run_curves(self, data, **options):
        with tempfile.TemporaryDirectory() as temp:
            result = write_convergence({'convergence': options}, data, temp)
            tables = {}
            for p in (Path(temp)/'reports/convergence').glob('*.dat'):
                with p.open() as handle:
                    tables[p.stem] = list(csv.DictReader(handle, delimiter='\t'))
            return result, tables

    def test_stationary_and_offset_replicas(self):
        rng = np.random.default_rng(32)
        a, b = rng.normal(size=(2, 8000))
        self.assertLess(abs(split_rhat([a,b])-1), .01)
        self.assertGreater(split_rhat([a,b+4]), 2)
        self.assertIsNone(split_rhat([np.ones(100), np.ones(100)]))
        self.assertIsNone(split_rhat([a]))

    def test_correlated_stationary_signal_has_lower_effective_sample_size(self):
        from mdwb.diagnostics import sampling_diagnostic
        rng = np.random.default_rng(11)
        noise = rng.normal(size=16000)
        correlated = np.zeros(len(noise))
        for i in range(1, len(noise)):
            correlated[i] = .9 * correlated[i-1] + noise[i]
        independent = sampling_diagnostic(noise)[0]
        dependent = sampling_diagnostic(correlated)[0]
        self.assertGreater(independent['ess'], len(noise)*.7)
        self.assertLess(dependent['ess'], independent['ess']*.15)
        self.assertGreater(dependent['g'], 10)

    def test_same_means_can_hide_different_distributions(self):
        a = np.tile([-1, 1], 500)
        b = np.tile([-5, 5], 500)
        _, tables = self.run_curves([series(a), series(b, 'B')])
        self.assertEqual(a.mean(), b.mean())
        self.assertEqual(float(tables['between_replicas'][-1]['js_distance']), 1)

    def test_resource_limits_record_omissions(self):
        with tempfile.TemporaryDirectory() as temp:
            result = write_convergence({'diagnostics':{'max_samples':10}}, [series(range(100))], temp)
            self.assertEqual(result['within_rows'], 0)
            self.assertEqual(result['omissions'], 1)
            self.assertIn('resource_limit', (Path(temp)/'reports/convergence/omissions.dat').read_text())

    def test_drift_with_equal_whole_replica_means(self):
        a = np.linspace(-4, 4, 1000)
        self.assertGreater(split_rhat([a, a]), 1.5)
        _, t = self.run_curves([series(a)])
        rows = [r for r in t['within_replica'] if r['scope'] == 'cumulative']
        self.assertGreater(abs(float(rows[0]['mean_error'])), 3)
        self.assertAlmostEqual(float(rows[-1]['mean_error']), 0)
        self.assertFalse(t['between_replicas'])

    def test_equal_replica_reference_not_frame_weighted(self):
        _, t = self.run_curves([series(np.zeros(40)), series(np.ones(400)*10, 'B')])
        ref = t['references'][0]
        self.assertEqual(float(ref['mean']), 5)
        self.assertEqual(float(ref['sd']), 5)
        self.assertEqual(max(float(r['elapsed']) for r in t['between_replicas']), 39)
        self.assertTrue(all(float(r['js_distance']) == 1 for r in t['between_replicas']))

    def test_tail_reference_uses_physical_ns_not_frame_count(self):
        x = np.arange(100)*.5
        _, t = self.run_curves([series(np.arange(100), x=x)], reference_mode='tail', tail_ns=5)
        self.assertEqual(float(t['references'][0]['mean']), 94.5)
        _, t = self.run_curves([series(np.arange(10))], reference_mode='tail', tail_ns=20)
        self.assertFalse(t['within_replica'])
        self.assertIn('full_requested_duration', t['omissions'][0]['reason'])

    def test_custom_reference_is_metric_and_unit_specific(self):
        refs = {'rmsd/protein/value': {'mean': 2, 'sd': 1, 'unit': 'angstrom'}}
        _, t = self.run_curves([series(np.ones(100)*3)], reference_mode='custom', references=refs)
        self.assertEqual(float(t['within_replica'][-1]['mean_error']), 1)
        self.assertEqual(t['within_replica'][-1]['reference_js_distance'], '')
        refs['rmsd/protein/value']['unit'] = 'degree'
        _, t = self.run_curves([series(np.ones(100))], reference_mode='custom', references=refs)
        self.assertEqual(t['within_replica'][-1]['reference_mean'], '')
        self.assertTrue(t['omissions'])

    def test_file_reference_matches_custom(self):
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp)/'ref.json'
            path.write_text(json.dumps({'rmsd/protein/value': {'mean': 7, 'sd': 0, 'unit': 'angstrom'}}))
            _, t = self.run_curves([series(np.ones(100)*7)], reference_mode='file', reference_path=str(path))
            self.assertEqual(float(t['within_replica'][-1]['mean_error']), 0)

    def test_circular_wrap_and_categories_do_not_get_linear_ess(self):
        mean, sd = moments(np.array([-179,179]), True)
        self.assertAlmostEqual(abs(mean), 180)
        self.assertLess(sd, 2)
        _, t = self.run_curves([series(np.tile([-179,179],50), circular=True, unit='degree')])
        self.assertEqual(t['within_replica'][-1]['ess_stationary'], '')
        _, t = self.run_curves([series(np.tile([0,1],50), analysis='cluster', discrete=True, shared_basis=True),
                                series(np.tile([1,0],50), 'B', analysis='cluster', discrete=True, shared_basis=True)])
        self.assertEqual(t['within_replica'][-1]['mean'], '')
        self.assertEqual(float(t['between_replicas'][-1]['js_distance']), 0)
        self.assertTrue(t['state_populations'])

    def test_nonfinite_irregular_profiles_and_independent_clusters(self):
        for s in (series([1,2,np.nan,4]), series(range(5), x=np.array([0,1,2,4,5]))):
            _, t = self.run_curves([s])
            self.assertFalse(t['within_replica'])
            self.assertTrue(t['omissions'])
        _, t = self.run_curves([series(range(10), kind='profile')])
        self.assertFalse(t['within_replica'])
        _, t = self.run_curves([series(range(10), analysis='cluster', discrete=True),
                                series(range(10), 'B', analysis='cluster', discrete=True)])
        self.assertIn('not_shared', t['omissions'][0]['reason'])

    def test_metric_filter_and_disabled(self):
        _, t = self.run_curves([series(range(100))], metrics=['sasa'])
        self.assertFalse(t['within_replica'])
        result, _ = self.run_curves([series(range(100))], enabled=False)
        self.assertFalse(result['enabled'])

    def test_config_rejects_invalid_targets_and_unknown_time(self):
        base = {'topology': 'a.top', 'replicas': [{'name':'A', 'trajectories':['a.nc']}]}
        for options in ({'reference_mode':'tail'}, {'window_frames':True}, {'tail_ns':float('nan')},
                        {'reference_mode':'custom'}, {'reference_mode':'file'},
                        {'references':{'a/b/c':{'mean':0,'sd':-1,'unit':'angstrom'}}}):
            with self.subTest(options=options), self.assertRaises(ValueError):
                normalize_config({**base, 'convergence':options}, check_files=False)

    def test_wizard_tail_asks_for_time_interval(self):
        cfg = normalize_config({'topology':'a.top', 'replicas':[{'name':'A','trajectories':['a.nc']}]}, check_files=False)
        with patch('builtins.input', side_effect=['y','n','tail','2','5','rmsd','40','10']):
            configure_convergence(cfg)
        self.assertEqual(cfg['frames']['dt_ps'], 2)
        self.assertEqual(cfg['convergence']['tail_ns'], 5)
        self.assertEqual(cfg['convergence']['metrics'], ['rmsd'])

    def test_geometry_clustering_uses_original_masks_and_shared_labels(self):
        cfg = normalize_config({'topology':'a.top', 'replicas':[{'name':'A','trajectories':['a.nc']}],
             'analysis':{'enabled':['cluster','distance','dihedral']},
             'monitors':[{'name':'bond','type':'distance','masks':['@99','@100']},
                         {'name':'torsion','type':'dihedral','masks':['@99','@100','@101','@102']}],
             'advanced':{'cluster_metric':'features','cluster_algorithm':'hieragglo',
                         'cluster_features':[{'monitor':'bond','weight':2},{'monitor':'torsion','weight':.01}]}}, check_files=False)
        batch = build_pooled_batches(cfg, SELECTIONS)[0]
        self.assertTrue(batch['commands'][0].startswith('distance A_MON_bond'))
        self.assertIn('@99', batch['commands'][0])
        self.assertTrue(batch['commands'][2].startswith('strip'))
        self.assertIn('data A_MON_bond,A_MON_torsion euclid wgt 2,0.01', batch['post_commands'][0])
        self.assertIn('hieragglo clusters', batch['post_commands'][0])
        self.assertEqual(len(batch['artifacts']), 2)
        self.assertTrue(batch['artifacts'][0]['metadata']['pooled'])

    def test_report_integration_writes_html_tables_and_svg(self):
        from mdwb.reporting import report_results
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            artifacts = []
            for replica in ('A', 'B'):
                path = root / f'{replica}.dat'
                path.write_text('#Frame value\n' + ''.join(f'{i+1} {np.sin(i/5)}\n' for i in range(80)))
                artifacts.append({'name':replica, 'status':'complete', 'artifacts':[
                    {'path':str(path), 'analysis':'rmsd', 'section':'protein',
                     'kind':'timeseries', 'unit':'angstrom'}]})
            cfg = {'plots':{'enabled':True, 'dpi':72}, 'frames':{'dt_ps':1000, 'start':1, 'stride':1},
                   'convergence':{'checkpoints':4}, 'structure_view':{'enabled':False}}
            summary = report_results(cfg, {'replicas':artifacts}, root)
            self.assertGreater(summary['convergence']['within_rows'], 0)
            figures = [root / name for name in summary['figures'] if 'convergence/' in name]
            self.assertEqual(len(figures), 1)
            self.assertIn('Between-replica JS distance', figures[0].read_text())
            html = (root/'report_pages/sampling.html').read_text()
            self.assertIn('convergence/within_replica.dat', html)
            self.assertIn('Convergence over time', html)

    def test_longer_replica_gets_its_full_within_curve(self):
        _, tables = self.run_curves([series(np.arange(40)), series(np.arange(100), 'B')])
        rows = tables['within_replica']
        self.assertEqual(max(float(r['elapsed']) for r in rows if r['replica'] == 'B'), 99)
        self.assertEqual(max(float(r['elapsed']) for r in rows if r['replica'] == 'A'), 39)
        self.assertEqual(max(float(r['elapsed']) for r in tables['between_replicas']), 39)

    def test_structural_cluster_metrics(self):
        for metric in ('rms', 'dme', 'srmsd'):
            cfg = {'analysis':{'enabled':['cluster']}, 'advanced':{'cluster_metric':metric}}
            batch = build_pooled_batches(cfg, SELECTIONS)[0]
            self.assertIn(f'{metric} @* sieve', batch['post_commands'][0])


if __name__ == '__main__':
    unittest.main()
