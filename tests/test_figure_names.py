"""Readable exports retain distinct identities and report links."""
import copy
import tempfile
import unittest
from pathlib import Path

from mdwb.config import DEFAULTS
from mdwb.figure_names import metric_figure_stem
from mdwb.reporting import report_results


class FigureNameTests(unittest.TestCase):
    def test_names_describe_metric_and_preserve_collisions(self):
        first = metric_figure_stem('rg', 'active site', 'A_RG_site')
        self.assertIn('radius-of-gyration__active-site', first)
        self.assertNotEqual(first, metric_figure_stem('rg', 'active-site', 'A_RG_site'))
        self.assertNotEqual(first, metric_figure_stem('rg', 'active site', 'A_RG_other'))
        self.assertNotIn('/', metric_figure_stem('rg', '../site', 'a/b'))
        self.assertLessEqual(len(metric_figure_stem('rg', 'a'*1000, 'b'*1000)), 110)

    def test_rendered_diagnostics_have_readable_names_and_lag_help(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            config = copy.deepcopy(DEFAULTS)
            config['frames'].update(dt_ps=10, stride=5)
            config['plots']['font'] = 'DejaVu Sans'
            replicas = []
            for name in ('rep01', 'rep02'):
                path = root/f'{name}.dat'
                path.write_text('#Frame A_RG_site\n' + ''.join(f'{i+1} {2+(i%9)/10}\n' for i in range(100)))
                replicas.append({'name': name, 'status': 'complete', 'artifacts': [
                    {'path': path.name, 'analysis': 'rg', 'section': 'active_site',
                     'kind': 'timeseries', 'unit': 'angstrom'}]})
            result = report_results(config, {'replicas': replicas}, root)
            figures = result['figures']
            for kind in ('autocorrelation', 'convergence-over-time', 'between-replica-distribution-distance',
                         'probability-distribution', 'distribution-and-running-mean', 'replica-comparison'):
                matches = [p for p in figures if p.endswith(f'__{kind}.svg')]
                self.assertTrue(matches, kind)
                for path in matches:
                    self.assertIn('radius-of-gyration__active-site', path)
                    self.assertTrue((root/path).is_file())
                    self.assertEqual(result['figure_analyses'][path], 'rg')
            pages = '\n'.join(p.read_text() for p in root.rglob('*.html'))
            self.assertIn('Lag is the separation', pages)
            self.assertIn('50 ps (0.05 ns)', pages)
            self.assertIn('does not prove', pages)
            self.assertFalse(any('rendering failed' in w for w in result['warnings']))
