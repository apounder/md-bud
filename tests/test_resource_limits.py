"""Resource prompts use the actual joined, trimmed replica frame counts."""
import contextlib
import copy
import io
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from mdwb.cli import configure_resource_limits, wizard
from mdwb.config import DEFAULTS
from test_config_topology import amber_metadata


class ResourceLimitTests(unittest.TestCase):
    def config(self, analyses):
        config = copy.deepcopy(DEFAULTS)
        config['analysis']['enabled'] = analyses
        return config

    def prompt(self, config, replicas, answers):
        output = io.StringIO()
        with patch('builtins.input', side_effect=answers), contextlib.redirect_stdout(output):
            configure_resource_limits(config, replicas)
        return output.getvalue()

    def replicas(self):
        return [{'name': f'rep{i}', 'lengths': [20000]} for i in range(3)]

    def test_raise_retains_frames_and_separate_memory_guard(self):
        config = self.config(['pca'])
        output = self.prompt(config, self.replicas(), ['', '', ''])
        self.assertIn('60,000 frames across 3 replicas', output)
        self.assertEqual(config['advanced']['max_frames'], 60000)
        self.assertEqual(config['frames']['stride'], 1)
        self.assertEqual(config['advanced']['max_memory_gb'], 4)

    def test_cluster_limit_preserves_all_analyses_and_existing_stride(self):
        for stride, total in ((1, 60000), (3, 20001)):
            config = self.config(['rmsd', 'pca', 'cluster'])
            config['frames']['stride'] = stride
            before = copy.deepcopy(config['frames'])
            output = self.prompt(config, self.replicas(), ['20000', '', '', ''])
            self.assertIn('Please enter a valid value', output)
            self.assertEqual(config['frames'], before)
            self.assertEqual(config['advanced']['max_frames'], total)
            self.assertEqual(config['advanced']['cluster_sieve'], 10)

    def test_cluster_sieve_does_not_change_pca_or_input_windows(self):
        from mdwb.analysis import build_pooled_batches
        from mdwb.runner import select_windows
        config = self.config(['pca', 'cluster'])
        selections = {'representative': '@CA'}
        before = build_pooled_batches(config, selections)
        windows = select_windows([5, 15], config['frames'])
        config['advanced']['cluster_sieve'] = 20
        after = build_pooled_batches(config, selections)
        self.assertEqual(before[0], after[0])  # Full-frame PCA is identical.
        self.assertEqual(before[1]['commands'], after[1]['commands'])
        self.assertIn('sieve 20 random', str(after[1]))
        self.assertEqual(windows, select_windows([5, 15], config['frames']))
        self.assertEqual(sum(w['count'] for w in windows if w), 20)

    def test_trim_and_stride_continue_across_chunks(self):
        config = self.config(['pairwise_rmsd'])
        config['frames'].update(start=4, stop=15, stride=4)
        replicas = [{'name': 'a', 'lengths': [5, 15]}, {'name': 'b', 'lengths': [20]}]
        output = self.prompt(config, replicas, ['', '', ''])
        self.assertIn('6 frames across 2 replicas', output)

    def test_unverified_lengths_allow_explicit_limits(self):
        config = self.config(['pca'])
        output = self.prompt(config, None, ['60000', '3000', 'nan', '8'])
        self.assertIn('unverified', output)
        self.assertEqual(config['advanced']['max_frames'], 60000)
        self.assertEqual(config['advanced']['max_matrix_atoms'], 3000)
        self.assertEqual(config['advanced']['max_memory_gb'], 8)

    def test_basic_skips_prompts_and_dccm_only_asks_atoms_and_memory(self):
        self.prompt(self.config(['rmsd']), None, [])
        config = self.config(['dccm'])
        self.prompt(config, None, ['1500', '2'])
        self.assertEqual(config['advanced']['max_frames'], 20000)
        self.assertEqual(config['advanced']['max_matrix_atoms'], 1500)

    def test_both_wizard_paths_save_resource_choice(self):
        for quick in (True, False):
            with self.subTest(quick=quick), tempfile.TemporaryDirectory() as temp:
                root = Path(temp)
                (root/'top.prmtop').write_text(amber_metadata([('ALA', ['N', 'CA', 'C', 'O'])]*3))
                (root/'rep.nc').write_bytes(b'placeholder')
                target = root/'config.json'
                def answer(prompt):
                    if prompt.startswith('Run CPPTRAJ now?'):
                        return 'n'
                    if prompt.startswith('Setup path:'):
                        return 'basic' if quick else 'guided'
                    return ''
                def extras(config, info):
                    config['analysis']['enabled'] = ['pca']
                with patch('builtins.input', side_effect=answer), \
                     patch('mdwb.cli.configure_frame_range', return_value=self.replicas()), \
                     patch('mdwb.cli.additional_analyses', side_effect=extras), \
                     contextlib.redirect_stdout(io.StringIO()):
                    wizard(root, target)
                saved = json.loads(target.read_text())
                self.assertEqual(saved['advanced']['max_frames'], 60000)
                self.assertEqual(saved['frames']['stride'], 1)


if __name__ == '__main__':
    unittest.main()
