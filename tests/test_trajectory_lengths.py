"""Read counts early; preserve joined-segment and physical-time trim semantics."""
import contextlib
import copy
import io
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch
from mdwb.config import DEFAULTS
from mdwb.cli import configure_frame_range
from mdwb.runner import read_replica_lengths, length_warning, preflight, RunError


class TrajectoryLengthTests(unittest.TestCase):
    def config(self):
        cfg = copy.deepcopy(DEFAULTS)
        cfg.update(topology='top.pdb', replicas=[{'name':'A','trajectories':['a1.nc','a2.nc']},
                                               {'name':'B','trajectories':['b.nc']}])
        return cfg

    def lengths(self, second=80):
        return [{'name':'A','trajectories':['a1.nc','a2.nc'],'lengths':[40,60]},
                {'name':'B','trajectories':['b.nc'],'lengths':[second]}]

    def test_counts_sum_segments_without_modifying_config(self):
        cfg = self.config()
        before = copy.deepcopy(cfg)
        with tempfile.TemporaryDirectory() as temp, patch('mdwb.runner.invoke', side_effect=['Frames: 40\n','Frames: 60\n','Frames: 80\n']) as invoke:
            replicas = read_replica_lengths('cpptraj', cfg, temp)
        self.assertEqual(cfg, before)
        self.assertEqual(replicas[0]['lengths'], [40,60])
        self.assertIn('A=100, B=80', length_warning(replicas))
        self.assertEqual(invoke.call_count, 3)
        self.assertTrue(all(c.args[1][-1] == '-tl' for c in invoke.call_args_list))
        self.assertIsNone(length_warning(self.lengths(100)))

    def test_invalid_count_is_not_claimed_equal(self):
        for log in ('Frames: 0', 'unknown reader'):
            with tempfile.TemporaryDirectory() as temp, patch('mdwb.runner.invoke', return_value=log), self.assertRaises(RunError):
                read_replica_lengths('cpptraj', self.config(), temp)

    def test_warning_precedes_range_prompt_and_can_trim_both_ends(self):
        cfg, out = self.config(), io.StringIO()
        def answer(prompt):
            if prompt.startswith('Saved-frame interval'):
                self.assertIn('TRAJECTORY LENGTH MISMATCH', out.getvalue())
            return next(answers)
        answers = iter(['1000','y','ns','10','60','2'])
        with patch('mdwb.runner.executable', return_value='cpptraj'), patch('mdwb.runner.read_replica_lengths', return_value=self.lengths()), patch('builtins.input', side_effect=answer), contextlib.redirect_stdout(out):
            configure_frame_range(cfg)
        self.assertEqual(cfg['frames'], {'start':11,'stop':61,'stride':2,'dt_ps':1000.0})
        self.assertIn('Analyze A: 26 selected frames', out.getvalue())
        self.assertIn('Analyze B: 26 selected frames', out.getvalue())

    def test_shortest_endpoint_default_and_frame_discard(self):
        cfg = self.config()
        with patch('mdwb.runner.executable', return_value='cpptraj'), patch('mdwb.runner.read_replica_lengths', return_value=self.lengths()), patch('builtins.input', side_effect=['0','y','11','','1']), contextlib.redirect_stdout(io.StringIO()):
            configure_frame_range(cfg)
        self.assertEqual(cfg['frames']['start'],11)
        self.assertEqual(cfg['frames']['stop'],80)

    def test_ns_bounds_round_inward_to_saved_frames(self):
        cfg = self.config()
        with patch('mdwb.runner.executable', return_value='cpptraj'), patch('mdwb.runner.read_replica_lengths', return_value=self.lengths(100)), patch('builtins.input', side_effect=['1000','ns','1.1','3.9','1']), contextlib.redirect_stdout(io.StringIO()):
            configure_frame_range(cfg)
        self.assertEqual((cfg['frames']['start'],cfg['frames']['stop']), (3,4))

    def test_missing_engine_reports_unverified_and_still_allows_selection(self):
        cfg, out = self.config(), io.StringIO()
        with patch('mdwb.runner.executable', side_effect=RunError('missing')), patch('builtins.input', side_effect=['0','5','90','3']), contextlib.redirect_stdout(out):
            configure_frame_range(cfg)
        self.assertIn('unverified', out.getvalue())
        self.assertNotIn('equal saved-frame', out.getvalue())
        self.assertEqual(cfg['frames']['stop'],90)

    def test_run_checks_lengths_before_any_atom_selection(self):
        with tempfile.TemporaryDirectory() as temp, patch('mdwb.runner.read_replica_lengths', side_effect=RunError('length probe first')), patch('mdwb.runner.invoke') as invoke:
            with self.assertRaisesRegex(RunError,'length probe first'):
                preflight('cpptraj',self.config(),{},[],[],Path(temp))
            invoke.assert_not_called()


if __name__ == '__main__':
    unittest.main()
