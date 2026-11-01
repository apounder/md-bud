"""Exercise the process boundary with CPPTRAJ-style argv re-tokenization."""
import contextlib
import io
import shlex
import subprocess
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch
from mdwb.runner import invoke, RunError
from mdwb.cli import configure_frame_range
from mdwb.config import DEFAULTS
import copy


class CpptrajPathTests(unittest.TestCase):
    def test_space_paths_survive_argv_and_trajectory_parsing(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            inputs = root / "good dye's inputs"
            inputs.mkdir()
            top, traj = inputs/'stripped rep.prmtop', inputs/'joined traj.nc'
            top.write_text('topology'); traj.write_text('trajectory')
            work = root/'results with spaces'; work.mkdir()
            def native(command, **kwargs):
                # CPPTRAJ first joins argv, then tokenizes, discarding OS argument boundaries.
                args = shlex.split(' '.join(command[1:]))
                t = Path(kwargs['cwd'])/args[args.index('-p')+1]
                self.assertEqual(t.read_text(),'topology')
                # TrajLength invokes AddInputTrajectory, which parses the value again.
                path = shlex.split(args[args.index('-y')+1])
                self.assertEqual(len(path),1)
                self.assertEqual((Path(kwargs['cwd'])/path[0]).read_text(),'trajectory')
                kwargs['stdout'].write('Frames: 123\n')
                return subprocess.CompletedProcess(command,0)
            with patch('mdwb.runner.subprocess.run',side_effect=native):
                result = invoke('cpptraj',['-p',str(top),'-y',str(traj),'-tl'],work,work/'length.log')
            self.assertIn('Frames: 123',result)
            self.assertIn(str(top),result)
            self.assertFalse(list(work.glob('.cpptraj-input-*')))
            self.assertEqual(traj.read_text(),'trajectory')

    def test_script_and_mask_paths_survive_without_changing_working_directory(self):
        with tempfile.TemporaryDirectory() as temp:
            work=Path(temp)/'analysis output'; work.mkdir()
            script=work/'input commands.in'; script.write_text('quit\n')
            top=work/'top file.prmtop'; top.write_text('topology')
            for args in (['-i',str(script)], ['-p',str(top),'-ms',"(:1@C4') | (:2@CA)"], ['-p',str(top),'-ms',":1@C4'"]):
                def native(command,**kwargs):
                    parsed=shlex.split(' '.join(command[1:]))
                    self.assertEqual(Path(kwargs['cwd']),work)
                    if '-i' in parsed:
                        self.assertEqual((work/parsed[1]).read_text(),'quit\n')
                    else:
                        self.assertEqual((work/parsed[1]).read_text(),'topology')
                        self.assertEqual(parsed[3],args[3])
                    return subprocess.CompletedProcess(command,0)
                with patch('mdwb.runner.subprocess.run',side_effect=native):
                    invoke('cpptraj',args,work,work/'process.log')

    def test_failure_retains_wizard_log_directory(self):
        cfg=copy.deepcopy(DEFAULTS)
        cfg.update(topology='top.pdb',replicas=[{'name':'A','trajectories':['a.nc']}])
        with tempfile.TemporaryDirectory() as temp:
            folder=Path(temp)/'lengths'; folder.mkdir()
            log=folder/'length_A_0001.log'
            def fail(*args):
                log.write_text('Error: reader failure\n')
                raise RunError(f'CPPTRAJ failed; see {log}')
            out=io.StringIO()
            with patch('mdwb.cli.tempfile.mkdtemp',return_value=str(folder)), patch('mdwb.runner.executable',return_value='cpptraj'), patch('mdwb.runner.read_replica_lengths',side_effect=fail), contextlib.redirect_stdout(out):
                with self.assertRaises(RunError):
                    configure_frame_range(cfg)
            self.assertTrue(log.is_file())
            self.assertIn(str(folder),out.getvalue())

    def test_native_failure_keeps_original_path_mapping(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp); top=root/'missing top.prmtop'
            def fail(command,**kwargs):
                kwargs['stdout'].write('Error: missing file\n')
                return subprocess.CompletedProcess(command,1)
            log=root/'failure.log'
            with patch('mdwb.runner.subprocess.run',side_effect=fail), self.assertRaises(RunError):
                invoke('cpptraj',['-p',str(top),'-tl'],root,log)
            self.assertIn(str(top),log.read_text())
            self.assertIn('Error: missing file',log.read_text())
            self.assertFalse(list(root.glob('.cpptraj-input-*')))


if __name__ == '__main__':
    unittest.main()
