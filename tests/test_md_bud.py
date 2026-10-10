"""Offline integration checks. No cluster jobs or AMBER simulations are run."""
import fcntl
import json
import os
from pathlib import Path
import shutil
import runpy
import subprocess
import tempfile
import unittest
from unittest.mock import patch

CLI = Path(__file__).resolve().parents[1] / "md-bud"
MOCK = r'''#!/usr/bin/env python3
import json, os, pathlib, sys
command = pathlib.Path(sys.argv[0]).name
store = pathlib.Path(os.environ['MOCK_STATE'])
state = json.loads(store.read_text())
args = sys.argv[1:]
if command == 'sbatch':
    if os.environ.get('FAIL_SBATCH'):
        print('submission rejected', file=sys.stderr); sys.exit(1)
    state['counter'] += 1
    job = str(state['counter'])
    state['jobs'][job] = 'PENDING'
    state['calls'].append(args)
    store.write_text(json.dumps(state))
    print('not-a-job-id' if os.environ.get('BAD_SBATCH') else job)
elif command == 'sacct':
    for job, status in state['jobs'].items(): print(job + '|' + status + '|')
elif command == 'squeue':
    for job, status in state['jobs'].items():
        if status in ['PENDING', 'RUNNING']: print(job + '|' + status)
elif command == 'scancel':
    for job in args: state['jobs'][job] = 'CANCELLED'
    store.write_text(json.dumps(state))
elif command == 'module': pass
elif command == 'srun': os.execvp(args[0], args)
elif command in ['pmemd.cuda', 'pmemd.MPI']:
    if os.environ.get('FAIL_ENGINE'): sys.exit(1)
    for flag in ['-o', '-r', '-x', '-inf']:
        if flag in args: pathlib.Path(args[args.index(flag)+1]).write_text('mock data\n')
elif command == 'cpptraj':
    if os.environ.get('FAIL_CPPTRAJ'): sys.exit(1)
    if '-tl' in args:
        print('CPPTRAJ 25.0: topology has 456 atoms')
        print('Frames: ' + os.environ.get('MOCK_FRAMES', '4'))
        print('Elapsed time: 0.123')
    else:
        pathlib.Path('combined.nc').write_text('combined data\n')
        print('combined')
'''


class CampaignTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix="md-bud-tests-")
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.project = self.root / "campaign with spaces"
        self.topology = self.root / "personal topology.prmtop"
        self.coordinates = self.root / "initial coords.rst7"
        self.topology.write_text("test topology\n")
        self.coordinates.write_text("test coordinates\n")
        self.store = self.root / "scheduler.json"
        self.store.write_text(json.dumps({"counter": 100, "jobs": {}, "calls": []}))
        commands = self.root / "bin"
        commands.mkdir()
        for name in ["sbatch", "sacct", "squeue", "scancel", "module", "srun",
                     "pmemd.cuda", "pmemd.MPI", "cpptraj"]:
            script = commands / name
            script.write_text(MOCK)
            script.chmod(0o755)
        self.env = os.environ | {"PATH": str(commands) + ":" + os.environ["PATH"],
                                 "MOCK_STATE": str(self.store)}
        self.env.pop("SLURM_JOB_ID", None)

    def cli(self, *args, expected=0, **env):
        result = subprocess.run([str(CLI), *map(str, args)], env=self.env | env,
                                cwd=self.root, capture_output=True, text=True)
        self.assertEqual(result.returncode, expected, result.stdout + result.stderr)
        return result

    def generate(self, *args):
        return self.cli("generate", "--output", self.project, "--topology", self.topology,
                        "--coordinates", self.coordinates, "--account", "test-account",
                        "--ns", "0.02", "--ntwx", "2500", "--segments", "2",
                        "--replicas", "2", "--combine", "--seed", "123", *args)

    def manifest(self):
        return json.loads((self.project / "md-bud.json").read_text())

    def engine(self, name, replica=1, expected=0, spool=False, **env):
        directory = self.project / f"replica-{replica:03d}"
        script = directory / f"{name}.sh"
        if spool:
            copy = self.root / "spooled-script.sh"
            shutil.copyfile(script, copy)
            script = copy
            env["SLURM_JOB_ID"] = "777"
        result = subprocess.run(["bash", str(script)], env=self.env | env, cwd=directory,
                                capture_output=True, text=True)
        self.assertEqual(result.returncode, expected, result.stdout + result.stderr)
        return result

    def test_help_is_read_only(self):
        self.cli("--help")
        self.cli("generate", "--help")
        self.assertFalse(self.project.exists())

    def test_generation_and_bash_syntax(self):
        self.generate()
        data = self.manifest()
        self.assertNotIn(str(self.topology), json.dumps(data))
        self.assertEqual(len(data["replicas"]), 2)
        self.assertEqual([stage["name"] for stage in data["replicas"][0]["stages"]],
                         ["min", "heat", "eq", "md-0001", "md-0002", "combine"])
        for script in self.project.glob("replica-*/*.sh"):
            subprocess.run(["bash", "-n", str(script)], check=True)
        first = (self.project / "replica-001/heat-01.in").read_text()
        second = (self.project / "replica-002/heat-01.in").read_text()
        self.assertNotEqual(first, second)
        self.assertNotIn("nmropt", (self.project / "replica-001/md-0001.in").read_text())

    def test_no_overwrite_and_bad_values(self):
        self.generate()
        self.cli("generate", "--output", self.project, "--topology", self.topology,
                 "--coordinates", self.coordinates, "--account", "test-account", expected=2)
        self.cli("generate", "--output", self.root / "invalid", "--topology", self.topology,
                 "--coordinates", self.coordinates, "--account", "bad\naccount", expected=2)
        self.cli("generate", "--output", self.root / "invalid", "--topology", self.topology,
                 "--coordinates", self.coordinates, "--account", "test", "--replicas", "0", expected=2)
        self.assertFalse((self.root / "invalid").exists())

    def test_workflow_choices(self):
        for workflow, names in [("prep", ["min", "heat", "eq"]),
                                ("production", ["md-0001"])]:
            destination = self.root / workflow
            self.cli("generate", "--output", destination, "--topology", self.topology,
                     "--coordinates", self.coordinates, "--account", "test", "--workflow", workflow)
            data = json.loads((destination / "md-bud.json").read_text())
            self.assertEqual([s["name"] for s in data["replicas"][0]["stages"]], names)

    def test_dependencies_and_duplicate_prevention(self):
        self.generate()
        self.cli("submit", self.project, "--dry-run")
        self.assertEqual(self.manifest()["jobs"], {})
        self.cli("submit", self.project)
        calls = json.loads(self.store.read_text())["calls"]
        self.assertEqual(len(calls), 12)
        self.assertFalse(any("--dependency" in arg for arg in calls[0]))
        self.assertIn("--dependency=afterok:101", calls[1])
        self.assertFalse(any("--dependency" in arg for arg in calls[6]))
        self.cli("submit", self.project)
        self.assertEqual(len(json.loads(self.store.read_text())["calls"]), 12)

    def test_rejected_and_ambiguous_submissions(self):
        self.generate()
        self.cli("submit", self.project, expected=2, FAIL_SBATCH="1")
        self.assertEqual(self.manifest()["jobs"], {})
        self.cli("submit", self.project, expected=2, BAD_SBATCH="1")
        self.assertIsNone(self.manifest()["jobs"]["replica-001/min"]["id"])
        self.cli("submit", self.project, expected=2)
        self.cli("reset", self.project, "--from", "min", "--apply", expected=2)

    def test_failure_exit_and_spool_directory(self):
        self.generate()
        self.engine("min", expected=1, FAIL_ENGINE="1")
        self.assertFalse((self.project / "replica-001/min.done").exists())
        for stage in ["min", "heat", "eq", "md-0001", "md-0002", "combine"]:
            self.engine(stage, spool=True)
            self.assertTrue((self.project / f"replica-001/{stage}.done").exists())
        self.engine("md-0001", expected=1)

    def test_reset_preserves_inputs_previous_stages_other_replica(self):
        self.generate()
        for stage in ["min", "heat", "eq", "md-0001", "md-0002", "combine"]:
            self.engine(stage)
        replica = self.project / "replica-001"
        other = self.project / "replica-002/md-0002.nc"
        other.write_text("other replica")
        before = (replica / "md-0002.nc").read_bytes()
        self.cli("reset", self.project, "--replica", "1", "--from", "md-2")
        self.assertEqual((replica / "md-0002.nc").read_bytes(), before)
        self.cli("reset", self.project, "--replica", "1", "--from", "md-2", "--apply")
        self.assertTrue((replica / "md-0001.nc").exists())
        self.assertTrue((replica / "eq.done").exists())
        self.assertFalse((replica / "md-0002.nc").exists())
        self.assertFalse((replica / "combined.nc").exists())
        self.assertTrue((replica / "md-0002.in").exists())
        self.assertTrue(other.exists())
        self.assertTrue((self.project / "inputs/topology.prmtop").exists())
        archived = list((self.project / "archive").glob("*/replica-001/md-0002.nc"))
        self.assertEqual(len(archived), 1)
        self.assertEqual(archived[0].read_bytes(), before)

    def test_cancel_before_reset_and_resubmit(self):
        self.generate()
        self.cli("submit", self.project)
        self.cli("reset", self.project, "--replica", "1", "--from", "eq", "--apply", expected=2)
        self.cli("reset", self.project, "--replica", "1", "--from", "eq", "--apply", "--cancel")
        self.assertNotIn("replica-001/eq", self.manifest()["jobs"])
        self.assertIn("replica-001/heat", self.manifest()["jobs"])
        self.assertIn("replica-002/eq", self.manifest()["jobs"])
        self.cli("submit", self.project, "--replica", "1")
        calls = json.loads(self.store.read_text())["calls"]
        self.assertIn("--dependency=afterok:102", calls[12])

    def test_reset_refuses_unknown_scheduler_state(self):
        self.generate()
        self.cli("submit", self.project)
        state = json.loads(self.store.read_text())
        state["jobs"] = {}
        self.store.write_text(json.dumps(state))
        self.cli("reset", self.project, "--from", "min", "--apply", "--cancel", expected=2)
        self.assertFalse((self.project / "archive").exists())

    def test_manual_stage_lock_prevents_reset(self):
        self.generate()
        with (self.project / "replica-001/.min.lock").open("a") as lock:
            fcntl.flock(lock, fcntl.LOCK_EX)
            self.cli("reset", self.project, "--replica", "1", "--from", "min", "--apply", expected=2)
        self.assertFalse((self.project / "archive").exists())

    def test_frame_parser_and_failures_do_not_move_data(self):
        self.generate()
        for replica in [1, 2]:
            for segment in [1, 2]:
                (self.project / f"replica-{replica:03d}/md-{segment:04d}.nc").write_text("trajectory")
        result = self.cli("frames", self.project)
        self.assertIn("Total: 16 frames", result.stdout)
        self.cli("frames", self.project, expected=1, MOCK_FRAMES="3")
        self.cli("frames", self.project, expected=2, FAIL_CPPTRAJ="1")
        self.assertEqual(len(list(self.project.glob("replica-*/*.nc"))), 4)
        self.assertFalse((self.project / "sum.txt").exists())
        (self.project / "replica-001/md-0001.nc").unlink()
        self.cli("frames", self.project, expected=1)

    def test_inputs_cannot_silently_change(self):
        self.generate()
        (self.project / "inputs/topology.prmtop").write_text("different topology")
        self.cli("submit", self.project, expected=2)
        self.assertEqual(self.manifest()["jobs"], {})

    def test_partial_archive_failure_restores_outputs(self):
        self.generate()
        replica = self.project / "replica-001"
        (replica / "md-0002.nc").write_text("trajectory")
        (replica / "md-0002.out").write_text("output")
        module = runpy.run_path(str(CLI))
        manifest = self.manifest()
        original = Path.rename

        def simulated_failure(source, target):
            if source == replica / "md-0002.out":
                raise OSError("simulated archive failure")
            return original(source, target)

        with patch.object(Path, "rename", simulated_failure):
            with self.assertRaises(OSError):
                module["archive_outputs"](self.project, manifest, ["replica-001/md-0002"])
        self.assertEqual((replica / "md-0002.nc").read_text(), "trajectory")
        self.assertEqual((replica / "md-0002.out").read_text(), "output")
        self.assertEqual(self.manifest()["archives"], [])

    def test_failed_job_requires_reset(self):
        self.generate()
        self.cli("submit", self.project)
        state = json.loads(self.store.read_text())
        state["jobs"]["101"] = "FAILED"
        self.store.write_text(json.dumps(state))
        self.cli("submit", self.project, expected=2)
        self.assertEqual(len(json.loads(self.store.read_text())["calls"]), 12)

    def test_optional_restraints_and_nvt(self):
        restraints = self.root / "distance restraints.dist"
        restraints.write_text("&rst iat=1,2, r1=1, r2=2, r3=3, r4=4, rk2=1, rk3=1, /\n")
        self.generate("--restraints", restraints, "--ensemble", "nvt", "--modules", "")
        content = (self.project / "replica-001/md-0001.in").read_text()
        self.assertIn("nmropt=1", content)
        self.assertIn("DISANG=../inputs/restraints.dist", content)
        self.assertIn("ntp=0", content)
        self.assertNotIn("pres0", content)
        self.assertNotIn("module load", (self.project / "replica-001/md-0001.sh").read_text())

    def test_inconsistent_downstream_completion_blocks_submission(self):
        self.generate()
        replica = self.project / "replica-001"
        for extension in ["nc", "rst7", "out", "done"]:
            (replica / f"md-0002.{extension}").write_text("stale downstream output")
        self.cli("submit", self.project, expected=2)
        self.assertEqual(self.manifest()["jobs"], {})


if __name__ == "__main__":
    unittest.main()
