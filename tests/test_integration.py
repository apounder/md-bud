"""Campaign-to-analysis integration; fixtures are synthetic and never run AMBER."""
import contextlib
import copy
import io
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

from mdbud import campaign
from mdbud.cli import main
from mdbud.project import campaign_config
from mdwb.config import load_config
from test_config_topology import amber_metadata


class IntegrationTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.project = self.root / "campaign with spaces"
        self.top = self.root / "system.prmtop"
        self.coords = self.root / "initial.rst7"
        self.top.write_text(amber_metadata([("ALA", ["N", "CA", "C", "O"])] * 3, periodic=1))
        self.coords.write_text("SYNTHETIC restart placeholder\n")

    def cli(self, *args, expected=0):
        with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()) as errors:
            result = main(list(map(str, args)))
        self.assertEqual(result, expected, errors.getvalue())

    def generate(self, *args):
        self.cli("generate", "-i", self.project, "-p", self.top, "-s", self.coords,
                 "-A", "test", "-R", "2", "-n", "2", "-t", "0.02", "--ntwx", "2500", *args)

    def finish(self):
        manifest = campaign.load(self.project)
        for rep in manifest["replicas"]:
            directory = self.project / rep["name"]
            for stage in rep["stages"]:
                if stage["name"].startswith("md-"):
                    for output in stage["outputs"]:
                        (directory / output).write_text("SYNTHETIC output, deliberately not a real MD trajectory\n")
                    (directory / (stage["name"] + ".done")).write_text("test\n")

    def test_generation_and_handoff_preserve_replica_order_and_saved_time(self):
        self.generate("--combine")
        cfg = campaign_config(self.project, check_files=False)
        self.assertEqual(cfg["frames"]["dt_ps"], 5.0)
        self.assertEqual(cfg["analysis"]["level"], "basic")
        self.assertEqual([r["name"] for r in cfg["replicas"]], ["replica-001", "replica-002"])
        for rep in cfg["replicas"]:
            self.assertEqual([Path(p).name for p in rep["trajectories"]], ["md-0001.nc", "md-0002.nc"])
        self.assertEqual(len(cfg["notes"]["campaign_manifest_sha256"]), 64)

    def test_trim_and_stride_do_not_confuse_integration_dt_with_saved_frames(self):
        self.generate()
        cfg = campaign_config(self.project, check_files=False, discard_ns=0.01, stride=3)
        self.assertEqual(cfg["frames"], {"start": 3, "stop": -1, "stride": 3, "dt_ps": 5.0})
        with self.assertRaisesRegex(ValueError, "every expected"):
            campaign_config(self.project, check_files=False, discard_ns=0.04)

    def test_analysis_refuses_missing_or_unfinished_production_before_engine(self):
        self.generate()
        with patch("mdwb.runner.run_workflow") as run:
            self.cli("analyze", self.project, expected=2)
        run.assert_not_called()
        self.finish()
        (self.project / "replica-002/md-0002.done").unlink()
        with self.assertRaisesRegex(ValueError, "not complete"):
            campaign_config(self.project)

    def test_export_and_offline_plan_work_before_trajectories_exist(self):
        self.generate()
        config = self.project / "analysis_config.json"
        self.cli("analyze", self.project, "--save-config", config)
        exported = json.loads(config.read_text())
        self.assertFalse(Path(exported["topology"]).is_absolute())
        self.cli("analyze", self.project, "--plan", "--output", self.project / "preview")
        self.assertTrue((self.project / "preview/PLAN.md").is_file())
        self.assertFalse((self.project / "analysis_results").exists())
        self.cli("analyze", self.project, "--save-config", config, expected=2)
        self.cli("analyze", self.project, "--plan", "--output", self.project / "preview", expected=2)

    def test_exported_configuration_survives_campaign_relocation(self):
        self.generate()
        self.cli("analyze", self.project, "--save-config", self.project / "analysis.json")
        moved = self.root / "moved campaign"
        self.project.rename(moved)
        cfg = load_config(moved / "analysis.json", check_files=False)
        self.assertTrue(Path(cfg["topology"]).is_relative_to(moved))
        self.assertTrue(Path(cfg["output"]).is_relative_to(moved))
        self.assertTrue(all(Path(p).is_relative_to(moved) for r in cfg["replicas"] for p in r["trajectories"]))

    def test_changed_input_and_bad_segment_order_are_rejected(self):
        self.generate()
        target = self.project / "inputs/topology.prmtop"
        original = target.read_bytes()
        target.write_bytes(original + b"\n")
        with self.assertRaisesRegex(ValueError, "input changed"):
            campaign_config(self.project, check_files=False)
        target.write_bytes(original)
        manifest = campaign.load(self.project)
        stages = manifest["replicas"][0]["stages"]
        stages[-1], stages[-2] = stages[-2], stages[-1]
        campaign.save(self.project, manifest)
        with self.assertRaisesRegex(ValueError, "consecutive"):
            campaign_config(self.project, check_files=False)

    def test_prep_only_campaign_does_not_manufacture_production(self):
        self.generate("-f", "prep")
        self.cli("analyze", self.project, "--plan", expected=2)
        self.assertFalse((self.project / "analysis_results").exists())

    def test_template_keeps_custom_monitor_and_uses_campaign_identity(self):
        self.generate()
        template = {"topology": "wrong", "replicas": [], "frames": {"dt_ps": 1, "start": 2},
                    "analysis": {"enabled": ["distance"]},
                    "monitors": [{"name": "pair", "type": "distance", "masks": [":1@CA", ":2@CA"]}]}
        saved = copy.deepcopy(template)
        cfg = campaign_config(self.project, template=template, check_files=False)
        self.assertEqual(template, saved)
        self.assertEqual(cfg["monitors"][0]["name"], "pair")
        self.assertEqual(cfg["frames"]["dt_ps"], 5.0)
        self.assertEqual(cfg["frames"]["start"], 2)
        self.assertEqual(cfg["topology"], str(self.project / "inputs/topology.prmtop"))

    def test_profile_settings_are_validated_and_cli_overrides_them(self):
        profile = self.root / "cluster.json"
        self.cli("configure", profile, "--account", "profile-account", "--modules", "", "--mem", "8G")
        self.generate("--profile", profile, "--mem", "12G")
        settings = campaign.load(self.project)["config"]
        self.assertEqual(settings["account"], "test")
        self.assertEqual(settings["mem"], "12G")
        self.assertEqual(settings["modules"], [])
        self.cli("configure", profile, "--account", "other", expected=2)

    def test_profiles_cannot_submit_jobs_or_override_input_paths(self):
        profile = self.root / "bad.json"
        profile.write_text('{"submit":true}')
        self.cli("generate", "--profile", profile, expected=2)
        self.assertFalse(self.project.exists())

    def test_defaults_record_mc_barostat_and_unrestrained_equilibration(self):
        self.generate()
        rep = self.project / "replica-001"
        self.assertIn("barostat=2", (rep / "md-0001.in").read_text())
        self.assertIn("barostat=1", (rep / "eq-01.in").read_text())
        final = (rep / "eq-06.in").read_text()
        self.assertIn("ntr=0", final)
        self.assertIn("nstlim=50000", final)
        self.assertNotIn("restraintmask", final)
        self.assertIn("-c eq-06.rst7", (rep / "md-0001.sh").read_text())
        self.assertNotIn(b"\r\n", (rep / "md-0001.sh").read_bytes())

    def test_nonuniform_segment_sampling_rejected(self):
        self.cli("generate", "-i", self.project, "-p", self.top, "-s", self.coords,
                 "-A", "test", "--ns", "0.02", "--ntwx", "3000", expected=2)
        self.assertFalse(self.project.exists())

    def test_select_single_replica_and_reject_invalid_selection(self):
        self.generate()
        cfg = campaign_config(self.project, check_files=False, replica=2)
        self.assertEqual([r["name"] for r in cfg["replicas"]], ["replica-002"])
        self.cli("analyze", self.project, "-R", "3", "--plan", expected=2)

    def test_generation_requires_no_numerical_dependencies(self):
        code = "from mdbud.cli import main; raise SystemExit(main(['doctor']))"
        result = subprocess.run([sys.executable, "-S", "-c", code], cwd=Path(__file__).resolve().parents[1],
                                text=True, capture_output=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(json.loads(result.stdout)["missing"], [])

    def test_end_to_end_campaign_analysis_with_synthetic_cpptraj(self):
        from test_runner import SyntheticCpptraj
        self.generate()
        self.finish()
        with patch("mdwb.runner.executable", return_value="synthetic"), \
             patch("mdwb.runner.invoke", side_effect=SyntheticCpptraj(length_output="Frames: 4\n")):
            self.cli("analyze", self.project, "--no-plots", "--no-progress")
        report = self.project / "analysis_results"
        self.assertTrue((report / "index.html").is_file())
        manifest = json.loads((report / "manifest.json").read_text())
        self.assertEqual(manifest["status"], "complete")
        self.assertEqual([r["name"] for r in manifest["replicas"]], ["replica-001", "replica-002"])
        self.assertEqual([r["frame_count"] for r in manifest["replicas"]], [8, 8])
        self.assertEqual(len(manifest["inputs"]), 5)  # One topology and four segments.

    def test_incomplete_frame_counts_fail_even_when_file_markers_exist(self):
        from test_runner import SyntheticCpptraj
        self.generate()
        self.finish()
        with patch("mdwb.runner.executable", return_value="synthetic"), \
             patch("mdwb.runner.invoke", side_effect=SyntheticCpptraj(length_output="Frames: 3\n")):
            self.cli("analyze", self.project, "--no-plots", "--no-progress", expected=2)
        manifest = json.loads((self.project / "analysis_results/manifest.json").read_text())
        self.assertEqual(manifest["status"], "failed")

    def test_scheduler_chain_dry_run_does_not_submit(self):
        self.generate()
        with patch("mdbud.campaign.run") as external, contextlib.redirect_stdout(io.StringIO()) as out:
            campaign.submit(type("Args", (), {"project": self.project, "replica": 1, "dry_run": True})())
        external.assert_not_called()
        self.assertIn("--dependency=afterok:", out.getvalue())
        self.assertEqual(campaign.load(self.project)["jobs"], {})


if __name__ == "__main__":
    unittest.main()
