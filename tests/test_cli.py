"""Wizard configuration and offline planning exercise no native MD engine."""
import contextlib
import io
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from mdwb.cli import _indices, main, wizard
from mdwb.config import load_config, normalize_config
from test_config_topology import amber_metadata


class CliTests(unittest.TestCase):
    def test_wizard_standard_defaults_save_and_plan_without_cpptraj(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            topology = root / "protein.prmtop"
            topology.write_text(amber_metadata([("ALA", ["N", "CA", "C", "O"])] * 3), encoding="ascii")
            (root / "rep1.nc").write_bytes(b"Synthetic placeholder, not MD")
            config_path = root / "config.json"
            # Accept setup defaults, but explicitly defer execution.
            def answer(prompt):
                if prompt.startswith("Setup path:"):
                    return "guided"
                return "n" if prompt.startswith("Run CPPTRAJ now?") else ""
            with patch("builtins.input", side_effect=answer), contextlib.redirect_stdout(io.StringIO()):
                self.assertEqual(wizard(root, config_path), 0)
            config = json.loads(config_path.read_text())
            self.assertIsNone(config["analysis"]["enabled"])
            self.assertIsNone(config["frames"]["dt_ps"])
            self.assertEqual(config["analysis"]["level"], "standard")
            self.assertEqual(config["replicas"][0]["trajectories"], [str(root / "rep1.nc")])
            plan_dir = root / "plan"
            with contextlib.redirect_stdout(io.StringIO()):
                self.assertEqual(main(["plan", str(config_path), "--output", str(plan_dir)]), 0)
            plan = (plan_dir / "PLAN.md").read_text(encoding="utf-8")
            self.assertIn("rms", plan)
            self.assertIn("offline plan", plan)
            self.assertIn("secstruct", plan)

    def test_new_configuration_runs_by_default(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            (root / "protein.prmtop").write_text(amber_metadata([("ALA", ["N", "CA", "C", "O"])] * 3))
            (root / "rep1.nc").write_bytes(b"placeholder")
            target = root / "config.json"
            with patch("builtins.input", return_value=""), \
                 patch("mdwb.runner.run_workflow", return_value={"status": "complete"}) as run, \
                 contextlib.redirect_stdout(io.StringIO()):
                self.assertEqual(wizard(root, target), 0)
            run.assert_called_once()
            self.assertEqual(run.call_args.args[0], load_config(target))

    def test_existing_configuration_runs_without_setup_and_preserves_json(self):
        with tempfile.TemporaryDirectory() as temp:
            target = Path(temp) / "saved study.json"
            original = '{"saved": "unchanged"}\n'
            target.write_text(original)
            config = {"output": str(Path(temp) / "results")}
            for status, code in (("complete", 0), ("failed", 2)):
                with self.subTest(status=status), patch("builtins.input", return_value="") as prompt, \
                     patch("mdwb.cli.load_config", return_value=config) as load, \
                     patch("mdwb.cli.discover") as discover, \
                     patch("mdwb.runner.run_workflow", return_value={"status": status}) as run, \
                     contextlib.redirect_stdout(io.StringIO()):
                    self.assertEqual(main(["wizard", "--config", str(target)]), code)
                prompt.assert_called_once_with("Run CPPTRAJ now? (y/n) [y]: ")
                load.assert_called_once_with(target)
                run.assert_called_once_with(config)
                discover.assert_not_called()
                self.assertEqual(target.read_text(), original)

    def test_declining_existing_configuration_prints_portable_run_command(self):
        with tempfile.TemporaryDirectory() as temp:
            target = Path(temp) / "saved study.json"
            target.write_text("{}")
            output = io.StringIO()
            with patch("builtins.input", return_value="n"), \
                 patch("sys.argv", ["mdworkbench.pyz"]), \
                 patch("mdwb.cli.load_config") as load, \
                 patch("mdwb.runner.run_workflow") as run, contextlib.redirect_stdout(output):
                self.assertEqual(wizard(temp, target), 0)
            load.assert_not_called()
            run.assert_not_called()
            self.assertEqual(target.read_text(), "{}")
            self.assertIn(f"python mdworkbench.pyz run '{target}'", output.getvalue())
            self.assertIn("Run it later", output.getvalue())

    def test_default_startup_discovers_custom_configs_and_ignores_unrelated_json(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            for name in ("study1.json", "study2.json"):
                (root / name).write_text(json.dumps({"topology": "top.prmtop", "replicas": []}))
            (root / "manifest.json").write_text('{"status": "complete"}')
            (root / "broken.json").write_text("invalid JSON")
            with patch("pathlib.Path.cwd", return_value=root), \
                 patch("builtins.input", return_value="2"), \
                 patch("mdwb.cli.offer_run", return_value=0) as offer, \
                 patch("sys.stdin.isatty", return_value=True), contextlib.redirect_stdout(io.StringIO()):
                self.assertEqual(main([]), 0)
            offer.assert_called_once_with(root / "study2.json")

    def test_explicit_reordered_file_choices_remain_ordered(self):
        self.assertEqual(_indices("3,1-2", 3), [2, 0, 1])
        with self.assertRaises(ValueError):
            _indices("1,1-2", 3)

    def test_cli_missing_inputs_returns_error_code(self):
        with contextlib.redirect_stderr(io.StringIO()):
            self.assertEqual(main(["run", "does-not-exist-config.json"]), 2)

    def test_replica_names_cannot_collide_on_case_insensitive_filesystems(self):
        config = {"topology": "system.prmtop", "replicas": [
            {"name": "RepA", "trajectories": ["a.nc"]},
            {"name": "repa", "trajectories": ["b.nc"]}]}
        with self.assertRaisesRegex(ValueError, "Duplicate replica"):
            normalize_config(config, check_files=False)


if __name__ == "__main__":
    unittest.main()
