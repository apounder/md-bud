"""Workflow tests using synthetic CPPTRAJ responses, never a chemistry engine.

The smoke tests exercise real configuration, topology inspection, input-file
generation, manifest writing and numerical reporting. Only the subprocess
boundary is replaced. Their success is not evidence of native CPPTRAJ execution.
"""
import csv
import copy
import json
import shlex
import subprocess
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import numpy as np

from mdwb import runner
from mdwb.config import DEFAULTS, normalize_config


def write_pdb(path):
    lines, index = [], 0
    for residue in range(1, 4):
        for atom, element in (("N", "N"), ("CA", "C"), ("C", "C"), ("O", "O")):
            index += 1
            lines.append(f"ATOM  {index:5d} {atom:^4s} ALA A{residue:4d}    "
                         f"{float(index):8.3f}{float(residue):8.3f}{0.0:8.3f}"
                         f"{1.0:6.2f}{0.0:6.2f}          {element:>2s}\n")
    path.write_text("".join(lines) + "END\n", encoding="utf-8")


def make_config(root, *, enabled=None, two_reps=True):
    top = root / "topology.pdb"
    write_pdb(top)
    paths = [root / name for name in ("rep1_part1.nc", "rep1_part2.nc", "rep2.nc")]
    for path in paths:
        path.write_bytes(b"Synthetic fixture: deliberately not a real NetCDF trajectory.\n")
    replicas = [{"name": "rep1", "trajectories": [str(path) for path in paths[:2]]}]
    if two_reps:
        replicas.append({"name": "rep2", "trajectories": [str(paths[2])]} )
    data = {"topology": str(top), "output": str(root / "results"),
            "replicas": replicas, "fit_mask": ":1-3@CA",
            "sections": [{"name": "custom_backbone", "mask": ":1-3@N,CA,C,O"}],
            "imaging": {"mode": "off"}, "analysis": {"level": "basic"},
            "frames": {"start": 2, "stop": -1, "stride": 2, "dt_ps": 5.0}}
    if enabled is not None:
        data["analysis"]["enabled"] = enabled
    return normalize_config(data, root)


class SyntheticCpptraj:
    """Write tiny deterministic tables as a controllable external-process double."""
    def __init__(self, *, fail_batch=None, omit_alignment=False, mask_output=None, length_output=None):
        self.calls = []
        self.scripts = []
        self.fail_batch = fail_batch
        self.omit_alignment = omit_alignment
        self.mask_output = mask_output
        self.length_output = length_output

    def __call__(self, binary, args, cwd, log_path):
        cwd, log_path = Path(cwd), Path(log_path)
        self.calls.append(list(args))
        log_path.parent.mkdir(parents=True, exist_ok=True)
        if "--version" in args:
            response = "SYNTHETIC_CPPTRAJ_TEST_DOUBLE 1.0\n"
        elif "-ms" in args:
            response = self.mask_output if self.mask_output is not None else "Selected= 1 2 3 4\n"
        elif "-mr" in args:
            response = "Selected= 1 2 3\n"
        elif "-tl" in args:
            trajectory = Path(args[args.index("-y") + 1]).name
            length = {"rep1_part1.nc": 5, "rep1_part2.nc": 4, "rep2.nc": 7}.get(trajectory, 8)
            response = self.length_output if self.length_output is not None else f"Reader diagnostic\nFrames: {length}\n"
        elif "-i" in args:
            script = Path(args[args.index("-i") + 1]).read_text(encoding="utf-8")
            self.scripts.append((cwd, script))
            parsed = [shlex.split(line, posix=True) for line in script.splitlines()
                      if line.strip() and not line.startswith("#")]
            nframes = sum((int(tokens[3]) - int(tokens[2])) // int(tokens[4]) + 1
                          for tokens in parsed if tokens[0] == "trajin" and len(tokens) == 5)
            for tokens in parsed:
                if tokens[0] == "trajout":
                    (cwd / tokens[1]).write_text("SYNTHETIC COMMON REFERENCE\n", encoding="utf-8")
                elif "out" in tokens:
                    target = cwd / tokens[tokens.index("out") + 1]
                    if target.name == "alignment_check.dat" and self.omit_alignment:
                        continue
                    is_profile = tokens[0] == "atomicfluct"
                    nrows = 3 if is_profile else nframes
                    offset = 2.0 if "rep2" in cwd.parts else 0.0
                    values = np.arange(1, nrows + 1, dtype=float) * 0.25 + offset
                    if cwd.name == self.fail_batch:
                        values[:] = 9999.0
                    np.savetxt(target, np.column_stack([np.arange(1, nrows + 1), values]),
                               fmt="%.8g", header=f"{'Res' if is_profile else 'Frame'} {tokens[1]}")
                elif tokens[0] == "writedata":
                    target = cwd / tokens[1]
                    np.savetxt(target, np.column_stack([np.arange(1, nframes + 1), np.ones(nframes)]),
                               fmt="%.8g", header="Frame DEFERRED")
            if cwd.name == self.fail_batch:
                log_path.write_text("Error: SYNTHETIC intentional batch failure\n", encoding="utf-8")
                raise runner.RunError(f"SYNTHETIC batch failure; see {log_path}")
            response = f"SYNTHETIC wrote {nframes} selected frames\n"
        else:
            raise AssertionError(f"Unexpected mocked CPPTRAJ arguments: {args}")
        log_path.write_text(response, encoding="utf-8")
        return response


class FrameWindowTests(unittest.TestCase):
    def test_stride_continues_across_segment_boundaries(self):
        windows = runner.select_windows([5, 4, 6], {"start": 4, "stop": 14, "stride": 3})
        self.assertEqual([(w["start"], w["stop"], w["count"]) for w in windows],
                         [(4, 4, 1), (2, 2, 1), (1, 4, 2)])

    def test_many_chunkings_equal_selection_on_joined_replica(self):
        for lengths in ([3, 2, 4], [1, 1, 7], [9], [4, 5]):
            for start in range(1, 10):
                for stride in range(1, 6):
                    for stop in (-1, 9):
                        windows = runner.select_windows(lengths, {"start": start, "stop": stop, "stride": stride})
                        actual = [w["offset"] + frame for w in windows if w
                                  for frame in range(w["start"], w["stop"] + 1, w["stride"])]
                        self.assertEqual(actual, list(range(start, 10, stride)))

    def test_unselected_segments_and_single_frame_boundary(self):
        windows = runner.select_windows([5, 4, 6], {"start": 6, "stop": 6, "stride": 1})
        self.assertIsNone(windows[0])
        self.assertEqual(windows[1]["start"], 1)
        self.assertEqual(windows[1]["count"], 1)
        self.assertIsNone(windows[2])

    def test_short_replica_ranges_fail_before_subprocess_actions(self):
        for frames in ({"start": 10, "stop": -1, "stride": 1}, {"start": 1, "stop": 10, "stride": 1}):
            with self.assertRaises(runner.RunError):
                runner.select_windows([5, 4], frames)


class PreflightTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.config = make_config(self.root)

    def _preflight(self, fake):
        with patch.object(runner, "invoke", side_effect=fake):
            return runner.preflight("synthetic", self.config, {}, [], [], self.root / "preflight")

    def test_parses_frame_counts_and_retains_replica_boundaries(self):
        selected, replicas = self._preflight(SyntheticCpptraj())
        self.assertEqual(selected[self.config["fit_mask"]], [1, 2, 3, 4])
        self.assertEqual([r["lengths"] for r in replicas], [[5, 4], [7]])
        self.assertEqual([r["frame_count"] for r in replicas], [4, 3])
        self.assertNotIn("windows", self.config["replicas"][0])

    def test_empty_or_unrecognized_masks_are_rejected(self):
        for output in ("Selected=\n", "Some banner without selected indices\n"):
            with self.subTest(output=output), self.assertRaises(runner.RunError):
                self._preflight(SyntheticCpptraj(mask_output=output))

    def test_too_few_fit_atoms_are_rejected(self):
        with self.assertRaisesRegex(runner.RunError, "at least three"):
            self._preflight(SyntheticCpptraj(mask_output="Selected= 1 2\n"))

    def test_unknown_or_empty_trajectory_length_is_rejected(self):
        for output in ("Frames: 0\n", "Frame count unknown\n", "Frames: -1\n"):
            with self.subTest(output=output), self.assertRaisesRegex(runner.RunError, "frame count"):
                self._preflight(SyntheticCpptraj(length_output=output))

    def test_overlapping_interface_selections_are_rejected(self):
        self.config["interactions"] = [{"name": "overlap", "mask1": ":1", "mask2": ":1-2"}]
        with self.assertRaisesRegex(runner.RunError, "overlap"):
            self._preflight(SyntheticCpptraj())

    def test_contact_matrix_cap_uses_residue_span_not_selected_count(self):
        fake = SyntheticCpptraj()
        def residue_gap(binary, args, cwd, log_path):
            if "-mr" in args:
                return "Selected= 2 9000\n"
            return fake(binary, args, cwd, log_path)
        batch = {"name": "sparse_interface", "contact_masks": [":2", ":9000"]}
        with patch.object(runner, "invoke", side_effect=residue_gap):
            with self.assertRaisesRegex(runner.RunError, "residue span 8999"):
                runner.preflight("synthetic", self.config, {}, [batch], [], self.root / "preflight")


class ResourceTests(unittest.TestCase):
    def setUp(self):
        self.config = copy.deepcopy(DEFAULTS)

    def test_pca_requires_enough_frames_for_requested_modes(self):
        batch = {"name": "pca_pooled", "coordinate_mask": ":1-10", "pooled": True}
        with self.assertRaisesRegex(runner.RunError, "more frames"):
            runner.resource_check(self.config, batch, 10, 3)

    def test_cluster_sieving_cannot_leave_fewer_frames_than_clusters(self):
        batch = {"name": "cluster_pooled", "coordinate_mask": ":1-10", "pooled": True}
        with self.assertRaisesRegex(runner.RunError, "Too few sieved frames"):
            runner.resource_check(self.config, batch, 10, 20)

    def test_pairwise_respects_frame_and_memory_limits(self):
        batch = {"name": "pairwise_pooled", "coordinate_mask": ":1-10", "pooled": True}
        with self.assertRaisesRegex(runner.RunError, "max_frames"):
            runner.resource_check(self.config, batch, 10, self.config["advanced"]["max_frames"] + 1)
        self.config["advanced"]["max_memory_gb"] = 0.001
        with self.assertRaisesRegex(runner.RunError, "memory estimate"):
            runner.resource_check(self.config, batch, 10, 1000)

    def test_atom_matrix_cap_precedes_allocation(self):
        batch = {"name": "dccm", "coordinate_mask": ":1-5000"}
        with self.assertRaisesRegex(runner.RunError, "max_matrix_atoms"):
            runner.resource_check(self.config, batch, 5000, 100)


class ProcessBoundaryTests(unittest.TestCase):
    def test_nonzero_exit_preserves_subprocess_log(self):
        with tempfile.TemporaryDirectory() as temp:
            log = Path(temp) / "failure.log"
            source = Path(temp) / "a path/input.in"
            source.parent.mkdir()
            source.write_text("SYNTHETIC input\n")
            def failed_process(command, **kwargs):
                self.assertEqual(command[:2], ["cpptraj", "-i"])
                self.assertEqual((Path(kwargs["cwd"]) / command[2]).read_text(), source.read_text())
                self.assertNotIn("shell", kwargs)
                kwargs["stdout"].write("Details before failure\nError: invalid atom mask\n")
                return subprocess.CompletedProcess(command, 2)
            with patch.object(runner.subprocess, "run", side_effect=failed_process):
                with self.assertRaisesRegex(runner.RunError, "exit 2"):
                    runner.invoke("cpptraj", ["-i", "a path/input.in"], temp, log)
            self.assertIn("invalid atom mask", log.read_text())

    def test_native_error_is_detected_even_with_zero_exit_code(self):
        with tempfile.TemporaryDirectory() as temp:
            def failed_process(command, **kwargs):
                kwargs["stdout"].write("  Error: setup failed\n")
                return subprocess.CompletedProcess(command, 0)
            with patch.object(runner.subprocess, "run", side_effect=failed_process):
                with self.assertRaisesRegex(runner.RunError, "setup failed"):
                    runner.invoke("cpptraj", [], temp, Path(temp) / "failure.log")


class SyntheticWorkflowTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)

    def run_mocked(self, config, fake):
        with patch.object(runner, "executable", return_value="synthetic-cpptraj"), \
             patch.object(runner, "invoke", side_effect=fake):
            return runner.run_workflow(config, make_plots=False, progress=lambda message: None)

    def test_basic_workflow_generates_manifest_and_joint_numeric_tables(self):
        config = make_config(self.root)
        fake = SyntheticCpptraj()
        manifest = self.run_mocked(config, fake)
        output = Path(config["output"])
        self.assertEqual(manifest["status"], "complete")
        self.assertEqual(manifest["report"]["replicas_with_numeric_data"], 2)
        self.assertEqual(manifest["report"]["figures"], [])
        self.assertTrue(all(rep["status"] == "complete" for rep in manifest["replicas"]))
        self.assertIn("sha256", manifest["inputs"][0])
        self.assertEqual(json.loads((output / "manifest.json").read_text())["status"], "complete")
        with (output / "reports" / "all_replicas.dat").open() as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual({row["replica"] for row in rows}, {"rep1", "rep2"})
        self.assertEqual({row["analysis"] for row in rows}, {"rmsd_global", "rmsd_local", "rg", "rmsf"})
        rep1_map = np.loadtxt(output / "replicas" / "rep1" / "frame_map.dat")
        np.testing.assert_array_equal(rep1_map[:, 3], [2, 4, 6, 8])
        np.testing.assert_array_equal(rep1_map[:, 4], [0, 10, 20, 30])
        references = [line for _, script in fake.scripts for line in script.splitlines() if line.startswith("reference ")]
        self.assertEqual(len(references), 2)
        self.assertEqual(len(set(references)), 1)
        self.assertTrue(references[0].endswith("[COMMON]"))
        prep = [script for cwd, script in fake.scripts
                if cwd.name == "reference" and "trajout common.rst7 restart" in script]
        self.assertEqual(len(prep), 1)
        self.assertIn('" 2 2\n', prep[0])
        self.assertNotIn("windows", config["replicas"][0])

    def test_failed_batch_leftovers_are_excluded_from_combined_data(self):
        config = make_config(self.root, enabled=["rg", "sasa"])
        manifest = self.run_mocked(config, SyntheticCpptraj(fail_batch="sasa"))
        self.assertEqual(manifest["status"], "partial")
        self.assertTrue(all(rep["status"] == "partial" for rep in manifest["replicas"]))
        self.assertTrue(all({art["analysis"] for art in rep["artifacts"]} == {"rg"} for rep in manifest["replicas"]))
        output = Path(config["output"])
        self.assertTrue(list(output.glob("replicas/*/sasa/sasa*.dat")))
        combined = (output / "reports" / "all_replicas.dat").read_text()
        self.assertNotIn("9999", combined)
        self.assertNotIn("\tsasa\t", combined)
        self.assertIn("SYNTHETIC intentional batch failure", (output / "replicas" / "rep1" / "sasa" / "cpptraj.log").read_text())

    def test_missing_alignment_check_marks_batch_partial(self):
        config = make_config(self.root, enabled=["rg"])
        manifest = self.run_mocked(config, SyntheticCpptraj(omit_alignment=True))
        self.assertEqual(manifest["status"], "partial")
        self.assertTrue(all(not rep["artifacts"] for rep in manifest["replicas"]))

    def test_preflight_failure_persists_failed_manifest(self):
        config = make_config(self.root)
        with self.assertRaises(runner.RunError):
            self.run_mocked(config, SyntheticCpptraj(mask_output="Selected=\n"))
        manifest = json.loads((Path(config["output"]) / "manifest.json").read_text())
        self.assertEqual(manifest["status"], "failed")
        self.assertIn("zero atoms", manifest["error"])
        self.assertTrue(list(Path(config["output"]).glob("preflight/mask*.log")))

    def test_existing_results_are_never_overwritten(self):
        config = make_config(self.root)
        output = Path(config["output"])
        output.mkdir()
        sentinel = output / "keep.txt"
        sentinel.write_text("existing result")
        fake = SyntheticCpptraj()
        with self.assertRaisesRegex(runner.RunError, "not empty"):
            self.run_mocked(config, fake)
        self.assertEqual(sentinel.read_text(), "existing result")
        self.assertEqual(fake.calls, [])

    def test_post_commands_are_written_after_run_before_quit(self):
        config = make_config(self.root, enabled=["rg"], two_reps=False)
        artifact = {"path": "deferred.dat", "kind": "timeseries", "analysis": "deferred", "section": "test", "unit": "count", "metadata": {}}
        batch = {"name": "deferred", "description": "Synthetic deferred output", "commands": ["createcrd TEST"],
                 "post_commands": ["writedata deferred.dat TEST"], "artifacts": [artifact]}
        fake = SyntheticCpptraj()
        with patch.object(runner, "build_batches", return_value=[batch]):
            manifest = self.run_mocked(config, fake)
        self.assertEqual(manifest["status"], "complete")
        script = next(script for cwd, script in fake.scripts if cwd.name == "deferred")
        lines = script.splitlines()
        self.assertLess(lines.index("createcrd TEST"), lines.index("run"))
        self.assertLess(lines.index("run"), lines.index("writedata deferred.dat TEST"))
        self.assertLess(lines.index("writedata deferred.dat TEST"), lines.index("quit"))

    def test_pooled_split_uses_verified_counts_and_local_frame_numbers(self):
        output = self.root
        source = output / "pooled.dat"
        source.write_text("#Frame A_PC:1\n1 10\n2 20\n3 30\n4 40\n5 50\n")
        artifact = {"path": str(source), "kind": "timeseries", "analysis": "pca", "section": "representative", "unit": "angstrom", "metadata": {"pooled": True}}
        replicas = [{"name": "rep1", "frame_count": 2, "artifacts": []}, {"name": "rep2", "frame_count": 3, "artifacts": []}]
        runner._split_pooled({"name": "pca_pooled"}, [artifact], replicas, output)
        split = np.loadtxt(replicas[1]["artifacts"][0]["path"])
        np.testing.assert_array_equal(split, [[1, 30], [2, 40], [3, 50]])
        self.assertTrue(replicas[1]["artifacts"][0]["metadata"]["shared_basis"])
        self.assertIn("A_PC:1", Path(replicas[1]["artifacts"][0]["path"]).read_text())

    def test_pooled_split_rejects_missing_or_reordered_frames(self):
        for content in ("#Frame PC\n1 5\n3 6\n", "#Frame PC\n2 5\n1 6\n"):
            source = self.root / "pooled.dat"
            source.write_text(content)
            artifact = {"path": str(source), "kind": "timeseries"}
            replicas = [{"name": "rep1", "frame_count": 2, "artifacts": []}]
            with self.assertRaisesRegex(runner.RunError, "safely split"):
                runner._split_pooled({"name": "pca_pooled"}, [artifact], replicas, self.root)
            self.assertEqual(replicas[0]["artifacts"], [])


if __name__ == "__main__":
    unittest.main()
