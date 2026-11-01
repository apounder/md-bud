"""Meaningful reporting tests using small synthetic cpptraj ASCII fixtures."""
import csv
import importlib.util
import tempfile
import unittest
from pathlib import Path

import numpy as np

from mdwb.reporting import (
    NumericTable, Series, TableFormatError, _coordinate, _matrix_from_table,
    _nastruct_series, circular_statistics, concatenate_series,
    descriptive_statistics, parse_numeric_table, report_results,
)


class ReportingTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.config = {"plots": {"enabled": False}, "stats": {"block_size": 2}, "frames": {"dt_ps": None, "start": 1, "stride": 1}}

    def tearDown(self):
        self.temp.cleanup()

    def write(self, name, content):
        path = self.root / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(content, encoding="utf-8")
        return path

    def read_tsv(self, path):
        with path.open(newline="", encoding="utf-8") as handle:
            return list(csv.DictReader(handle, delimiter="\t"))

    def artifact(self, name, kind="timeseries", **extra):
        return {"path": name, "analysis": "rmsd", "section": "protein", "kind": kind, "unit": "angstrom", **extra}

    def test_numeric_headers_fortran_nan_and_duplicate_columns(self):
        path = self.write("data.dat", "# description of this file\n#Frame value value\n1 2D+0 nan\n2 3d-1 8\n")
        table = parse_numeric_table(path)
        self.assertEqual(table.columns, ["Frame", "value", "value_2"])
        self.assertAlmostEqual(table.data[1, 1], 0.3)
        self.assertTrue(np.isnan(table.data[0, 2]))

    def test_ragged_or_text_rows_are_not_silently_lost(self):
        for content in ("#Frame A\n1 2\n2 3 4\n", "#Frame A\n1 2\n2 absent\n"):
            with self.assertRaises(TableFormatError):
                parse_numeric_table(self.write("bad.dat", content))

    def test_time_converts_saved_frame_interval_with_stride_and_start(self):
        table = NumericTable(np.array([[1, 4], [2, 5]], float), ["Frame", "rmsd"], [])
        x, unit, index = _coordinate(table, {}, {"frames": {"dt_ps": 10, "start": 11, "stride": 5}}, "timeseries")
        np.testing.assert_allclose(x, [0.1, 0.15])
        self.assertEqual(unit, "time (ns)")
        self.assertEqual(index, 0)
        x, unit, _ = _coordinate(table, {}, self.config, "timeseries")
        np.testing.assert_array_equal(x, [1, 2])
        self.assertEqual(unit, "analyzed frame")

    def test_unknown_coordinate_is_rejected_in_multicolumn_timeseries(self):
        table = NumericTable(np.ones((3, 2)), ["something", "value"], [])
        with self.assertRaises(TableFormatError):
            _coordinate(table, {}, self.config, "timeseries")

    def test_concatenation_marks_boundaries_without_cross_replica_steps(self):
        a = Series("A", "rmsd", "protein", "r", "timeseries", np.array([4, 6, 8]), np.array([1, 2, 3]), "angstrom", "ns", "a")
        b = Series("B", "rmsd", "protein", "r", "timeseries", np.array([4, 6]), np.array([10, 11]), "angstrom", "ns", "b")
        x, y, boundaries = concatenate_series([a, b])
        np.testing.assert_array_equal(x, [0, 2, 4, 6, 8])
        np.testing.assert_array_equal(y, [1, 2, 3, 10, 11])
        self.assertEqual(boundaries, [5])

    def test_block_means_preserve_missing_sample_locations(self):
        stats = descriptive_statistics(np.array([1, 1, np.nan, 9, 5, 5, 99]), 2)
        self.assertEqual(stats["n_complete_blocks"], 2)
        self.assertAlmostEqual(stats["block_mean_sd"], np.std([1, 5], ddof=1))
        self.assertEqual(stats["n_finite"], 6)
        self.assertNotIn("standard_error", stats)

    def test_unequal_lengths_equal_replica_weight_and_missing_replica(self):
        self.write("a.dat", "#Frame RMSD\n1 0\n2 0\n3 0\n4 0\n")
        self.write("b.dat", "#Frame RMSD\n1 10\n")
        manifest = {"replicas": [{"name": "A", "status": "complete", "artifacts": [self.artifact("a.dat")]}, {"name": "B", "status": "complete", "artifacts": [self.artifact("b.dat")]}, {"name": "C", "status": "failed", "artifacts": [self.artifact("missing.dat")]}]}
        summary = report_results(self.config, manifest, self.root)
        rows = self.read_tsv(self.root / "reports" / "replicate_statistics.dat")
        self.assertEqual(float(rows[0]["equal_replica_mean"]), 5.0)
        self.assertEqual(rows[0]["n_replicas"], "2")
        self.assertEqual(summary["replicas_with_numeric_data"], 2)
        tidy = self.read_tsv(self.root / "reports" / "all_replicas.dat")
        self.assertEqual(len(tidy), 5)
        self.assertEqual({row["replica"] for row in tidy}, {"A", "B"})
        self.assertTrue(any("missing" in warning for warning in summary["warnings"]))
        wide = self.read_tsv(next((self.root / "reports" / "combined").glob("*.dat")))
        self.assertEqual(wide[1]["B"], "")

    def test_one_replica_does_not_claim_between_replica_sd(self):
        self.write("a.dat", "#Frame RMSD\n1 1\n2 3\n")
        report_results(self.config, {"replicas": [{"name": "A", "artifacts": [self.artifact("a.dat")]}]}, self.root)
        row = self.read_tsv(self.root / "reports" / "replicate_statistics.dat")[0]
        self.assertEqual(row["between_replica_mean_sd"], "")

    def test_matrix_xyz_maps_axes_and_rejects_duplicate_cells(self):
        table = NumericTable(np.array([[2, 10, 0.2], [1, 20, 0.3], [1, 10, 1], [2, 20, 1]]), ["x", "y", "value"], [])
        matrix = _matrix_from_table("A", {"metadata": {"matrix_format": "xyz"}}, table, "source")
        np.testing.assert_array_equal(matrix.x, [1, 2])
        np.testing.assert_array_equal(matrix.y, [10, 20])
        np.testing.assert_allclose(matrix.values, [[1, 0.2], [0.3, 1]])
        table.data[-1] = table.data[0]
        with self.assertRaises(TableFormatError):
            _matrix_from_table("A", {"metadata": {"matrix_format": "xyz"}}, table, "source")

    def test_cpptraj_indexed_matrix_preserves_noncontiguous_atom_axes(self):
        table = parse_numeric_table(self.write("dccm.dat", "#Atom-Atom 1 2\n1 1 0.5\n2 0.5 1\n"))
        matrix = _matrix_from_table("A", {"analysis": "dccm", "metadata": {"atom_ids": [12, 97]}}, table, "source")
        np.testing.assert_array_equal(matrix.x, [12, 97])
        np.testing.assert_array_equal(matrix.y, [12, 97])
        self.assertEqual(matrix.matrix_type, "correlation")
        self.assertEqual(matrix.xunit, "topology atom number")
        table = parse_numeric_table(self.write("contact.dat", "#Residue-Residue 2 7 9\n3 1 0 2\n8 0 1 3\n"))
        matrix = _matrix_from_table("A", {}, table, "source")
        np.testing.assert_array_equal(matrix.x, [2, 7, 9])
        np.testing.assert_array_equal(matrix.y, [3, 8])
        self.assertEqual(matrix.values.shape, (2, 3))

    def test_matrix_means_require_matching_axes(self):
        self.write("a.dat", "1 2\n3 4\n")
        self.write("b.dat", "5 6\n7 8\n")
        artifacts = [self.artifact(name, "matrix", metadata={"matrix_format": "dense", "x_labels": labels}) for name, labels in (("a.dat", [1, 2]), ("b.dat", [1, 3]))]
        manifest = {"replicas": [{"name": str(i), "artifacts": [artifact]} for i, artifact in enumerate(artifacts)]}
        summary = report_results(self.config, manifest, self.root)
        self.assertTrue(any("axes or shapes differ" in warning for warning in summary["warnings"]))
        self.assertFalse((self.root / "reports" / "matrices").exists())
        artifacts[1]["metadata"]["x_labels"] = [1, 2]
        report_results(self.config, manifest, self.root)
        rows = self.read_tsv(next((self.root / "reports" / "matrices").glob("*.dat")))
        self.assertEqual(float(rows[0]["equal_replica_mean"]), 3.0)
        self.assertAlmostEqual(float(rows[0]["between_replica_sd"]), np.sqrt(8), places=8)

    def test_circular_means_do_not_average_179_and_minus179_to_zero(self):
        stats = circular_statistics(np.array([179.0, -179.0]))
        self.assertAlmostEqual(abs(stats["circular_mean_degree"]), 180)
        self.assertGreater(stats["resultant_length"], 0.99)
        undefined = circular_statistics(np.array([0.0, 180.0]))
        self.assertIsNone(undefined["circular_mean_degree"])
        self.write("torsion.dat", "#Frame phi\n1 179\n2 -179\n")
        manifest = {"replicas": [{"name": "A", "artifacts": [self.artifact("torsion.dat", unit="degree", metadata={"circular": True})]}]}
        report_results(self.config, manifest, self.root)
        self.assertEqual(self.read_tsv(self.root / "reports" / "descriptive_statistics.dat"), [])
        row = self.read_tsv(self.root / "reports" / "circular_statistics.dat")[0]
        self.assertAlmostEqual(abs(float(row["circular_mean_degree"])), 180)

    def test_dssp_states_are_heatmap_and_fractions_not_numeric_means(self):
        self.write("dssp.dat", "#Frame DSSP_12 DSSP_13\n1 4 0\n2 4 2\n3 0 2\n")
        manifest = {"replicas": [{"name": "A", "artifacts": [self.artifact("dssp.dat", analysis="dssp", unit="state", metadata={"categorical": True, "palette": "dssp"})]}]}
        summary = report_results(self.config, manifest, self.root)
        self.assertEqual(summary["matrix_count"], 1)
        self.assertEqual(self.read_tsv(self.root / "reports" / "descriptive_statistics.dat"), [])
        rows = self.read_tsv(self.root / "reports" / "categorical_fractions.dat")
        self.assertEqual(len(rows), 4)
        self.assertAlmostEqual(sum(float(row["fraction"]) for row in rows if row["series"] == "DSSP_12"), 1)

    def test_nastruct_pair_identity_units_missing_groove(self):
        path = self.write("BPstep.na.dat", "#Frame BP1 BP2 Shift Slide Rise Tilt Roll Twist Zp Major Minor\n1 1-24 2-23 0 1 3.4 2 3 34 1 ---- 5\n1 2-23 3-22 0 2 3.3 1 2 35 1 12 6\n2 1-24 2-23 0 1 3.5 2 3 33 1 11 5\n")
        series = _nastruct_series(path, "A", self.artifact(path.name, "table", analysis="nastruct", metadata={"format": "nastruct"}), self.config)
        rise = next(item for item in series if item.name == "Rise[1-24/2-23]")
        np.testing.assert_allclose(rise.values, [3.4, 3.5])
        self.assertEqual(rise.unit, "angstrom")
        twist = next(item for item in series if item.name == "Twist[1-24/2-23]")
        self.assertTrue(twist.circular)
        major = next(item for item in series if item.name == "Major[1-24/2-23]")
        self.assertTrue(np.isnan(major.values[0]))

    def test_unknown_text_table_is_retained_and_indexed(self):
        path = self.write("hbond.dat", "#Acceptor DonorH Donor Frames Frac\n:1@O :2@H :2@N 90 0.9\n")
        result = report_results(self.config, {"replicas": [{"name": "A", "artifacts": [self.artifact(path.name, "table")]}]}, self.root)
        self.assertTrue(path.exists())
        self.assertEqual(result["status"], "complete_with_warnings")
        row = self.read_tsv(self.root / "reports" / "artifact_index.dat")[0]
        self.assertEqual(row["parse_status"], "unparsed")

    def test_known_hbond_table_preserves_contact_identity_and_column_units(self):
        path = self.write("hbond.dat", "#Acceptor DonorH Donor Frames Frac AvgDist AvgAng\n:1@O :2@H :2@N 90 0.9 2.8 155\n")
        result = report_results(self.config, {"replicas": [{"name": "A", "artifacts": [self.artifact(path.name, "table", metadata={"format": "hbond_occupancy"})]}]}, self.root)
        self.assertEqual(result["status"], "complete")
        rows = self.read_tsv(self.root / "reports" / "all_replicas.dat")
        self.assertEqual(len(rows), 4)
        distance = next(row for row in rows if row["series"].endswith("|AvgDist"))
        self.assertIn("Acceptor=:1@O", distance["series"])
        self.assertEqual(distance["unit"], "angstrom")
        fraction = next(row for row in rows if row["series"].endswith("|Frac"))
        self.assertEqual(fraction["unit"], "fraction")

    def test_pooled_matrices_reported_without_counting_as_replica(self):
        self.write("covariance.dat", "#Coord-Coord 1 2\n1 1 -0.2\n2 -0.2 2\n")
        self.write("modes.dat", "Eigenvector block format retained intentionally\n1 3.4\n0.2 0.3 0.4\n")
        manifest = {"replicas": [], "pooled": [{"name": "pca_pooled", "status": "complete", "artifacts": [self.artifact("covariance.dat", "matrix", analysis="pca_covariance"), self.artifact("modes.dat", "table", metadata={"format": "eigenvectors"})]}]}
        result = report_results(self.config, manifest, self.root)
        self.assertEqual(result["status"], "complete")
        self.assertEqual(result["matrix_count"], 1)
        self.assertEqual(result["replicas_with_numeric_data"], 0)
        self.assertEqual(result["pooled_groups_requested"], 1)
        rows = self.read_tsv(self.root / "reports" / "all_replicas.dat")
        self.assertEqual({row["replica"] for row in rows}, {"POOLED:pca_pooled"})
        indexed = self.read_tsv(self.root / "reports" / "artifact_index.dat")
        self.assertEqual(indexed[1]["parse_status"], "retained_native_format")

    def test_missing_matplotlib_is_explicitly_unavailable(self):
        from unittest.mock import patch
        self.config["plots"]["enabled"] = True
        with patch("mdwb.reporting._plotting", side_effect=ImportError("no matplotlib")):
            result = report_results(self.config, {"replicas": []}, self.root)
        self.assertEqual(result["plots_status"], "unavailable")
        self.assertEqual(result["status"], "complete_with_warnings")

    def test_demo_is_explicitly_synthetic_and_has_unequal_replicas(self):
        from mdwb.demo import create_demo
        import json
        result = create_demo(self.root, render=False)
        self.assertTrue(result["synthetic"])
        self.assertEqual(result["replicas_with_numeric_data"], 3)
        self.assertEqual(result["matrix_count"], 3)
        self.assertIn("SYNTHETIC", (self.root / "SYNTHETIC_DEMO_README.md").read_text())
        manifest = json.loads((self.root / "synthetic_manifest.json").read_text())
        self.assertEqual({replica["frame_count"] for replica in manifest["replicas"]}, {240, 300, 360})
        unrelated = self.root / "unrelated"
        unrelated.mkdir()
        (unrelated / "file.dat").write_text("existing data")
        with self.assertRaises(ValueError):
            create_demo(unrelated, render=False)

    @unittest.skipUnless(importlib.util.find_spec("matplotlib"), "Matplotlib is not installed; no packages are installed by this test")
    def test_svg_smoke_and_appended_warning(self):
        self.config["plots"]["enabled"] = True
        self.write("a.dat", "#Frame RMSD\n1 1\n2 2\n3 3\n")
        self.write("b.dat", "#Frame RMSD\n1 2\n2 3\n")
        result = report_results(self.config, {"replicas": [{"name": name, "artifacts": [self.artifact(f"{name.lower()}.dat")]} for name in ("A", "B")]}, self.root)
        self.assertGreaterEqual(len(result["figures"]), 5)
        svg = next((self.root / "figures" / "concatenated").glob("*.svg")).read_text(encoding="utf-8")
        self.assertIn("not a continuous trajectory", svg)
        self.assertIn("<text", svg)


if __name__ == "__main__":
    unittest.main()
