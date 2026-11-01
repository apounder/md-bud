import tempfile
import unittest
from pathlib import Path

from mdbud.qc import inspect_log


class QCTests(unittest.TestCase):
    def read(self, text):
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / "md.out"
            path.write_text(text)
            return inspect_log(path)

    def test_summary_excludes_averages_and_keeps_units(self):
        result = self.read(" NSTEP = 100 TEMP(K) = 300.0 PRESS = 1.0\n Etot = -1.0D+02 EPtot = -150.0\n VOLUME = 1000.0\n Density = 1.0\n NSTEP = 200 TEMP(K) = 310.0\n A V E R A G E S   O V E R\n TEMP(K) = 999.0\n")
        self.assertEqual(result["samples"], 2)
        self.assertEqual(result["metrics"]["temperature_K"]["mean"], 305.0)
        self.assertEqual(result["metrics"]["total_energy_kcal_mol"]["mean"], -100.0)
        self.assertEqual(result["metrics"]["density_g_cm3"]["last"], 1.0)
        self.assertEqual(result["errors"], [])

    def test_numeric_and_engine_failures_are_reported(self):
        result = self.read(" NSTEP = 10 TEMP(K) = NaN\n Etot = ********\n SHAKE cannot be accomplished\n")
        self.assertEqual(len(result["errors"]), 3)
        self.assertNotIn("temperature_K", result["metrics"])

    def test_banners_and_empty_logs_do_not_create_measurements(self):
        result = self.read("******************************\n no sampled results\n")
        self.assertEqual(result, {"samples": 0, "metrics": {}, "errors": []})
