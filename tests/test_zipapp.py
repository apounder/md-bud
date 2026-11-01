"""Single-file cluster deployment must preserve assets and failure exit codes."""
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from zipfile import ZipFile


class ZipappTests(unittest.TestCase):
    def test_portable_archive_runs_outside_checkout_and_preserves_exit_codes(self):
        source = Path(__file__).resolve().parents[1]
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            archive = root / "md-bud.pyz"
            build = subprocess.run([sys.executable, str(source / "tools/build_zipapp.py"), "--output", str(archive)],
                                   cwd=root, capture_output=True, text=True)
            self.assertEqual(build.returncode, 0, build.stderr)
            result = subprocess.run([sys.executable, "-S", str(archive), "doctor"], cwd=root,
                                    text=True, capture_output=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            failure = subprocess.run([sys.executable, str(archive), "validate", "absent.json"], cwd=root,
                                     capture_output=True, text=True)
            self.assertEqual(failure.returncode, 2, failure.stdout + failure.stderr)
            with ZipFile(archive) as package:
                self.assertIn("mdwb/report_template.html", package.namelist())
                self.assertIn("mdwb/vendor/3Dmol-LICENSE.txt", package.namelist())
                self.assertIn("LICENSE", package.namelist())
