"""Exercise an installed wheel and portable archive away from the source tree."""
import argparse
import json
from pathlib import Path
import subprocess
import sys
import tempfile


def main():
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--archive", required=True)
    args = cli.parse_args()
    archive = Path(args.archive).resolve()
    results = []
    with tempfile.TemporaryDirectory(prefix="md-bud-release-") as temp:
        root = Path(temp)
        environment_check = subprocess.run([sys.executable, "-c", "import mdbud; print(mdbud.__file__)"],
                                            cwd=root, text=True, capture_output=True, check=True)
        if "site-packages" not in environment_check.stdout:
            raise RuntimeError("Install the built wheel into this interpreter before testing.")
        for label, prefix in (("installed_wheel", [sys.executable, "-m", "mdbud"]),
                              ("portable_archive", [sys.executable, str(archive)])):
            for action in (["--version"], ["doctor"], ["demo", "--no-plots", "--output", str(root / label)]):
                process = subprocess.run([*prefix, *action], cwd=root, text=True, capture_output=True)
                if process.returncode:
                    raise RuntimeError(f"{label} {action}: {process.stdout}\n{process.stderr}")
            report = root / label / "index.html"
            if not report.is_file() or "__REPORT_DATA__" in report.read_text(encoding="utf-8"):
                raise RuntimeError(f"Report template/assets failed: {label}")
            failed = subprocess.run([*prefix, "validate", "absent.json"], cwd=root, capture_output=True)
            if failed.returncode != 2:
                raise RuntimeError(f"Failure exit code lost: {label}")
            results.append({"artifact": label, "version_doctor_demo_assets_failure_exit": "passed"})
    print(json.dumps({"installed_module": environment_check.stdout.strip(), "checks": results}, indent=2))


if __name__ == "__main__":
    main()
