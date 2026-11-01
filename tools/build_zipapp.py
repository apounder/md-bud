"""Build a portable MD-bud executable archive using only the standard library."""
import argparse
from pathlib import Path
import tempfile
import shutil
import zipapp


def main():
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--output", default="dist/md-bud.pyz")
    args = cli.parse_args()
    source = Path(__file__).resolve().parents[1]
    target = Path(args.output).resolve()
    if target.exists():
        cli.error("Output exists; choose a new filename.")
    target.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory() as temp:
        root = Path(temp)
        for package in ("mdbud", "mdwb"):
            shutil.copytree(source / package, root / package,
                            ignore=shutil.ignore_patterns("__pycache__", "*.pyc"))
        for name in ("LICENSE", "THIRD_PARTY.md", "CITATION.md"):
            shutil.copyfile(source / name, root / name)
        (root / "__main__.py").write_text("from mdbud.cli import main\nraise SystemExit(main())\n", encoding="utf-8")
        zipapp.create_archive(root, target, compressed=True,
                              interpreter="/usr/bin/env python3")
    print(target)


if __name__ == "__main__":
    main()
