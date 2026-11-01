"""Rebuild the offline example with NumPy; no CPPTRAJ, Matplotlib or network use."""
from pathlib import Path
import argparse
import json
import shutil
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from mdwb.demo import create_demo
from mdwb.dashboard import write_dashboard


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=ROOT / "examples" / "rendered")
    args = parser.parse_args()
    output = args.output.resolve()
    summary = create_demo(output, render=False)
    config = json.loads((output / "synthetic_config.json").read_text(encoding="utf-8"))
    manifest = json.loads((output / "synthetic_manifest.json").read_text(encoding="utf-8"))
    manifest["structure_view"] = {
        "status": "example", "selection": "Experimental crambin, PDB 1CRN",
        "message": "Viewer example: experimental crambin (PDB 1CRN). The synthetic plots are unrelated to this structure; no crambin simulation was performed."
    }
    (output / "reference").mkdir(exist_ok=True)
    shutil.copyfile(ROOT / "examples" / "1CRN.pdb", output / "reference" / "viewer.pdb")
    # Store portable example paths; report regeneration resolves these against output.
    replacements = []
    for replica in manifest["replicas"]:
        for artifact in replica["artifacts"]:
            original = artifact["path"]
            artifact["path"] = Path(original).relative_to(output).as_posix()
            replacements.append((original, artifact["path"]))
    summary["combined_data"] = "reports/all_replicas.dat"
    (output / "reports" / "report_summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    (output / "synthetic_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    write_dashboard(config, manifest, summary, output)
    # Remove machine-specific paths from the distributable demonstration only.
    for path in output.rglob("*"):
        if path.is_file() and path.suffix in {".html", ".dat", ".json", ".md"}:
            text = path.read_text(encoding="utf-8")
            for original, relative in replacements:
                text = text.replace(original, relative).replace(json.dumps(original)[1:-1], relative)
            path.write_text(text, encoding="utf-8")
    print(f"Open {output / 'index.html'} in a browser with JavaScript and WebGL.")


if __name__ == "__main__":
    main()
