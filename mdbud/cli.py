"""One public CLI for both engines; heavy dependencies load only on demand."""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import shutil
import sys

from . import __version__

CAMPAIGN = {"generate", "submit", "status", "reset", "frames"}
ANALYSIS = {"wizard", "discover", "validate", "catalog", "plan", "run", "report", "demo", "bundle"}


def parser():
    cli = argparse.ArgumentParser(prog="md-bud", description="MD-bud: prepare AMBER campaigns, manage Slurm jobs, analyze replicas and share reports.",
                                  epilog="Start: md-bud generate --help | md-bud analyze --help | md-bud wizard --quick")
    cli.add_argument("--version", action="version", version=f"MD-bud {__version__}")
    commands = cli.add_subparsers(dest="command")
    descriptions = {"generate": "Generate AMBER inputs and Slurm jobs (--profile FILE supported)",
                    "submit": "Submit or preview a campaign dependency chain",
                    "status": "Inspect a campaign (--offline works without Slurm)",
                    "reset": "Preview or archive a stage and downstream outputs",
                    "frames": "Check production frame counts with CPPTRAJ",
                    "analyze": "Analyze a campaign directly, or export/preview its analysis",
                    "qc": "Summarize production temperature/energy/density logs and numerical failures",
                    "configure": "Save reusable cluster settings as a local JSON profile",
                    "doctor": "Check prerequisites for generation, analysis or cluster execution",
                    "wizard": "Guided analysis setup for existing trajectories",
                    "discover": "Find existing topology and trajectory files",
                    "validate": "Validate an existing analysis JSON",
                    "catalog": "List available scientific analyses",
                    "plan": "Preview an analysis JSON without running CPPTRAJ",
                    "run": "Run an existing analysis JSON",
                    "report": "Regenerate reports from analysis outputs",
                    "demo": "Create a clearly labelled synthetic report",
                    "bundle": "Package analysis outputs for sharing"}
    for name, description in descriptions.items():
        commands.add_parser(name, help=description)
    return cli


def generate(argv):
    from . import campaign
    pre = argparse.ArgumentParser(add_help=False)
    pre.add_argument("--profile", default=os.environ.get("MDBUD_PROFILE"))
    args, remainder = pre.parse_known_args(argv)
    if "--help" in remainder or "-h" in remainder:
        print("Reusable settings: --profile FILE (or MDBUD_PROFILE). Explicit flags override the profile.\n")
        return campaign.main(["generate", *remainder])
    options = []
    if args.profile:
        data = json.loads(Path(args.profile).expanduser().read_text(encoding="utf-8-sig"))
        if not isinstance(data, dict):
            raise ValueError("Profile must be a JSON object of generate option names and values.")
        gen = campaign.parser()._subparsers._group_actions[0].choices["generate"]
        actions = {a.dest: a for a in gen._actions if a.option_strings and a.dest not in
                   {"help", "output", "topology", "coordinates", "restraints", "submit"}}
        for key, value in data.items():
            if key not in actions:
                raise ValueError(f"Unsupported profile setting: {key}")
            action = actions[key]
            flag = next(s for s in action.option_strings if s.startswith("--"))
            if isinstance(action, argparse._StoreTrueAction):
                if type(value) is not bool:
                    raise ValueError(f"Profile {key} must be true or false.")
                if value:
                    options.append(flag)
            elif type(value) in (str, int, float):
                options.append(f"{flag}={value}")
            else:
                raise ValueError(f"Profile {key} must be a string or number.")
    return campaign.main(["generate", *options, *remainder])


def configure(argv):
    from .project import write_new_json
    cli = argparse.ArgumentParser(prog="md-bud configure", description="Save cluster settings; existing files are never overwritten.")
    cli.add_argument("file", help="Local JSON profile, e.g. cluster.json")
    for flag in ("account", "modules", "gpu", "mem", "prep-mem", "partition", "mail", "walltime"):
        cli.add_argument("--" + flag)
    args = vars(cli.parse_args(argv))
    path = args.pop("file")
    data = {k: v for k, v in args.items() if v is not None}
    if not data:
        cli.error("Supply at least one setting, e.g. --account your-allocation")
    print(f"Profile saved: {write_new_json(path, data)}")
    print(f"Use: md-bud generate --profile {path} -i runs/study -p system.prmtop -s system.inpcrd")
    return 0


def analyze(argv):
    cli = argparse.ArgumentParser(prog="md-bud analyze", description="Analyze all ordered production segments in a campaign. Defaults: basic RMSD/RMSF/Rg, all frames, per-replica identities retained.")
    cli.add_argument("project", help="Campaign directory containing md-bud.json")
    mode = cli.add_mutually_exclusive_group()
    mode.add_argument("--plan", action="store_true", help="Write an offline plan; trajectories may not exist yet")
    mode.add_argument("--save-config", metavar="FILE", help="Export a movable analysis JSON without running; refuse overwrite")
    cli.add_argument("--config", help="Analysis JSON template for custom selections/monitors; campaign files and dt_ps take precedence")
    cli.add_argument("--level", choices=["basic", "standard", "advanced"], help="Preset; overrides template analysis choices (default: basic)")
    cli.add_argument("--output", help="New results/plan directory (default: PROJECT/analysis_results)")
    from .campaign import positive
    cli.add_argument("-R", "--replica", type=positive, help="Analyze one replica (default: all)")
    def nonnegative(raw):
        import math
        value = float(raw)
        if not math.isfinite(value) or value < 0:
            raise argparse.ArgumentTypeError("must be finite and nonnegative")
        return value
    cli.add_argument("--discard-ns", type=nonnegative, help="Exclude source frames at or before this elapsed production time in each replica (default: 0)")
    cli.add_argument("--stride", type=positive, help="Keep every Nth frame after discarding (default: 1)")
    cli.add_argument("--cpptraj", default=None, help="CPPTRAJ executable override")
    cli.add_argument("--no-plots", action="store_true", help="Tables/HTML without Matplotlib")
    cli.add_argument("--no-progress", action="store_true")
    args = cli.parse_args(argv)
    from .campaign import locked
    # Keep reset/submission from changing this campaign during direct analysis.
    with locked(Path(args.project).expanduser().resolve()):
        return _analyze_campaign(args)


def _analyze_campaign(args):
    from .project import campaign_config, portable_config, write_new_json
    from mdwb.config import load_config
    template = load_config(args.config, check_files=False) if args.config else None
    config = campaign_config(args.project, template=template, level=args.level,
                             discard_ns=args.discard_ns, stride=args.stride, output=args.output,
                             replica=args.replica, check_files=not (args.plan or args.save_config))
    if args.cpptraj:
        config["cpptraj"] = args.cpptraj
    if args.no_progress:
        config["progress"]["enabled"] = False
    if args.no_plots:
        config["plots"]["enabled"] = False
    if args.save_config:
        path = write_new_json(args.save_config, portable_config(config, args.save_config))
        print(f"Analysis configuration saved: {path}\nRun when production is complete: md-bud run {path}")
        return 0
    from mdwb.runner import make_plan, run_workflow
    if args.plan:
        print(f"Offline plan: {make_plan(config)}")
        return 0
    manifest = run_workflow(config, make_plots=not args.no_plots)
    print(f"Analysis {manifest['status']}: {config['output']}")
    return 0 if manifest["status"] == "complete" else 2


def doctor(argv):
    from mdwb.community import environment_info
    cli = argparse.ArgumentParser(prog="md-bud doctor", description="Check the dependencies for the work you want to do; no simulation is executed.")
    cli.add_argument("--task", choices=["generation", "analysis", "cluster", "all"], default="generation")
    cli.add_argument("--cpptraj", default="cpptraj")
    cli.add_argument("--engine", choices=["cpu", "cuda"], default="cuda", help="Production engine to check for cluster tasks (default: cuda)")
    cli.add_argument("--no-plots", action="store_true", help="Do not require optional Matplotlib")
    args = cli.parse_args(argv)
    info = environment_info(args.cpptraj)
    missing = []
    if sys.version_info < (3, 11):
        missing.append("Python >=3.11")
    info.update(package="MD-bud", version=__version__, task=args.task)
    if args.task in {"analysis", "all"}:
        for name in (["numpy"] if args.no_plots else ["numpy", "matplotlib"]):
            if not info["meets_minimum_version"][name]:
                missing.append(name)
        if not info["cpptraj_executable"]:
            missing.append("CPPTRAJ")
    if args.task in {"cluster", "all"}:
        names = ["bash", "flock", "realpath", "sbatch", "squeue", "sacct", "scancel", "srun", "pmemd.MPI"]
        if args.engine == "cuda":
            names.append("pmemd.cuda")
        info["cluster_executables"] = {n: shutil.which(n) for n in names}
        if os.name != "posix":
            missing.append("Linux/WSL cluster environment")
        missing.extend(n for n in names if not info["cluster_executables"][n])
    info["missing"] = missing
    print(json.dumps(info, indent=2))
    return 2 if missing else 0


def qc(argv):
    from .campaign import positive
    from .qc import inspect_campaign
    from .project import write_new_json
    cli = argparse.ArgumentParser(prog="md-bud qc", description="Inspect AMBER production logs without numerical packages or CPPTRAJ.")
    cli.add_argument("project")
    cli.add_argument("-R", "--replica", type=positive)
    cli.add_argument("--output", help="Save JSON to a new file; existing files are never overwritten")
    args = cli.parse_args(argv)
    result = inspect_campaign(args.project, args.replica)
    if args.output:
        write_new_json(args.output, result)
    print(json.dumps(result, indent=2, allow_nan=False))
    return 1 if result["issues"] else 0


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    try:
        if not argv:
            parser().print_help()
            return 0
        command, *rest = argv
        if command == "generate":
            return generate(rest)
        if command in CAMPAIGN:
            from .campaign import main as campaign_main
            return campaign_main(argv)
        if command in ANALYSIS:
            from mdwb.cli import main as analysis_main
            return analysis_main(argv)
        if command in {"analyze", "configure", "doctor", "qc"}:
            return globals()[command](rest)
        parser().parse_args(argv)
        return 0
    except (ValueError, OSError, RuntimeError, ImportError, KeyError, TypeError) as exc:
        print(f"md-bud: {exc}", file=sys.stderr)
        return 2
    except (KeyboardInterrupt, EOFError):
        print("md-bud: stopped; existing files and logs retained.", file=sys.stderr)
        return 130
