# Merge and migration

MD-bud 2.0.0 was initially assembled in a separate workspace, then transferred
into the MD-file-generation checkout for review and publication. Its existing
Git history is retained. The MD-analyser source checkout remains unchanged.

| Source | Imported revision | Destination |
|---|---|---|
| [MD-file-generation](https://github.com/apounder/MD-file-generation) | `40b32d554a3649c734cc7eab69a81cb38aa85096` | `mdbud/campaign.py`, original regression tests and campaign reference |
| [MD-analyser](https://github.com/apounder/MD-analyser) | `4276792ac5a41c9f03347959e437ff5ebe814ecb` | `mdwb/`, analysis tests, documentation, examples and bundled viewer notices |

Both MIT copyright notices are preserved in the root license. Third-party
viewer assets and example-structure attribution remain in THIRD_PARTY.md.
No hosted GitHub repository has been created, overwritten, pushed or archived.

## Existing users

- Install this project once; the installed distribution is `md-bud`.
- Existing generation commands retain their flags. Existing campaign manifests
  retain schema version `1.0.0`, independently of the package's `2.0.0` version.
- Existing campaigns keep their original inputs and protocol. Merely opening or
  analyzing a campaign never regenerates its AMBER input files.
- Replace `python mdworkbench.py ...` with `md-bud ...`, or continue using the
  compatibility script/entry point. Analysis config schema remains version 1.
- `md-bud analyze CAMPAIGN` is the new direct handoff. Generic trajectory projects
  continue to use `wizard` and `run`. Existing complete analysis configurations
  can serve as `analyze --config FILE` templates.
- Source-only generation now uses `python3 md-bud` plus the two package folders.
  To distribute one file, build the unified `md-bud.pyz`; the old standalone
  generator and old analysis-only pyz are no longer needed.
- `md-bud` without arguments prints unified help. Use `md-bud wizard --quick`
  for the analysis wizard. Campaign analysis starts with the basic preset;
  the legacy wizard retains its own guided choices.

## Deliberate protocol changes for newly generated campaigns

1. NPT production explicitly uses the Monte Carlo barostat (`barostat=2`). Use
   `--barostat berendsen` to request the inherited behavior. Restrained equilibration
   explicitly uses Berendsen; the final unrestrained step uses the production choice.
2. A final unrestrained equilibration step (`eq-06`) follows the five inherited
   restrained steps. Its default duration is 0.1 ns, set with `--equilibration-ns`.
   This is a fixed protocol duration, not an automatic equilibrium test.
3. Production steps must divide exactly by `ntwx`. Otherwise segment boundaries
   introduce a longer gap between saved frames and a single frame interval would
   misrepresent the time axis. Old irregular campaigns can still be managed, but
   direct campaign analysis refuses them until explicitly reprocessed.

Other inherited settings include 310 K, a 2 fs production timestep, 10 Å cutoff,
the heating schedule, restraint schedule, Slurm resource defaults and Fir modules.
Review these against the intended physical model. Three replicas in examples
are a workflow illustration, not a universal sampling sufficiency claim.

## Architecture

`mdbud.cli` routes public commands. `mdbud.campaign` owns generation and scheduler
operations. `mdbud.project` maps manifest records into the existing `mdwb` analysis
configuration and checks input hashes, stage completion and segment order.
`mdbud.qc` streams production logs without extra Python dependencies. `mdwb`
continues to own selections, CPPTRAJ orchestration, scientific calculations and
reporting. Keeping these modules avoids duplicating trajectory or selection logic.

Campaign and analysis manifests serve different purposes and are retained:
the former records jobs and simulation inputs, the latter records analysis
choices and outcomes. Combining them into one mutable file would make scheduler
recovery and analysis reproducibility unnecessarily interdependent.
