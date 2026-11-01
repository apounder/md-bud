# Contributing

This merged source bundle has not been published by this task. Add the actual
repository/maintainer details when publishing; do not invent contact addresses
or software DOIs. Include MD-bud and native CPPTRAJ versions in reports.

For a bug report, provide the smallest shareable configuration, expected/observed
behavior, `manifest.json`, relevant CPPTRAJ input/log and environment versions.
Remove confidential paths and structures. A short, legally shareable trajectory
with its matching topology is preferable to an entire production simulation.

Campaign/Slurm operations live in `mdbud/campaign.py`; the public CLI in
`mdbud/cli.py`; manifest-to-analysis integration in `mdbud/project.py`; lightweight
production-log inspection in `mdbud/qc.py`.
Analysis generators live in `mdwb/analysis.py`; orchestration in `runner.py`;
numeric interpretation and SVGs in `reporting.py`; sampling/distribution methods in
`diagnostics.py`; configuration and prompts in `config.py` and `cli.py`.
Keep CPPTRAJ command references beside new analysis definitions. Specify units,
periodic-boundary behavior, atom correspondence, angular wrapping and failure modes.
Do not apply linear statistics to state labels or circular torsions.

Before releasing changes, run the regression suite in an existing environment:

```text
python -m unittest discover -s tests -v
```

Linux is the primary target. The suite includes Bash/Slurm process doubles,
campaign-to-analysis integration and numerical/reporting regressions. See
`docs/VALIDATION.md` for observed runs. Add native small protein and nucleic-acid
acceptance tests with suitable redistribution rights before claiming engine-level
validation. Keep native engine results separate from mocked process boundaries.

Document observed validation separately from proposed checks. Do not add a passing
CI badge until a real repository workflow has passed. Use the MIT license for this
original workbench code; external software retains its own licensing and citations.
