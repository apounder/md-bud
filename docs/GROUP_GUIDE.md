# Using the workbench in a research group

Version 1.4 embeds numerical previews and statistics in `index.html`; copying that
file alone is sufficient for offline viewing of its embedded content, including
when publication SVG rendering was disabled on the cluster. See the
[self-contained report guide](HTML_REPORT.md) for the new navigation and limits.

## First use

Use an existing environment with Python >=3.10, NumPy >=1.24, Matplotlib >=3.7
for plots, and CPPTRAJ for trajectory calculations. Nothing is installed automatically.
Run commands from the package directory, or replace `mdworkbench.py` with the path
to the supplied `mdworkbench.pyz` file.

```text
python mdworkbench.py doctor
python mdworkbench.py catalog
python mdworkbench.py wizard --root /path/to/simulations --config study.json
python mdworkbench.py validate study.json
python mdworkbench.py plan study.json --output review_plan
python mdworkbench.py run study.json
```

`doctor` reports package versions and executable availability, without launching
CPPTRAJ or importing numerical libraries. Presence does not guarantee binary or
package compatibility. `validate` reads the topology/configuration and lists
resolved selections and applicable analysis batches; it does not read trajectory
frames or execute native masks. The run performs native preflight checks.

The wizard asks for a study title, replica grouping, sections, a preset and numbered
additional analyses. Distance, angle and dihedral monitors have geometry-specific
questions. SASA and NAStruct remain outside Standard; Advanced or explicit add-ons
enable them. Confirm that sequential segments of one simulation are grouped as
one replica and that independent restarts belong to separate replicas.

## Read and share results

Open `analysis_results/index.html` in a browser after a run. No server, account or
internet connection is required. Filter filenames/figures, inspect calculation and
reporting status, and browse the data and reproducibility files. Calculations retain
the existing terminal progress display; the offline viewer is a completed-results
browser, not a live progress service or a FINE integration.

The most useful files are:

- `reports/all_replicas.dat`: tidy numerical export with replica labels.
- `reports/replicate_statistics.dat`: equal-weight replica summaries.
- `reports/diagnostics/sampling.dat`: per-observable sampling estimates and reasons
  an estimate was omitted. ESS is not the number of independent simulation replicas.
- `reports/diagnostics/replica_distribution_distance.dat`: distribution comparisons
  and links to the exact common-bin probabilities used.
- `figures/diagnostics/`: overlaid ACF curves and replica-distance heatmaps.
- `METHODS.md`, `environment.json`, `resolved_config.json`, `manifest.json`, native
  inputs and logs: a record to review before writing a methods section.

```text
python mdworkbench.py bundle analysis_results --destination study_results.zip
```

Use a new archive filename outside the results directory. The ZIP includes generated
results, prepared reference, scripts, logs and provenance, with portable artifact
paths. Original input trajectories/topology are not copied. Do not place unrelated
files in the results directory; generated-output folders are bundled recursively.
Provenance may expose local input paths and study names, so review before sharing.

Recipients extract the ZIP and open `index.html`. In an existing suitable Python
environment they can regenerate figures or data without rerunning CPPTRAJ:

```text
python mdworkbench.py report extracted_results
python mdworkbench.py report extracted_results --no-plots
```

The original environment record is preserved; `report_environment.json` records
the environment used for regeneration. An old calculation failure remains visible
even if report regeneration succeeds. Bundling is not an archive of the original
simulation and does not enable rerunning its calculations without original inputs.

## Optional settings

These fields may be added to the saved JSON configuration:

```json
{
  "study": {"title": "Protein–DNA replicate comparison", "description": "Production window selected after equilibration review."},
  "diagnostics": {
    "enabled": true,
    "max_series": 100,
    "max_samples": 1000000,
    "max_lag": 10000,
    "histogram_bins": 32,
    "max_pairs": 10000
  }
}
```

This is a settings fragment, not a runnable input: retain topology and replicas.
`max_samples` bounds one ACF series or all samples in one comparison group.
`max_series` bounds ACF attempts and distribution groups separately; `max_pairs`
bounds the total number of pairwise comparisons. Omissions are recorded explicitly.
These controls do not reduce the memory needed by the existing base table parser.
The default limits may omit less frequently inspected residue-level series; review
the omission tables and raise limits deliberately for the observables of interest.

## Cluster jobs and troubleshooting

Generate and review the JSON interactively once, then submit a noninteractive
`run study.json --no-progress` command using your group's scheduler template.
Use a different output directory for each job. There is no automatic scheduler
submission, concurrent replica executor or resume command in this version.

For failures, inspect `manifest.json` and the corresponding CPPTRAJ log. Exit code 2
indicates failure/partial run for `run`, invalid configuration, or missing prerequisites
for `doctor`; report warnings also appear in the report summary/viewer. `report`
returns zero when regeneration finishes, including outputs containing warnings.
An interrupted invocation returns 130. Do not interpret an absent artifact as a
zero measurement. Keep the saved input scripts when reporting an issue.
