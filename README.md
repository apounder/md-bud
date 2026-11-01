<p align="center">
  <img src="./banner.png" alt="MD-bud banner" width="1200">
</p>

# MD-bud

One package for **AMBER simulation campaigns and CPPTRAJ analysis**, designed
primarily for **Linux Slurm clusters**. Generate inputs, manage replica jobs,
check production logs, and analyze results with one command-line tool.

MD-bud combines MD-file-generation and MD-analyser. Their campaign and analysis
formats remain supported. There is one installation and one public CLI; the
existing analysis engine is kept as an internal module rather than rewritten.

## Install on the cluster

Use an existing **Python 3.11+** environment. From this checkout:

```bash
# Generate/manage campaigns: standard library only
python3 -m pip install .

# Recommended: analysis, numerical tables, SVG figures and offline HTML reports
python3 -m pip install '.[plots]'

md-bud --help
md-bud doctor --task generation
md-bud doctor --task analysis
```

`'.[analysis]'` installs NumPy without Matplotlib; use `--no-plots` for analysis.
AMBER, AmberTools/CPPTRAJ, Bash, `flock` and Slurm are **external cluster software**;
pip does not install them. Load the site's modules before running relevant commands.
Generation and offline planning do not require AMBER, Slurm, NumPy or Matplotlib.

You can also run `python3 md-bud ...` or `python3 -m mdbud ...` from the checkout.
For a cluster without installation permissions, see [single-file deployment](docs/CLUSTER.md).

## A complete campaign

Start with a **prepared periodic explicit-solvent AMBER topology and coordinates**.
Structure preparation, protonation and force-field parameterization come first.

```bash
# Save your own allocation once; add module/resource overrides if needed
md-bud configure cluster.json --account your-allocation

# Three replicas, five sequential 10 ns segments each: 50 ns per replica
md-bud generate --profile cluster.json -i runs/study \
  -p system.prmtop -s system.inpcrd -R 3 -n 5 -t 10

# Review the generated .in/.sh files, then preview and submit
md-bud submit runs/study --dry-run
md-bud submit runs/study
md-bud status runs/study

# After production, in a suitable compute allocation
md-bud frames runs/study
md-bud qc runs/study
md-bud analyze runs/study
```

Open `runs/study/analysis_results/index.html` after copying the results to your
computer. The report works offline. Keep the whole output folder for linked data.
Analysis can consume substantial CPU/RAM; [submit it to a compute node](docs/CLUSTER.md).

Each replica follows `min → heat → eq → md-0001 → md-0002 → …`.
Dependencies require predecessor success. Inputs are copied and checksummed.
Generation does not submit jobs unless explicitly requested. Existing campaigns
and analysis outputs are protected against replacement; reset archives old outputs.

There is **no need to combine trajectories first**. Campaign analysis reads only
the manifest's ordered production segments, retaining each replica's identity.
Heating/equilibration trajectories and optional `combined.nc` are excluded.
The time interval comes from `dt × ntwx`, rather than from filenames or the
integration timestep alone. Campaign analysis requires all selected production
stages to be complete, verifies copied-input checksums and checks native frame
counts against the expected count for every segment before analysis.

## Easy tasks and sensible defaults

```bash
# Review analysis before trajectories exist; use a separate output directory
md-bud analyze runs/study --plan --output runs/study/analysis_preview

# Export a reusable, editable configuration without running
md-bud analyze runs/study --save-config runs/study/analysis_config.json

# Remove the first 5 ns of production from every replica, keep every fifth frame
md-bud analyze runs/study --discard-ns 5 --stride 5

# More structural analyses, or numerical results without Matplotlib
md-bud analyze runs/study --level standard
md-bud analyze runs/study --no-plots

# Trajectories from an existing project or a different engine
md-bud wizard --root /path/to/trajectories --quick
md-bud run analysis_config.json

# Review a reset, then archive the affected outputs and cancel affected jobs
md-bud reset runs/study -R 2 --from md-3
md-bud reset runs/study -R 2 --from md-3 --apply --cancel
```

Repeated analysis examples above are alternatives: choose a fresh `--output`
directory for each subsequent run. A preview directory also cannot be reused for
execution. `--save-config` creates a new file, with paths relative to that file,
so moving an entire campaign retains its internal links.

| Default | Why / what to review |
|---|---|
| Analysis: RMSD, RMSF, radius of gyration | A small structural starting point; larger presets/add-ons are available. |
| All production frames, stride 1 | No automatic discard or downsampling. Choose a discard after assessing relaxation. |
| Automatic protein/nucleic sections and representative fit | Modified residues, ligands and unusual complexes may need explicit masks. |
| Periodic imaging and shared reference | Review the anchor for multichain complexes; all replicas share the same reference. |
| 310 K, NPT, 1 bar, 2 fs, 10 Å cutoff | Retained biochemical starting settings, not a protocol for every system. |
| Monte Carlo production barostat | Explicit `barostat=2`; restrained density relaxation uses Berendsen. |
| Final 100 ps unrestrained equilibration | Added after the original restrained stages; duration is configurable and does not establish equilibration. |
| Saved coordinates every 50 ps | `ntwx=25000` at 2 fs. Use shorter intervals for fast motions; this is too sparse for many kinetic analyses. |
| One replica, one 10 ns segment | Kept for compatibility. Choose multiple independently equilibrated runs and study-specific lengths. |
| Saved, distinct stochastic seeds | Reproducible input choices; seeds alone do not prove replica independence. |

Cluster defaults retain the original Fir module stack. On another site, pass
`--modules "YOUR MODULE LIST"`, or `--modules ""` to inherit a loaded environment.
Profiles accept generation option names, e.g. `account`, `modules`, `mem`,
`temperature`, `replicas`, `ntwx`; explicit CLI options override profile values.
Use `--profile FILE` or set `MDBUD_PROFILE`. Personal settings belong outside Git.

## Available functions

| Area | Commands / capabilities |
|---|---|
| Campaigns | `generate`, `submit`, `status`, `reset`, `frames`, `configure` |
| Checks | `doctor --task generation\|analysis\|cluster\|all`, `qc` |
| Analysis | `analyze`, `wizard`, `discover`, `validate`, `catalog`, `plan`, `run` |
| Reports | `report`, `bundle`, `demo`; numerical tables, SVGs, offline HTML, bounded reference-structure preview |
| Structural add-ons | H bonds, contacts, secondary structure, SASA, nucleic geometry, torsions, distances/angles/dihedrals, pi-stacking, RDF and solvent shells |
| Sampling/structure comparisons | Replica summaries, distributions, descriptive convergence curves, correlation diagnostics, shared PCA and clustering, DCCM and pairwise RMSD with resource guards |

For custom analyses, export an analysis JSON and use `md-bud run FILE`, or pass
an existing complete configuration as `analyze --config FILE`. Campaign topology,
trajectories, output location and saved-frame interval take precedence;
`--level` replaces the template's analysis choices. `catalog` explains the
available calculations. Run-specific configurations, logs and methods records
retain the analysis choices.

## Scientific scope and remaining gaps

This is a conventional biomolecular workflow, not universal MD automation.
Stable RMSD, a readable trajectory and a completed scheduler job do not establish
physical validity or convergence. Independent replicas, correct atom identities,
unit-cell geometry and suitable system preparation remain essential.

`qc` summarizes production-log temperature, energy, volume, pressure and density
where present, flags numerical/engine errors, and reports missing logs. Its
20 K mean-temperature flag is a coarse heuristic. It does not determine
equilibration or test NVE energy drift. No density/pressure threshold is imposed.

The [computational chemistry review](docs/SCIENTIFIC_REVIEW.md) separates covered
workflows from remaining needs: preparation validation, adaptive equilibration,
membranes, enhanced sampling/free energies, transport/kinetics, and system-specific
chemical observables. [Analysis science notes](docs/SCIENCE_NOTES.md) document
selection, imaging, reference and statistical assumptions.

## Compatibility, testing and provenance

Existing `md-bud.json` campaign schema 1.0.0 and analysis JSON schema 1 are retained.
Old campaigns keep their original generated simulation inputs. New campaigns
record MD-bud 2.0.0 as producer. `mdworkbench` and `mdworkbench.py` remain analysis
compatibility entry points. See [migration](docs/MIGRATION.md) for protocol changes.

```bash
python3 -m pip install '.[dev]'
python3 -m unittest discover -s tests -v
python3 -m build
python3 tools/build_zipapp.py
```

See [validation](docs/VALIDATION.md) for observed results and native-engine limits.
Linux and Windows CI workflows are included; their presence does not imply they
have run on GitHub. AMBER/CPPTRAJ chemistry and real Slurm execution require a
site-specific smoke test. Source origins and preserved licenses are recorded in
[migration](docs/MIGRATION.md), [LICENSE](LICENSE) and [third-party notices](THIRD_PARTY.md).

Detailed reference: [campaign options](docs/CAMPAIGNS.md) ·
[analysis reference](docs/ANALYSIS.md) · [cluster deployment](docs/CLUSTER.md) ·
[convergence](docs/CONVERGENCE.md) · [stacking](docs/STACKING.md) · [citation](CITATION.md).
