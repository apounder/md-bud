# Analysis reference (imported from CPPTRAJ Workbench 1.7.0)

> This is the imported analysis feature reference. MD-bud now provides these
> commands in the same package as generation: use `md-bud wizard`, `md-bud run`,
> etc. The [root README](../README.md) and [migration notes](MIGRATION.md) describe
> current installation and campaign integration. Build the unified `md-bud.pyz`
> for one-file cluster deployment. Historical validation statements below refer
> to the source release; see [current validation](VALIDATION.md).

A guided, reproducible workflow for protein, DNA/RNA and mixed-system MD analysis
with CPPTRAJ. Version **1.7.0**.

## Quick start

Use an existing environment with **Python 3.10+**, **NumPy 1.24+** and **CPPTRAJ**
(AmberTools). **Matplotlib 3.7+** enables publication SVGs. Nothing is installed
automatically; run the source directly:

```text
python mdworkbench.py doctor
python mdworkbench.py wizard --root /path/to/trajectories --quick
python mdworkbench.py run analysis_config.json
```

The wizard discovers topology/trajectory files, asks how to group replicas and
segments, explains selections, and saves a reusable configuration. At “Run CPPTRAJ
now?”, press Enter to run (the default), or answer `n` to save it for later; the
wizard prints the command to run that file. Starting the wizard again detects an
existing configuration and offers to run it without repeating setup. It checks
`analysis_config.json` first, then other analysis JSON files in the current folder
and `--root` folder, offering a choice if several are found. Use `wizard --config
study.json` to run a specific saved configuration or create one with a new filename.
Omit `--quick` to choose between basic and guided setup.

For clusters, copy only `mdworkbench.pyz` and substitute that filename for
`mdworkbench.py` in these commands. The same environment prerequisites apply.

After analysis selection, both setup paths offer resource limits for PCA,
clustering and pairwise RMSD: maximum pooled frames, selected atoms and estimated
memory (GiB). The pooled total includes all replicas after trimming and stride.
If it exceeds the current cap, the wizard suggests a cap covering all selected
frames without changing the global stride. For three 20,000-frame replicas,
keeping every frame requires `advanced.max_frames` of at least 60,000.
Clustering has its own `advanced.cluster_sieve`: clusters are fitted to a random
1/N subset and remaining frames are assigned afterward. This leaves sampling
for other analyses unchanged. Keep the global analysis stride at 1 to analyze
every saved frame. Raising the frame cap does not bypass the memory guard.
DCCM offers atom and memory limits.
Existing saved JSON configurations can set these under `advanced` using
`max_frames`, `max_matrix_atoms` and `max_memory_gb`; running a saved configuration
does not repeat setup prompts.

## What it provides

- Any number of ordered trajectory segments per replica; combined `.dat` tables
  preserve replica identities and individual results.
- Basic presets plus optional hydrogen bonds, contacts, secondary structure,
  torsions, SASA, NAStruct, RDF, solvent shells, PCA, clustering and matrix analyses.
- Guided distance, angle and dihedral monitors; protein/backbone/nucleic selections,
  residue ranges and named regions.
- Individual, overlaid and concatenated SVGs where appropriate; coordinated
  palettes, configurable sizing and Arial preference with a recorded fallback.
- Probability distributions for scalar time series, with per-replica densities,
  matching histogram bins and numerical tables.
- Offline HTML with an analysis directory, trace/distribution/summary views, per-replica
  statistics and sampling guidance. Optional click-to-load **3Dmol.js** shows one
  reference structure, with atom/file limits and no trajectory animation.

Open `analysis_results/index.html` locally after a run. On a cluster, copy this
file to your computer: embedded previews work without a server or internet.
Keep the full output folder for linked raw data and supporting files.

```text
python mdworkbench.py report analysis_results --no-plots
python mdworkbench.py bundle analysis_results --destination study_results.zip
```

`--no-plots` retains numerical/HTML previews without Matplotlib. Open the included
[example report](examples/rendered/index.html) after downloading/cloning this
repository; it uses synthetic curves and an unrelated experimental structure.

See [convergence and geometry clustering](docs/CONVERGENCE.md) for cumulative/windowed
replica comparisons, pooled/last-ns/custom reference targets, and QM/MM geometry features.

[Pi-stacking from residue numbers](docs/STACKING.md) discovers candidate rings for
confirmation and naming, then reports per-pair/per-residue occupancy and replica mean/SD.

## Scope and guidance

One comparison requires a matching topology and the same atom identities/order
across replicas. You establish topology compatibility, replica independence and
scientific comparability; automatic discovery cannot prove these. Stable traces
do not prove convergence. Native CPPTRAJ and browser/WebGL behavior of this update
remain unverified; see the [validation record](docs/VALIDATION.md).

[HTML and 3D viewer](docs/HTML_REPORT.md) · [Group workflows](docs/GROUP_GUIDE.md) ·
[Figure styling](docs/PUBLICATION_STYLE.md) · [Configurations](examples/README.md) ·
[Citation](CITATION.md) · [Contributing](CONTRIBUTING.md) ·
[MIT license](LICENSE) · [Third-party attribution](THIRD_PARTY.md)
