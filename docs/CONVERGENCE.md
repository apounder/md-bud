# Convergence over time and geometry clustering

The workbench prepares and runs CPPTRAJ analyses for ordered trajectory segments
within replicas, retains replica boundaries, and reports numerical tables and
HTML/SVG views. PCA and clustering use a shared pooled basis. This update adds
observable-specific time-dependent convergence diagnostics to that workflow.

Both basic and guided setup ask about convergence and clustering after the
analysis preset/add-ons. Existing configurations get pooled convergence defaults.
Use `convergence.enabled: false` to disable these extra reports independently of
`diagnostics.enabled`. Saved JSON is reusable with `run` and `report`.

## Early trajectory length check and trimming

Immediately after grouping input trajectories, the wizard uses CPPTRAJ `-tl` to
read each segment's frame count and compares total joined lengths per replica.
An explicit **TRAJECTORY LENGTH MISMATCH** appears before range/analysis prompts.
Unequal individual segment sizes are shown, but do not imply unequal replica
lengths if their sums match. This reads trajectory metadata; no coordinate analysis
or trajectory rewriting is performed.

Then set the saved-frame interval and choose a frame or ns range. You can discard
the beginning, stop earlier, or choose the shortest replica's endpoint for all.
The same range and stride apply across each replica's joined segments. Frame 1 is
0 ns; inclusive time bounds round inward to actual saved frames. The selected
frame counts are displayed before analysis options. An explicit endpoint must
fit every replica; `-1` preserves each replica's own endpoint.

If CPPTRAJ is unavailable, setup explicitly says the lengths are unverified and
still saves your requested range. `run` rechecks lengths before atom-selection
validation or analysis, flags differences, and records the warning in the manifest.
It applies your saved range without introducing interactive prompts into batch
jobs. Equal frame counts mean equal durations only under the required common
saved-frame interval; they do not establish topology or ensemble compatibility.

## Reference choices

The optional **coordinate reference path** sets the existing top-level `reference`
for structural RMSD. A PDB/restart is not a statistical mean/SD reference.

`convergence.reference_mode` supports:

- `pooled`: full-run equal-replica mixture. Mean = average of replica means;
  variance = average of each replica's mean squared deviation from that mean.
  SD includes within- and between-replica fluctuations. It is population SD,
  not SD of replica means and not uncertainty in the average.
- `tail`: the same calculation using each replica's last `tail_ns` ns, with
  interval `(end - tail_ns, end]`. Requires known saved-frame spacing and the
  full requested duration in every replica. No silent fallback to frame counts.
- `custom`: metric-specific mean, SD and units supplied in `references`.
- `file`: a JSON file at `reference_path` containing that same references map.
  Relative paths resolve from the configuration's directory.

Pooled and tail references use the data being assessed, so agreement with them
is descriptive and cannot serve as independent validation. In particular,
a single replica's final cumulative mean necessarily equals its pooled target.

Example configuration fragment (use exact IDs and units from
`reports/all_replicas.dat`; the series name below is illustrative):

```json
{
  "convergence": {
    "enabled": true,
    "checkpoints": 20,
    "window_frames": 50,
    "metrics": ["rmsd", "rg", "monitor_distance", "cluster"],
    "reference_mode": "custom",
    "references": {
      "rmsd/protein/RMSD": {"mean": 1.8, "sd": 0.3, "unit": "angstrom"}
    }
  }
}
```

For a file reference, the file contains just the object inside `references`.
An empty `metrics` list includes all eligible time series. An entry can be an
analysis name or an exact `analysis/section/series` ID. A missing or unit-mismatched
custom target leaves reference comparisons unavailable and records an omission;
within-replica statistics and between-replica comparisons still run. Mean and SD
alone do not define a distribution, so custom/file targets do not get JS distances.

## Outputs and interpretation

`reports/convergence/` contains:

- `within_replica.dat`: cumulative and trailing-window mean/SD, differences from
  the reference, reference distribution distance, and stationary ESS/SEM estimates.
  Every replica's curve extends to its own endpoint. Windows contain up to
  `window_frames` analyzed samples; `n` reports their actual size.
- `between_replicas.dat`: cumulative and windowed pairwise Jensen–Shannon (JS)
  distances, plus classical split R-hat across all replicas of each linear metric.
  Comparisons end at the shortest elapsed duration. R-hat uses equal-length
  prefix halves; unmatched extra samples are excluded from that statistic only.
- `state_populations.dat`: cumulative/windowed fractions using shared categories.
- `references.dat` and `reference_histograms.dat`: actual target statistics,
  units and common histogram bins/probabilities used in comparisons.
- `omissions.dat`: missing/irregular data, insufficient duration, reference
  mismatches, non-shared PCA/cluster bases and resource-limit exclusions.

Open **Sampling and convergence review** in the HTML report. SVGs in
`figures/convergence/` show mean differences and within-/between-replica
JS distances over time and are included in the report gallery. Tables remain
available with `--no-plots`.

Linear time series such as RMSD, radius of gyration, distances, SASA, contact or
hydrogen-bond counts, solvent-shell counts and shared PCA projections can be
assessed. Static RMSF profiles, matrices, RDF curves and summary tables are not
time series and are excluded. Dihedrals use circular moments and fixed periodic
bins; categorical states use fractions. Neither receives linear ESS or R-hat.
Custom mean/SD targets currently apply to linear metrics only.

A fixed common histogram grid is retained throughout a metric's curves.
JS distance uses base-2 logarithms and ranges from 0 to 1. It is sensitive to bins
and sampling and is not a significance test. Missing data are not removed to
manufacture contiguous sampling. Resource caps reuse the diagnostics settings.

Classical split R-hat compares within- and between-half variances; constant
halves are reported as unavailable. This is **not** modern rank-normalized/folded
R-hat and can miss differences in distribution shape or scale. Inspect JS curves
and raw distributions alongside it. ESS and SEM assume stationary sampling.
Independent, comparable replicas sampling the same ensemble are required for
scientific interpretation. No threshold automatically declares convergence.

## Cluster choices for structural and QM/MM studies

`advanced.cluster_metric` supports `rms` (best-fit RMSD), `dme` (distance-RMSD),
`srmsd` (symmetry-corrected RMSD), and `features` (selected geometry datasets).
`advanced.cluster_algorithm` supports `kmeans` or `hieragglo` (average linkage).
Cluster count and sieve remain configurable. `advanced.cluster_sieve` fits a
random 1/N subset and assigns the remaining frames afterward, without changing
PCA or other analysis sampling. Keep `frames.stride` at 1 for all saved frames;
the resource-limit prompt never changes it. Coordinate clustering uses
`advanced.matrix_mask`, which can select a ligand, active site or QM/MM region.

Feature clustering chooses existing monitors, or creates them in the wizard:

- bond lengths, donor–acceptor and metal coordination distances: `distance`
  monitors with two ordered atom masks;
- bond/coordination angles: three masks, the middle one defining the vertex;
- reaction torsions and improper-dihedral geometry: four ordered masks.

Groups may also define center distances and angles. Specify center convention,
imaging and atom ordering deliberately. This measures geometry; it does not
infer chemical bonds, coordination numbers or QM energies.

```json
{
  "analysis": {"additional": ["cluster", "distance", "dihedral"]},
  "monitors": [
    {"name": "metal_ligand", "type": "distance", "masks": [":1@ZN", ":25@NE2"]},
    {"name": "reaction_torsion", "type": "dihedral",
     "masks": [":50@C1", ":50@C2", ":50@C3", ":50@O1"]}
  ],
  "advanced": {
    "cluster_metric": "features",
    "cluster_algorithm": "kmeans",
    "cluster_features": [
      {"monitor": "metal_ligand", "weight": 1.0},
      {"monitor": "reaction_torsion", "weight": 0.01}
    ]
  }
}
```

Weights are passed explicitly to CPPTRAJ's combined Euclidean metric. They are
not automatic standardization: degree and angstrom scales need deliberate
weighting. CPPTRAJ geometry actions retain the dihedral dataset type for its
periodic torsion metric. Features are measured before stripping so original
atom masks remain valid. Selected coordinates are retained for representatives.
Shared pooled labels are then split back into replica time series, enabling
comparable population curves. Cluster populations depend on the selected
features, weights, cluster count and sieve; inspect sensitivity to those choices.

Methods and command references:
[CPPTRAJ clustering](https://amberhub.chpc.utah.edu/cluster/),
[CPPTRAJ dataset types](https://amberhub.chpc.utah.edu/dataset/),
[Stan convergence analysis](https://mc-stan.org/docs/reference-manual/analysis.html).
Native CPPTRAJ execution of the new options still needs validation with a real
trajectory; tests here verify command planning and synthetic numerical behavior.

## What does lag mean?

Lag is the separation between observations of the same metric in one replica.
Lag 0 compares each value with itself; lag 1 compares neighboring analyzed
frames. If saved frames are 10 ps apart and the analysis stride is 5, one lag
step is 50 ps (0.05 ns), and 20 steps correspond to 1 ns. The autocorrelation
axis uses physical time when the saved-frame interval is known, otherwise
analyzed-frame units.

Positive autocorrelation means values at that separation tend to vary together.
Slow decay indicates fewer effectively independent samples. Near-zero
correlation at a lag does not establish convergence or agreement between
replicas. Lag is not elapsed simulation time or an equilibration cutoff.
The maximum lag limits the separation examined, not the trajectory frames used.

SVG exports now include the metric, region and plot type, such as
`radius-of-gyration__active-site__<short-id>__autocorrelation.svg` or
`radius-of-gyration__active-site__<short-id>__convergence-over-time.svg`.
The short ID prevents collisions between similar labels. Regenerate an existing
report with `python mdworkbench.pyz report analysis_results` to create the new
filenames and links without rerunning CPPTRAJ. Older SVG exports may remain in
the folder; the regenerated report indexes the current exports.
