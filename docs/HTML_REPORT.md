# Self-contained HTML reporting (1.4)

Generate the report on the cluster, then copy `index.html` to your computer and
open it in a browser with JavaScript enabled. All numerical previews and statistics
are embedded in this one file: no server, file fetching, internet, display server or
Matplotlib is required to view them. A compute node itself still needs a browser
and display if you want to view the report there; normally you view it locally.

```text
python mdworkbench.py run study.json --no-plots --no-progress
python mdworkbench.py report analysis_results --no-plots
```

The first command performs a new calculation; the second regenerates reporting
from existing numerical outputs. `--no-plots` skips Matplotlib publication SVGs
while retaining browser-native previews. NumPy is still required for reporting.
Neither command installs anything.

## Navigation

- **Study overview:** replica coverage, status and review-item count.
- **Explore results:** choose analysis, section/observable and one replica or an
  overlay of all replicas; the charts and exact statistics update together.
- **Replica statistics:** filter by analysis, replica and section/observable text;
  individual, between-replica, circular and categorical summaries stay separate.
- **Sampling review:** recorded evidence, conditional estimates and suggested
  inspection, with no automatic convergence verdict.
- **Full figures & files:** embedded publication SVGs when available, plus links
  to original outputs for users who have the complete results bundle.

The figure and data links require the bundle; copying the HTML alone preserves
the embedded content. Available publication SVGs are embedded up to a cumulative
20 MB source-size budget. Any additional SVGs remain linked in the bundle.

## Preview semantics and limits

Traces retain separate replica coordinates and use ordered min/max envelopes to
reduce file size, preserving extrema within fully finite buckets. A bucket with
missing coordinates or values becomes a gap rather than a false connecting line.
Circular traces break at the angle boundary. These reduced traces are inspection
previews, not replacements for the full numerical series or publication figures.

Histograms use all finite observations and a common 32-bin grid across included
replicas. Circular distributions use [-180,180); categorical data use their actual
state identifiers. State identifiers are displayed explicitly without inventing
biochemical names. Histogram heights are fractions, not probability densities.
For profiles, the distribution is over profile values rather than trajectory time.

Mean bars show full-data linear observable means with SDs in the nearby table.
Circular observables retain circular summaries. Pie charts are enabled only for
categorical fractions of a single selected replica, so independent trajectories
are not silently pooled. Hover on bars or pie sectors to inspect values.

Default limits are 500 replica/observable series and at most approximately 600
trace points per series. A complete observable group is omitted if including all
its replicas would exceed the limit, and the overview states the omitted count.
Categorical groups exceeding 100 distinct states are omitted from previews.
Matrix previews are not generated in the browser; matrix data and available
publication figures remain accessible. Exact statistics are independent of the
trace reduction and preview-series limit.

Change the preview limits in the configuration before generating reports:

```json
"reports": {"preview_series": 1000, "preview_points": 1200}
```

Larger limits increase HTML size and browser memory use. Statistical tables display
up to 500 matching rows at once; narrow the filters or use the full data files.
These UI limits do not change calculations. The base numerical parser still holds
individual tables in memory.

## Geometry setup

The numbered extra-analysis menu offers distance, angle and dihedral monitors.
Each supports guided topology residue numbers and atom names, numbered named
regions, or quick semicolon-separated selections:

```text
Distance: 25:CA; 41:CA
Angle:    25:CA; 41:CA; 60:CA
Dihedral: 25:N; 25:CA; 25:C; 26:N
```

These are illustrative selections, not assertions that those atoms exist in your
system. Guided mode lists the available atom names when topology metadata is
readable; indices refer to the topology, not necessarily original PDB numbering.
Angles use point 2 as the vertex; torsions use points 2–3 as the central axis.
Review the displayed ordered selections before saving. Group-center choices and
periodic imaging retain their explicit questions. Each setup question now has a
one-sentence explanation before the prompt.

## Scientific interpretation

Within-replica SD measures trajectory fluctuations; SD of replica means measures
the spread of available independent simulation averages with equal replica
weight. Neither is automatically a confidence interval. A single replica has no
between-replica SD, and two replicas give a fragile estimate of that spread.

Review rules flag limited replica coverage, missing numerical outputs, omitted
comparisons, unavailable sampling estimates, estimated ESS below 50, lag-limit
truncation and blocks below the ceil(5*g) heuristic. Thresholds are inspection
rules rather than validated convergence cutoffs. Half-trajectory mean differences
and distribution distances are descriptive; no significance threshold is applied.

The report concept draws on the combination of results, diagnostics and
reproducibility in [ROBERT](https://github.com/jvalegre/robert) and its
[2024 paper](https://doi.org/10.1002/wcms.1733). This is independent CPPTRAJ reporting
code and does not use ROBERT's software or its machine-learning scoring system.

The 1.4 synthetic design example was generated without publication plotting.
Browser appearance review was blocked by the browser security policy; no scientific
validation, native CPPTRAJ run or test-suite execution is claimed.


## Analysis menus and reference structure (1.6)

Explore results groups analyses into structure/flexibility, geometry/torsions,
interactions/solvation, secondary/nucleic structure, collective motion/matrices,
and other outputs. Select a family, then an analysis, observable and replica.
Replica statistics and full figures also use grouped analysis dropdowns;
the gallery collapses figures by analysis. Matrix records remain accessible via
the supporting tables and full figures rather than scalar trace previews.

New runs attempt to export one common reference frame as `reference/viewer.pdb`.
Only automatically detected heavy solute atoms are retained; solvent and ions
are omitted. Stripping may renumber residues: use original topology indices for
monitor setup. This static frame provides orientation, not evidence of sampled
states or a representation of every replica. Older reports without this export
show an unavailable message; regenerating HTML alone cannot create coordinates.

Click **Load 3D structure** on the overview. Drag to rotate and scroll to zoom;
choose cartoon, bond lines or sticks, and reset the view as needed. Use bond
lines for structures without a recognizable cartoon backbone. Rendering happens
on your local browser and requires JavaScript and WebGL, not a cluster display.

The pinned 3Dmol.js asset adds about 700 KB after base64 embedding, plus PDB text.
Default limits are 12,000 atoms and 1,500,000 input PDB bytes. Oversized structures
are omitted completely, and no trajectory movie or surface is generated. The
library is omitted when no suitable structure is available. Optional export or
WebGL failures leave numerical analysis usable. To disable or adjust embedding,
add this top-level object to your configuration:

```json
"structure_view": {"enabled": false, "max_atoms": 12000, "max_bytes": 1500000}
```

Re-enable with `true`. Increasing limits increases report size and browser cost.
No external asset request is needed when opening the report; third-party
attribution is embedded and recorded in `THIRD_PARTY.md`.


## Organized result views and probability distributions (1.7)

Use the overview's **Analysis directory** to select an analysis. In Explore results,
choose **Traces & profiles**, **Probability distributions**, or **Replica summaries**.
The region/observable and replica filters apply to all three views. The gallery
filters SVGs by analysis, figure type and a region/replica search; readable titles
appear above each image and the original filename remains in its caption.

Radius of gyration, RMSD, distances, angles, SASA and other numeric time series
receive empirical distributions automatically. Use the vertical-axis selector
to switch between probability density (integral = 1) and probability per bin
(sum = 1). Replicas share 32 equal-width bins over their joint finite range;
constant observables use one enclosing bin. Each replica is normalized separately,
so longer trajectories do not dominate a pooled curve. Missing samples are excluded
and counted. These distributions do not establish stationarity or convergence.

Circular torsions use 36 common bins on −180° to 180°; the endpoints are adjacent.
Categorical state IDs retain probability bars/pies rather than continuous densities.
RMSF, RDF and other spatial profiles retain their coordinate interpretation and
are not shown as trajectory probability distributions. An RDF is already a distinct
spatial correlation function and should not be confused with a normalized histogram.

With Matplotlib, standalone density SVGs appear in `figures/distributions/`, per
replica and as overlays. Shading is decorative; it is not a confidence interval.
The renderer uses [Matplotlib stairs](https://matplotlib.org/stable/api/_as_gen/matplotlib.axes.Axes.stairs.html)
to draw the exact histogram bins, with no smoothing or additional dependency.
`reports/probability_distributions.dat` records edges, counts, sample sizes,
probabilities, densities and units even with `--no-plots`.

Regenerate existing reports from their stored numerical outputs:

```text
python mdworkbench.py report analysis_results
```

Add `--no-plots` to update only numerical/HTML outputs. SVG titles update only when
the figures are regenerated with Matplotlib. Native data labels and source files
remain unchanged for reproducibility.
