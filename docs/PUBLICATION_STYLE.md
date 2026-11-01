# Publication styling and optional solvation analyses (1.5)

## Beginner and detailed setup

Launch the wizard and choose **basic** for automatic biomolecular regions and
RMSD/RMSF/radius of gyration, with an optional extras menu. Choose **guided** for
all selection, alignment, reference and specialist controls. Neither path changes
the definition of the existing Basic/Standard/Advanced presets. SASA and NAStruct
remain outside Standard; RDF and watershell are explicit add-ons only.

```text
python mdworkbench.py wizard --quick --root simulations --config study.json
```

The short path still asks you to verify topology, trajectory grouping, production
start and saved-frame spacing. The default fit uses detected representative atoms
and the first selected frame of replica 1; use guided mode for domain-specific
alignment or a different reference. Systems without recognized biomolecular
regions require guided setup with explicit masks.

Setup/run summaries now display aligned labels and values, separate stages and
timestamped progress. Redirected cluster logs contain no added ANSI coloring;
`NO_COLOR` disables terminal color explicitly. `report` prints a readable summary;
the complete machine-readable output remains in `reports/report_summary.json`.

## Rendered figure defaults

- Arial preferred, with Liberation Sans then DejaVu Sans as local fallbacks.
- 11 pt axis labels, 12 pt titles and 10 pt ticks/legends on a 180 mm canvas.
- Lagoon categorical colors, with Mineral and the established colorblind palette
  as alternatives; choose based on your journal and reproduction requirements.
- White backgrounds, stronger traces, light gridlines and no legend frames.
- Consistent replica colors for structural traces, their overlays/distributions,
  pooled PCA overlays and autocorrelation overlays.
- 600 DPI for raster layers embedded in SVGs; vector lines/text do not gain
  resolution from increasing DPI.
- Editable SVG text by default, or outlined glyphs for stable appearance on a
  machine without the rendering font.

The figure customization question exposes palette, font, size and text mode.
For saved configurations, adjust these fields:

```json
"plots": {
  "enabled": true,
  "font": "Arial",
  "font_size": 11,
  "width_mm": 180,
  "palette": "lagoon",
  "alpha": 0.85,
  "dpi": 600,
  "svg_text": "editable"
}
```

For a single-column figure, consider `width_mm: 90` and `font_size: 9`, then inspect
the actual exported figure at its intended size. Tight cropping can change the
exported bounding dimensions from the nominal canvas. Check dense labels, legends
and multi-panel figures against your journal's specifications; there is no single
universal “2026 publication standard.” The palette choices are design presets, not
a claim that every color combination is distinguishable under every vision condition.

`reports/figure_style.json` and the report summary record requested/resolved fonts
and style settings after rendering. Fonts are neither installed nor bundled. If
Arial is unavailable, the fallback is explicit; use an existing Arial-equipped
environment to render genuine Arial glyphs. With `svg_text: paths`, glyphs retain
their appearance but cease to be directly editable as text.

## MDAnalysis-informed additions

The [MDAnalysis analysis catalog](https://docs.mdanalysis.org/stable/documentation_pages/analysis_modules.html)
and [RDF module](https://docs.mdanalysis.org/stable/documentation_pages/analysis/rdf.html)
highlight spatial distributions as a useful complement to structural observables.
This update adds independent cpptraj command generation, not an MDAnalysis backend
or copied MDAnalysis implementation. No MDAnalysis installation is required.

**RDF:** specify two disjoint selections, bin spacing and maximum radius. The
workflow uses [cpptraj radial](https://amberhub.chpc.utah.edu/radial-rdf/) with
average-box-volume normalization and produces a per-replica profile, matching-axis
mean/SD tables, SVGs and browser profile previews. A valid periodic box is required.
Verify that the maximum radius remains below half the shortest periodic cell height
throughout the selected trajectory; this geometric bound is not automatically
measured by the wrapper. A heterogeneous solute environment need not approach
g(r)=1. No automatic coordination number or RDF-derived free energy is claimed.

**Solvent shells:** [cpptraj watershell](https://amberhub.chpc.utah.edu/watershell/)
produces inner and outer cumulative counts at user-defined cutoffs. The outer
count includes the inner count; subtracting them gives the intervening shell
population. These are geometric occupancies, not residence times. Select complete
one-residue solvent molecules, such as `:WAT`, for consistent CPU/CUDA behavior:
the [native implementation](https://github.com/Amber-MD/cpptraj/blob/master/src/Action_Watershell.cpp)
uses residue counts on CPU and molecule counts on CUDA, which also requires whole,
equal-size solvent molecules. Other solvent models need deliberate review.

Both analyses run with the imaged coordinate/box orientation retained: the RMS
check uses `nomod`, so it does not rotate coordinates before periodic distances.
Native preflight checks nonempty and disjoint selections. Neither analysis is
silently enabled for a beginner.

## Further extensions

MSD/diffusion, residence times, transport and reweighted free energies remain
specialist work: unwrapped coordinates, continuity, estimator definitions or bias
weights need additional design and validation. They have not been added as simple
buttons that could produce misleading results from the aligned workflow.

Version 1.5 was reviewed statically. Native solvation calculations and updated
Matplotlib rendering have not been executed in this environment, and no tests or
dependencies were installed/run for this update. Existing examples and historical
validation records do not validate the new features.
