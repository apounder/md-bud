# Changes

## 2.0.0 — MD-bud integration, 1 November 2026

- One installable package, public CLI and portable cluster zipapp for AMBER/Slurm
  generation and the full CPPTRAJ analysis/reporting engine.
- Direct campaign analysis with replica/segment mapping, copied-input checksums,
  expected native frame counts and correct saved-frame interval.
- Basic structural default, discard/stride, offline preview, movable exported
  configurations and support for complete custom analysis templates.
- Reusable cluster profiles, task-specific dependency checks and streaming
  production-log QC. Numerical/plotting dependencies are optional.
- Explicit MC NPT production pressure control and final configurable unrestrained
  equilibration for new campaigns; existing campaigns retain their inputs.
- Guard uniform segment output spacing and prevent resets during direct analysis.
- Preserve both licenses, origin revisions, original analysis API and legacy
  campaign/analysis formats. Fix inherited reporting/path regression failures.
- Linux/Windows test workflow, scientific scope review and cluster deployment guide.

## 1.7.0 — 2026-09-05

- Added an overview analysis directory, separate trace/distribution/summary views,
  and gallery filters by analysis, figure type and region/replica text.
- Added readable scientific SVG titles, axis labels and gallery titles while
  preserving native data identifiers and filenames.
- Added standalone probability-density SVGs per replica and as overlays, plus
  `reports/probability_distributions.dat` even without Matplotlib. HTML and SVG
  histograms share bin edges, use full finite samples, and normalize each replica
  separately. Categorical states remain probabilities; spatial profiles are
  excluded from temporal probability distributions.
- Rebuilt the synthetic example and checked Python/JavaScript syntax. No test
  suite, native CPPTRAJ calculations, Matplotlib rendering or browser checks.

## 1.6.0 — 2026-09-05

- Grouped report analyses into family dropdowns, including replica statistics
  and the full-figure gallery.
- Added an offline, click-to-load 3Dmol.js reference viewer with cartoon, line
  and stick styles. Export is optional and limited to one heavy-solute frame;
  default embedding limits are 12,000 atoms and 1.5 MB of PDB data.
- Bundled a pinned viewer asset with upstream licenses; no runtime CDN or new
  Python dependency. Structure-export failures do not stop numerical analyses.
- Added a rebuildable synthetic example with an explicitly unrelated experimental
  crambin structure, a concise README and repository ignore rules.
- Reviewed code and generated the example; no test suite, native CPPTRAJ run,
  Matplotlib rendering or browser/WebGL verification was performed.

## 1.5.0 — 2026-09-05

- Added a basic setup path (`wizard --quick`) with optional extras.
- Added opt-in RDF and solvent-shell analyses, explicit selections and cutoff
  validation, native mask preflight and box-orientation-preserving batches.
- Added formatted terminal summaries and timestamped run messages.
- Added Arial preference with recorded fallback, font/width controls, coordinated
  palettes, 600 DPI raster defaults and editable/outlined SVG text options.
- Static review only; no native solvation runs, rendering checks or tests executed.

## 1.4.0 — 2026-09-05

- Replaced the main HTML gallery with a self-contained, responsive results interface.
- Added embedded browser-native SVG previews from numerical data, including traces,
  common-bin histograms, replica mean bars and single-replica categorical pie charts.
- Added analysis/observable/replica selectors and filterable exact statistics.
- Embedded available publication figures with a 20 MB source-size budget.
- Added bounded trace previews and explicit omissions; statistics/histograms use
  full finite data of included observables.
- Added guided residue/atom picking, named-group picking, quick geometry syntax,
  an ordered-geometry review, and one-sentence help before setup questions.
- Generated a synthetic no-Matplotlib design example; browser appearance review
  was blocked by the browser's local-file URL policy, and no native runs or test
  suite were executed.

## 1.3.0 — 2026-09-05

- Added dedicated analysis, sampling-review and replica-comparison HTML pages.
- Added numerical overview tables, interpretation guidance and evidence-linked
  review findings, with explicit heuristic thresholds and no convergence score.
- Associated generated figures with their analyses directly during reporting.
- Static code review only; no tests, native runs or rendering checks executed.

## 1.2.0 — 2026-09-05

- Added `doctor`, `validate`, `catalog`, `bundle` and `report --no-plots`.
- Added an offline searchable HTML viewer for generated data, figures and logs.
- Added study titles, methods records and run/report environment provenance.
- Added per-replica ACFs, conditional sampling estimates and resource-limit records.
- Added shared histogram probabilities, replica Jensen–Shannon distances, ACF
  overlays and distance heatmaps.
- Added portable artifact paths in result archives and reporting after relocation.
- Added a literature review, research-group guide, contribution/citation guidance
  and MIT license.
- Static code review only for this version; no tests or calculations executed.

## 1.1.0

- Removed SASA and NAStruct from Standard; retained Advanced and explicit add-ons.
- Added a numbered extra-analysis menu and distance/angle/dihedral monitor prompts.
- Added terminal progress without additional dependencies.

## 1.0.0

- Initial replica-aware CPPTRAJ workflow, numerical reporting and SVG generation.
- See `docs/VALIDATION.md` for historical validation details.
