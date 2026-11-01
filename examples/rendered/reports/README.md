# Analysis reporting notes

- Trajectory-frame SD, quantiles, histograms and running averages are descriptive, not uncertainty estimates.
- Between-replica SD describes the spread of replica means; each available replica receives equal weight. One replica has no between-replica SD.
- Angles flagged circular use circular means and resultant lengths, with means undefined when the resultant is zero. They are excluded from arithmetic means and ordinary Pearson correlations.
- Block-mean SD uses complete contiguous blocks of the requested analyzed-frame size; it is not an effective sample-size estimate or a confidence interval. Choose blocks longer than the observable correlation time.
- Independent replicas are appended only for visualization, with boundary markers; the appended axis is not a continuous physical trajectory.
- Matrix and profile means require matching axes and shapes. Matrix cell means require finite values from every included replica.
- Matrix index compatibility does not independently establish identical atom identities or alignment; consistent selection, topology and reference setup are required.
- Replica groups with different lengths are overlaid on their own coordinates. Missing observations are not interpolated.
- Unknown table formats are retained as original artifacts and indexed. Numeric generic tables are exported by row without assigning physical meaning.
- POOLED entries label shared analysis outputs, not additional independent replicas. Pooled PCA covariance weights frames; longer replicas contribute more to the common basis. PCA scatterplots use reported projection coordinates without inferring explained variance or free energy.
- DSSP integer states are categories. Their fractions are reported; arithmetic means of state IDs are not calculated.

## Warnings

None.

## Main outputs

- all_replicas.dat: one tidy TSV with replica labels and source paths.
- descriptive_statistics.dat: per-replica descriptive statistics and contiguous block diagnostics.
- replicate_statistics.dat: equal-weight replica mean and between-replica SD.
- circular_statistics.dat and circular_replicate_statistics.dat: angular descriptions in degrees.
- categorical_fractions.dat: state fractions within each replica.
- combined/: aligned per-observable tables; blank values represent missing coordinates.
- probability_distributions.dat: common bin edges, counts, probabilities and densities per replica; scalar time series only.
- artifact_index.dat: every requested artifact, its source, and parse status.
- ../figures/: editable-text SVGs; dense matrix cells and PCA points are embedded rasters at the configured DPI.
