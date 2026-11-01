# Exploratory sampling diagnostics

ACFs are computed separately for finite, regularly sampled, nonconstant linear observables with at least 32 samples. Missing values are not removed to manufacture contiguous sampling. Circular and categorical series are excluded from ESS. FFT autocovariances use a fixed N denominator. Initial-positive monotone pairs start at lags 0/1, 2/3, etc.; g = max(1, -1 + 2 sum(pairs)), capped at N. ESS=N/g and SEM=sqrt(sample variance*g/N) assume stationarity. Neither is proof of equilibration or a calibrated confidence interval. Truncation at the lag limit needs review. The suggested block size ceil(5*g) is a heuristic; no configuration or data is changed.

Jensen-Shannon DISTANCE uses base-2 logs and is the square root of the divergence (0 to 1). Each replica histogram is normalized independently; common bins are retained alongside results. Circular bins cover [-180,180); categories use shared states. Values are descriptive, sensitive to binning and sample size, and are not significance tests or convergence thresholds. Identical sampled distributions do not establish adequate exploration.

Methods: https://www.stat.umn.edu/geyer/mcmc/library/mcmc/html/initseq.html
Context: https://pymbar.readthedocs.io/en/latest/timeseries.html and https://github.com/drorlab/pensa
