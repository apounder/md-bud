# Literature and software review

Reviewed 5 September 2026. This is a focused comparison of primary papers, official
documentation and repositories, not a systematic review or a performance benchmark.
The workbench remains a CPPTRAJ driver; it does not execute the reviewed libraries.

## What established projects contribute

| Project and primary references | Relevant strength | Workbench response in 1.2 |
| --- | --- | --- |
| [MDAnalysis](https://github.com/MDAnalysis/mdanalysis), [Michaud-Agrawal et al., 2011](https://doi.org/10.1002/jcc.21787), [Gowers et al., 2016](https://doi.org/10.25080/Majora-629e541a-00e) | Reusable trajectory/selection abstractions and a broad analysis ecosystem | Preserve explicit topology and atom-order requirements; improve offline configuration validation and document interoperability boundaries. A filename cannot establish correspondence. |
| [MDTraj](https://github.com/mdtraj/mdtraj), [McGibbon et al., 2015](https://doi.org/10.1016/j.bpj.2015.08.015) | Accessible trajectory operations and efficient numerical analysis | Keep reusable configurations and independent analysis batches. Format breadth remains dependent on the user's CPPTRAJ build; no universal-format guarantee. |
| [mdacli](https://github.com/MDAnalysis/mdacli), [official introduction](https://www.mdanalysis.org/2021/12/01/mdacli/) | Makes MDAnalysis analyses accessible through a CLI | Add discoverable `catalog`, `doctor` and `validate` commands alongside the existing guided wizard. |
| [MDAKits requirements](https://mdakits.mdanalysis.org/about.html) | Documentation, identifiable maintainers, licensing, versioning and regression tests make research software reusable | Add group-use guidance, contribution instructions, release notes and license. Existing tests remain supplied; this update has not been tested. This CPPTRAJ package is not an MDAKit and is not registered. |
| [PyMBAR](https://github.com/choderalab/pymbar), [time-series documentation](https://pymbar.readthedocs.io/en/latest/timeseries.html) | Correlated samples require more care than frame-count-based errors; equilibration selection is a separate decision | Add independent NumPy ACF, statistical inefficiency, effective sample-size and conditional SEM estimates per replica. No automatic equilibration removal. These are not PyMBAR outputs. |
| [PENSA](https://github.com/drorlab/pensa), [Vögele et al., J. Chem. Phys. 162, 014101 (2025)](https://doi.org/10.1063/5.0235544) | Systematic ensemble comparisons using feature distributions and collective descriptions | Add common-bin replica distributions and pairwise Jensen–Shannon distances with heatmaps. Existing pooled PCA uses a shared basis. PENSA's information-flow and advanced feature methods are not reproduced. |
| [ProLIF](https://github.com/cbouy/ProLIF), [Bouysset & Fiorucci, 2021](https://doi.org/10.1186/s13321-021-00548-6) | Chemistry-aware interaction fingerprints | Defer an optional adapter. Geometric contacts alone cannot justify labels such as aromatic stacking or hydrophobic interaction without chemistry-aware definitions. |
| [alchemlyb](https://github.com/alchemistry/alchemlyb) | Dedicated alchemical free-energy analysis | Defer integration: state energies, thermodynamic states and estimator assumptions are needed. Coordinate trajectories alone do not supply these. |

## Implemented scientific improvements

**Sampling diagnostics:** finite, regularly sampled linear observables with at least
32 observations receive an FFT autocorrelation calculation. The estimator uses a
fixed-N autocovariance denominator and initial-positive, monotone paired lag sums,
following the [Geyer initial-sequence description](https://www.stat.umn.edu/geyer/mcmc/library/mcmc/html/initseq.html)
([1992 paper](https://doi.org/10.1214/ss/1177011137)). The underlying stationary,
reversible sampling assumptions need scientific review for the trajectory/observable.
Statistical inefficiency is conservatively restricted to [1,N]; ESS=N/g and
SEM=sqrt(sample variance*g/N). No negative-correlation sampling benefit is claimed.
Lag-limit truncation and a heuristic block-size recommendation are recorded. Missing
values are not compacted, and circular/categorical signals do not receive these ESS
estimates. No equilibration time, p-value or confidence interval is inferred.

**Distribution comparison:** independently normalized histograms share a common
bin grid across replicas. Linear values use the observed pooled range, circular
angles use [-180,180), and categorical signals use shared states. Reported values
are the square root of base-2 Jensen–Shannon divergence, bounded by zero and one.
This is a distance, not divergence, and not a significance test. Binning and finite
sampling affect the result; circular bins also depend on the selected angular origin.
Zero distance cannot establish adequate conformational exploration. No distance
matrix is generated with fewer than two comparable replicas.

**Research-group usability:** the offline HTML viewer links SVGs, raw data and logs;
run-specific methods/environment records support review; a results-only ZIP can be
shared and its reports regenerated without the original simulation. These usability
choices are our synthesis of reproducible workflow needs, not methods claimed by
any one paper. Original input paths remain in provenance; bundling is not anonymization.

## Useful next extensions, deliberately not implemented

1. Chemistry-aware ProLIF fingerprints through an optional adapter with explicit
   hydrogen, protonation and bond-order handling.
2. Comparative studies across mutants/ligands using an explicit atom/residue map
   and condition labels. Current runs compare replicas sharing one topology.
3. Contact/H-bond lifetimes with definitions for continuous versus intermittent
   survival, censoring and periodicity; simple occupancy is not a lifetime.
4. Enhanced-sampling reweighting, free energies and kinetic models with supplied
   weights/state information and method-specific validation. Histogram populations
   do not automatically give unbiased free energies or kinetic rates.
5. Streaming reporting for very large datasets and restartable batches guarded
   by complete input/configuration fingerprints. Current diagnostic limits bound
   additional work, but the base parser still loads individual data tables in memory.

Before public release, identify actual authors/maintainers, add the real repository
URL, run the supplied regression suite and new numerical checks, and validate
representative native CPPTRAJ runs and SVG rendering. No registry submission,
publication, dependency installation or execution validation was performed here.
