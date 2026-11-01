# Scientific interpretation and scope

> These imported analysis notes retain historical source-release details. For the
> combined MD-bud package, see [current validation](VALIDATION.md) and the
> [computational chemistry review](SCIENTIFIC_REVIEW.md).

This workbench automates common structural MD measurements. It cannot establish that a simulation is physically valid, equilibrated, converged, or chemically comparable to another simulation. The recommendations below are analysis policy inferred from the cited methods and CPPTRAJ documentation; they are not claims that every feature is implemented. Consult the main README for implemented options.

## File discovery is not scientific compatibility

CPPTRAJ reads Amber NetCDF, Amber coordinate/restart files, CHARMM/NAMD DCD, GROMACS XTC/TRR, PDB and several other formats. Actual availability depends on the executable and how it was built. A filename extension is a discovery hint; CPPTRAJ performs format detection, and `trajin ... as FORMAT` can override it. Every trajectory must be associated with a suitable topology. A readable format does not imply that a requested analysis is meaningful for it. A PDB-only topology may lack reliable bonds, atom types or radii required by some calculations. See the official [format table](https://amberhub.chpc.utah.edu/cpptraj/trajectory-file-commands/) and [`trajin` documentation](https://amberhub.chpc.utah.edu/trajin/).

For combined analysis, use the same topology and atom order for all replicas. Matching atom counts are necessary but insufficient: coordinates can have the correct count and the wrong identity/order. Stripped and solvated trajectories generally need different topologies. Mutants, different protonation states, changed residue numbering and different parameterizations belong in separate projects unless a deliberate atom/residue correspondence has been established. Ordinary `.nc` files do not reliably encode that correspondence. CPPTRAJ's [RMSD tutorial](https://amberhub.chpc.utah.edu/amber-hub/start-here-rmsd-analysis-in-cpptraj/) explicitly requires matching atom counts and ordering.

## Replicas, segments and time

- **Independent replicas:** separate simulations used as independent observations. Give each its own replica name.
- **Sequential segments:** continuation files from one simulation. Put their paths, in temporal order, in the same replica. Check overlap, missing chunks and duplicate restart frames.
- **Replica exchange or biased sampling:** an ordinary collection of input files is not sufficient. State demultiplexing, ensemble selection and possibly reweighting must precede interpretation. CPPTRAJ has dedicated `remdtraj` options, described in [`trajin`](https://amberhub.chpc.utah.edu/trajin/).

Natural filename sorting cannot establish any of these scientific relationships. Concatenating independent replicas is useful for distributions or a marked display, but the boundary is a discontinuity. Never estimate diffusion, transition rates, lifetimes, autocorrelation or time-lagged models across that boundary.

Do not infer physical time from frame count, filenames or a nominal integration timestep. The saved-frame interval includes the output frequency. If the interval is unknown, use frame/sample indices. If an interval is declared, account for selection stride and make clear whether the origin is simulation time or elapsed analysis time. A NetCDF `time` variable may be missing, empty or restart at segment boundaries; the official [RMSD example](https://amberhub.chpc.utah.edu/amber-hub/start-here-rmsd-analysis-in-cpptraj/) itself contains an empty time variable.

In this workbench, `frames.start`, `frames.stop` and `frames.stride` apply to the ordered, concatenated source frames **within each replica**, then restart for the next independent replica. They do not restart at each segment. The runner obtains segment lengths from CPPTRAJ and translates this global selection into per-file windows; the stride phase is preserved across segment boundaries. `dt_ps` is the interval between source saved frames, before applying this stride.

## Residues and automatic selections

Ordinary masks such as `:10-35` refer to **one-based residue indices in the loaded topology**, not automatically a paper's numbering or a PDB's original residue IDs. `:;` selects preserved original numbering where supported; chain masks such as `::A` require chain metadata. Check the residue inventory and resolved selections before using a biological interpretation.

Examples: `:10-35` selects a region; `:10-35@N,CA,C,O` selects its protein backbone; `:10-35&!@H=` selects names not beginning with H. Name-based hydrogen exclusions can miss unusual naming; elemental selection is preferable when element metadata is reliable. Restrict atom-name backbone masks to identified protein/nucleic-acid residues. Standard residue-name detection is a convenience, and modified residues, ligands and coarse-grained models need custom masks. Distance masks are generally evaluated against a reference, not updated each frame, except in the `mask` action. See [mask syntax](https://amberhub.chpc.utah.edu/atom-mask-selection-syntax/).

## Imaging and alignment

For periodic structural analysis, inspect representative frames after imaging. `autoimage` groups coordinates by molecule and normally anchors the first molecule; a protein–DNA complex or multiple chains may require an explicit central anchor. Imaging requires meaningful unit-cell/molecule information. Applying a default anchor to every system is not universally correct. See [`autoimage`](https://amberhub.chpc.utah.edu/autoimage/).

Use a shared reference and alignment mask across replicas. An RMS fit transforms coordinates for subsequent actions. A section RMSD after the common fit with `nofit` retains domain or partner motion; fitting each section separately measures its internal deformation. Those answer different questions and should have different labels. RMSF and Cartesian correlation need a defined fit to remove overall translation/rotation. Diffusion instead requires a continuous, appropriately unwrapped trajectory, without a structural RMS fit. See [action ordering](https://amberhub.chpc.utah.edu/cpptraj/action-commands/).

## Choosing and comparing measurements

| Measurement | Interpretation to retain |
|---|---|
| RMSD, radius of gyration, distances | Selection, fit/reference, mass or geometric weighting, and distance imaging convention. A stable RMSD alone does not prove convergence. |
| Residue RMSF | Selection and common fit; one C-alpha or phosphate per residue differs from averaging atomic fluctuations over a residue. |
| Protein secondary structure | Categorical assignments and occupancy fractions; arithmetic means of structure codes have no physical meaning. See [`secstruct`](https://amberhub.chpc.utah.edu/secstruct/). |
| H bonds | Geometric criterion and donor/acceptor chemistry. CPPTRAJ defaults to donor–acceptor distance <3 Å and A–H–D angle >135°. Automatic F/O/N rules are approximate. Count both donor directions at a protein–DNA interface. Imaging is optional and off by default. See [`hbond`](https://amberhub.chpc.utah.edu/hbond/). |
| Contacts | Which atom pairs, cutoff, exclusion of close sequence neighbors, reference/native definition and normalization. A mean distance matrix is not a contact probability matrix. |
| Nucleic-acid geometry | Common residue-pair identities are essential. `nastruct` defaults to pairing from the first frame; replicas can otherwise assign different pairs to the same numbered column. Modified bases may require `resmap`/`baseref`. Pucker needs explicit dataset output; it is not included automatically in `naout`. See [`nastruct`](https://amberhub.chpc.utah.edu/nastruct/). |
| Torsions and pucker | Periodic angular data need circular statistics; 179° and −179° average near 180°, not 0°. |

## Matrices, PCA and clustering

For directly comparable PCA overlays, fit all selected coordinates consistently and project every replica onto **one common basis**, usually computed from pooled compatible frames. Independently calculated PC1/PC2 axes do not share a coordinate system; signs and nearly degenerate subspaces may differ. Preserve the mean structure, eigenvalues, eigenvectors, atom selection, reference and replica/frame mapping. A pooled basis weights replicas by contributed frames unless sampling is deliberately balanced. With N selected atoms, Cartesian covariance is 3N by 3N; rank is at most the number of frames minus one. Long trajectories and large selections can therefore be expensive and undersampled. PCA visual separation is not a kinetic model. See the official [combined DNA PCA tutorial](https://amber.utah.edu/AMBER-workshop/London-2015/pca/).

An equal-weight mean of per-replica DCCMs is a descriptive mean of correlations; it is not the DCCM obtained from pooled coordinates. The latter also reflects different replica means. Matrix differences require identical row/column identities. Use fixed −1 to +1 color limits for correlation matrices and a symmetric diverging scale for differences. Thresholded graph summaries depend on the chosen cutoff. Distance, correlation and covariance matrices are separate objects in [`matrix`](https://amberhub.chpc.utah.edu/matrix-2/).

Joint clustering supports shared cluster identities across replicas. Cluster 1 from two separately clustered trajectories need not represent the same structure. Representative structures and occupancy differences are useful; transitions over independent-replica boundaries are invalid.

## Uncertainty and publication figures

Frames are correlated observations. Prefer independent-replica summaries for uncertainty between simulations; when only one replica exists, use suitable block/autocorrelation diagnostics and acknowledge their assumptions. Do not calculate a naive standard error as frame SD divided by the square root of all frames. An SEM across replica means assumes independent, comparable replicas; with few replicas it is imprecise. Replica SD, SEM and confidence intervals must be labelled separately. Exclude equilibration by a scientifically justified rule and retain that choice in the configuration. Apparent precision does not account for force-field/model error. These recommendations follow [Grossfield et al., uncertainty and sampling guidance](https://pmc.ncbi.nlm.nih.gov/articles/PMC6286151/).

For overlays use matching limits, units and color identity. Mark boundaries in concatenated displays, retain individual replicas alongside any average, and state whether a pooled statistic gives equal weight to replicas or to frames. Smoothing is a display operation, not evidence of convergence. SVG output should retain editable text and explicit axes/units; a finished figure still needs system-specific biological labels and a caption documenting analysis choices.

## Integration-test data and local verification limits

The official [RMSD tutorial](https://amberhub.chpc.utah.edu/amber-hub/start-here-rmsd-analysis-in-cpptraj/) links a small trpzip2 topology and NetCDF trajectory (220 atoms, 13 residues, 1,201 frames). This is useful for an actual CPPTRAJ smoke test of RMSD/RMSF/radius of gyration, protein secondary structure and H bonds. Splitting its frames can test segment handling; duplicating the trajectory can test file aggregation but does **not** create independent simulations. It is nonperiodic and cannot test periodic imaging. The official [CPPTRAJ test suite](https://github.com/Amber-MD/cpptraj/tree/master/test) includes dedicated imaging, DSSP, matrix and nucleic-acid tests suitable for broader validation on a machine with CPPTRAJ.

During development in the supplied Windows workspace, `cpptraj`, `python` and `python3` were absent from PATH. WSL enumeration failed with `E_ACCESSDENIED`. A bundled Python runtime with NumPy was present, but initial imports found no Matplotlib, SciPy, netCDF4 or pytest. These observations establish the initial environment only; they do not validate generated CPPTRAJ commands. The parent README/test report records any subsequent checks and installation work. No successful live CPPTRAJ calculation is implied by a dry run or synthetic-data test.

Sources checked 5 September 2026. For an installed CPPTRAJ version, inspect `help COMMAND` and preserve the reported version in the run log.
