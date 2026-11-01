# Computational chemistry review

The combined package covers a useful conventional AMBER workflow. It does not
cover all molecular dynamics. The priorities below are engineering/scientific
judgment informed by the linked primary sources, not assertions that every item
has been implemented or scientifically validated.

## What the merge now handles

| Need | Implemented behavior | Remaining responsibility |
|---|---|---|
| Consistent simulation → analysis handoff | Copied-topology checksums, explicit production segment lists, replica identities, saved-frame interval from `dt × ntwx` | Verify physical model, atom identities and native trajectory integrity |
| Ordinary periodic solution simulations | Minimize, heat, restrained equilibration, final unrestrained equilibration, NPT/NVT production | Assess density/energy/geometry relaxation and extend preparation as needed |
| Production ensemble | Explicit Monte Carlo barostat for new NPT production; Langevin thermostat and constrained 2 fs protocol | Verify engine support and suitability for the system |
| Multiple runs | Distinct saved seeds, separate directories and ordered continuations | Independently equilibrate starts; seed differences alone do not prove independence |
| Structural interpretation | Shared fit/reference, replica-aware tables, circular angular statistics, distribution and correlation diagnostics | Choose selections/references consistent with the scientific question |
| Numerical/log checks | `frames` checks counts; `qc` summarizes production temperature/energy/pressure/density and flags errors | Inspect native logs; QC does not establish equilibrium, stationarity or correct ensemble |
| Sampling | Descriptive convergence curves, ESS diagnostics and replica comparisons inherited from the analyser | Quantify uncertainty per observable; finite stable traces cannot prove convergence |
| Reproducibility | Simulation seeds/inputs, analysis configuration, logs, methods record, campaign-manifest hash | Preserve force-field versions, preparation decisions and site/engine environment |

The production-barostat change avoids silently relying on AMBER's inherited
weak-coupling default. The primary [stochastic cell-rescaling paper](https://arxiv.org/abs/2006.09250)
discusses weak coupling's incorrect volume fluctuations. The
[MC barostat membrane study](https://pmc.ncbi.nlm.nih.gov/articles/PMC8776593/)
also shows why an NPT choice suitable for solution simulations cannot be applied
uncritically to membranes. Implementation options should be checked against the
[AMBER manual](https://ambermd.org/doc12/Amber24.pdf) and the installed engine.

## Highest-priority gaps before broad scientific adoption

1. **Native cluster acceptance tests.** Run a small periodic solvated protein and
   a nucleic-acid system using the actual site modules, GPU/CPU engines and scheduler.
   Compare generated inputs, measured frames and basic observables against reviewed
   manual CPPTRAJ commands. Offline process doubles cannot replace this.
2. **Preparation validation.** Add robust topology/restart/box/velocity validation,
   charge/ion inventories and clear provenance for protonation, waters, metals and
   parameterization. The generator currently checks that inputs exist and saves
   checksums; that does not establish compatibility or chemical correctness.
3. **Preparation diagnostics.** Extend log inspection to all min/heat/eq substages,
   density relaxation, restraint-release effects and study-specific stationarity
   assessments. Do not make a universal automatic 'equilibrated' switch. Production
   currently proceeds after successful preparation jobs, without an adaptive decision.
4. **Real scientific regression datasets.** Add small redistributable periodic
   datasets and native expected results, including multichain imaging, modified
   nucleotides, ligands and production restarts. Current regression fixtures are
   deliberately synthetic.

These priorities follow the need for physical preparation and observable-specific
sampling assessment described in the primary
[foundations guide](https://pmc.ncbi.nlm.nih.gov/articles/PMC6884151/) and
[uncertainty/sampling guide](https://pmc.ncbi.nlm.nih.gov/articles/PMC6286151/).

## Further scope that needs dedicated protocols

| System/question | Why a dedicated workflow is needed |
|---|---|
| Membranes/interfaces | Semi-isotropic pressure control, surface tension and suitable nonbonded treatment; isotropic solution defaults are insufficient |
| Implicit solvent, vacuum, coarse-grained or polarizable models | Different topology/constraint/electrostatics/integrator assumptions |
| Metals, covalent ligands and modified residues | Parameterization and chemical identity must be explicit; generic residue detection is insufficient |
| Free energies, umbrellas, metadynamics, alchemical changes | Bias/thermodynamic-state metadata, overlap, reweighting and uncertainty estimators |
| Replica exchange | State demultiplexing and ensemble-specific analysis before comparing trajectories |
| Diffusion, kinetics and residence times | Unwrapping, suitable output frequency and contiguous physical time; never cross independent-replica boundaries |
| QM/MM reactions | Reaction-coordinate and electronic-structure provenance; geometry monitors are useful but do not implement a QM/MM protocol |
| Ligand binding and biochemical function | System-specific contacts, chemical observables and experimental comparisons beyond RMSD |

## Defaults should expose scientific choices

The basic analysis preset is intentionally inexpensive and interpretable.
There is no automatic production discard: choose it after reviewing relaxation.
`--discard-ns X` removes saved frames at or before X ns from each replica;
`--stride N` then samples every Nth remaining frame. Sampling diagnostics use the
selected data, so discarding/downsampling must remain part of the methods record.
The first stored production frame is ordinarily one output interval after the
production start; display coordinates may use elapsed analyzed time from zero.

A 50 ps output interval is adequate for many slow structural summaries but loses
fast motions. Choose the interval around the observable before running, since
information not saved cannot be recovered later. Keep the timestep at or below
2 fs for this non-HMR protocol; longer steps require an explicitly different model.

The [CPPTRAJ `trajin` documentation](https://amberhub.chpc.utah.edu/trajin/) defines
frame selection; the [imaging documentation](https://amberhub.chpc.utah.edu/autoimage/)
explains why condensed/multicomponent systems may need a custom anchor. The
package retains its existing shared-reference and replica-boundary protections.
Unknown time intervals in generic projects stay unknown rather than being inferred.

Sources reviewed on 1 November 2026. The linked methods support the rationale;
they do not certify this implementation or any generated simulation.
