# Validation record

## Version 1.7 update

Reviewed shared histogram edges, per-replica normalization, circular wrapping,
finite-sample counts, profile/categorical separation, presentation labels and
HTML navigation. Python and report JavaScript syntax were checked. The synthetic
example was regenerated without Matplotlib. No test suite, native CPPTRAJ run,
browser/WebGL check or Matplotlib figure rendering was performed; visual appearance
remains unverified. No packages were installed.

## Version 1.6 update

Code review covered dropdown grouping, bounded PDB embedding, optional native
reference export and local third-party asset packaging. The synthetic report
example was regenerated with NumPy and a separately attributed experimental
structure. This is artifact generation, not validation of molecular calculations.
No test suite, native CPPTRAJ run, Matplotlib rendering or browser/WebGL check was
performed. Browser appearance remains unverified because the available browser
tool rejected local-file access. No packages were installed.

## Version 1.5 update

Native RDF/watershell commands were checked against official cpptraj documentation
and watershell source. Configuration, beginner flow, terminal output and figure-style
changes were reviewed by reading the code. No tests, native solvation calculations,
Matplotlib rendering or browser checks were run, and no dependencies were installed.
Earlier synthetic demonstrations do not validate the new figure rendering or analyses.

## Version 1.4 update

A synthetic report was generated with publication plotting disabled, exercising
report assembly and embedded numerical payload generation without Matplotlib or
CPPTRAJ. The synthetic example includes continuous, circular and categorical
observables and is explicitly labelled as artificial data. Browser visual review
was attempted but blocked by the browser's local-file URL security policy; the
JavaScript charts and interactions have not been visually verified. Code and setup
changes were reviewed statically. No test suite or native calculation was run, and
no dependencies were installed. Synthetic generation is not scientific validation.

## Version 1.3 update

The multi-page HTML report, evidence-linked review rules and semantic figure
indexing were inspected by reading the code. No tests, calculations, demos or
browser/rendering checks were run. No dependencies were installed. Earlier
validation results do not validate the new report pages.

## Version 1.2 update

Literature-informed diagnostics, diagnostic figures, offline HTML reporting,
environment/configuration checks, result bundling and methods records were reviewed
by reading the code only. No tests, calculations, demos or rendering checks were
run, and no dependencies were installed. Historical results below do not validate
these new features. Contributor guidance lists the remaining acceptance checks.

## Version 1.1 update

The preset/add-on menu, custom geometry monitors, configuration validation, and
terminal progress changes were reviewed by reading the code. No tests, demo,
CPPTRAJ calculation, or rendering checks were run for this update, at the user's
request. No dependencies were installed. The results below describe version 1.0
only and must not be treated as validation of version 1.1.

## Version 1.0 historical results

Validation date: 2026-09-05.

The final Python test suite ran **88 tests: 87 passed and 1 was skipped**.
Run the same suite from the source directory with:

```text
python -m unittest discover -s tests -v
```

The tests cover fixed-width Amber metadata, modified/terminal residues,
selection safety, file ordering and duplicate detection, global frame windows
across chunks, common-reference construction, resource guards, failure handling,
pooled projection splitting, numerical/circular/categorical summaries, contact
and nucleic parameter tables, compatible matrix axes, and the interactive wizard.

Workflow smoke tests use a clearly identified simulated CPPTRAJ subprocess.
They execute the actual Python orchestration/reporting code with synthetic
outputs. A separate three-replica synthetic demo exercises unequal lengths.
All three example JSON configurations passed schema/path normalization with
placeholder file existence checks disabled.

**Not validated in this environment:** real CPPTRAJ execution on molecular
trajectories, native command compatibility with a particular installed version,
and the appearance of rendered SVG files. CPPTRAJ and Matplotlib were absent;
the Matplotlib rendering smoke test was skipped explicitly. No dependencies
were installed. The implementation was checked against official CPPTRAJ
documentation and selected upstream source for output formats and commands.

The packaged `.pyz` is a standard Python executable ZIP containing the same
application modules. Its CLI/version and numerical-demo entry points were
smoke-tested without installing it. It still requires NumPy for numerical
reporting, Matplotlib for figures, and CPPTRAJ for trajectory calculations.

For an initial real-system check, select a short frame window, run the basic
preset, and inspect the generated input, log, reference, selected atoms and
numerical results before expanding the requested analyses. Verify complex
imaging and atom identities with your own topology and trajectories.

## Convergence and early trajectory-length checks (2026-09-08)

Executed using the existing `mace-agent` Python environment (NumPy 2.3.5 and
Matplotlib) with `MPLCONFIGDIR=/tmp/mdwb-mpl`:

- 18 convergence tests passed: stationary/shifted/drifting signals, correlated
  sampling, equal means with different distributions, equal-replica references,
  unequal durations, ns tails, custom/file targets and units, circular/categorical
  handling, missing/irregular data, resource caps, geometry cluster planning, and
  HTML/SVG reporting integration.
- 7 trajectory-length tests passed: joined segment counts, malformed metadata,
  warning before range prompts, ns/frame trimming, shortest endpoint, unavailable
  engine disclosure and length probing before atom-selection validation.
- Full suite: 113 tests, 3 failures, 1 skip. The same three failures were reproduced
  in a clean `git archive HEAD` baseline (88 tests, 3 failures, 1 skip):
  `test_wizard_standard_defaults_save_and_plan_without_cpptraj`,
  `test_svg_smoke_and_appended_warning`, and
  `test_basic_workflow_generates_manifest_and_joint_numeric_tables`.
  They concern the existing basic/standard wizard default, an SVG wording
  expectation, and reference-preparation script count respectively.
- Python compilation and `git diff --check` passed. The portable archive's Python
  modules were verified byte-for-byte against source; its dependency-free
  `python3 mdworkbench.pyz catalog` command succeeded.

CPPTRAJ is unavailable in this environment. Native length-reader behavior and
new clustering commands have not been executed against a real MD trajectory.
The length-reader tests use controlled CPPTRAJ output; cluster tests verify
commands and shared-label planning, not native clustering results.

## Pi-stacking workflow (2026-09-08)

- All 22 new stacking tests passed in the existing NumPy environment. Coverage
  includes connectivity-based discovery with arbitrary atom labels, fused rings,
  saturated-carbon exclusion, lone-pair oxygen rings, guanidinium candidates,
  PDB inference and MOL2 atom ordering, all-pair planning/batching, sign-independent
  normal angles, boundary cutoffs, equal-replica means and sample SD, any-partner
  unions, missing/invalid geometry, wizard confirmation/naming, and report/HTML/
  convergence integration. Replot tests reject changed centers/atom membership.
- Full suite: 135 tests, the same 3 previously documented failures, and 1 skip.
  No additional failures were introduced by this feature.
- Python compilation and whitespace validation passed. The portable zipapp was
  refreshed and its Python entries checked against the current source files.
- Geometry commands were reviewed against official CPPTRAJ vector/vectormath
  documentation. CPPTRAJ remains unavailable here: native trajectory execution,
  periodic imaging behavior and ring-plane calculation require real-data validation.
  Synthetic geometry tables are not evidence of a native MD calculation.

## CPPTRAJ paths containing spaces

Fixed the user-reported `good dye` failure during the initial trajectory-length
check. CPPTRAJ reconstructs and re-tokenizes argv, and parses `-y` again as a
trajectory expression. File arguments with whitespace/quotes now use temporary,
relative, space-free symlinks inside the process working directory. No trajectory
copy or input modification is performed. Original path mappings are written in
logs, and failed wizard length-check logs are retained instead of deleted.

Four process-boundary regression tests pass, simulating both native parsing stages
and covering spaced/apostrophe paths, script loading, masks, cleanup and failure
log retention. Full suite: 139 tests, the same three known failures and one skip.
This is not native validation on Dawn; replace the remote zipapp and retry there.
Source inspection: https://raw.githubusercontent.com/Amber-MD/cpptraj/master/src/Cpptraj.cpp
and https://raw.githubusercontent.com/Amber-MD/cpptraj/master/src/CpptrajState.cpp.

## Consistent bare residue-range selections

Bare residue selections such as `6-8` and `6-8,12` now use one shared normalizer
in wizard prompts, saved-configuration loading, selection resolution and direct
analysis planning. Alignment, advanced/PCA/clustering, imaging, region, interface,
geometry-monitor and solvation masks all accept the shorthand. Explicit atom masks
such as `@6-8` retain their meaning. Malformed/descending bare ranges are rejected.

Six new regression tests pass, including an advanced plan from a saved config
containing bare ranges. Full suite: 145 tests, the same three known failures and
one skip. The portable zipapp includes the new selections module.

## Interactive advanced resource limits

Both setup paths now expose pooled-frame, selected-atom and estimated-memory
limits after analysis selection. Previously these were JSON-only settings.
Known trajectory counts are reused from initial loading, with joined-segment
trimming and stride applied exactly. Above the cap, users can raise it or change
the shared stride. Memory validation remains separate and unchanged.

Six new tests pass: raising the cap for 60,000 pooled frames, per-replica stride
rounding, trimming across chunks, unknown lengths and invalid memory input,
analysis-specific prompts, and saving limits through both wizard paths.
Full suite: 151 tests, the same three known failures and one skip. CPPTRAJ is
unavailable here; native execution on Dawn has not been validated.

## Preserve full sampling when configuring resource limits

The resource prompt no longer offers to change global stride. It suggests a
frame cap covering the selected total, preserving all sampling settings.
Clustering continues to use its separate random sieve and assigns remaining
frames afterward. Help and runtime errors now direct users to resource budgets
and the clustering sieve, instead of suggesting global thinning for PCA limits.
Existing configurations are not silently changed; reset frames.stride to 1 if a
previous setup increased it unintentionally. No video-generation workflow exists.

Seven resource tests pass, including preservation of stride in mixed analyses,
unchanged PCA batches when cluster sieve changes, and both wizard paths. Full
suite: 152 tests, the same three known failures and one skip. Portable sources
were rebuilt; native CPPTRAJ execution remains unverified here.

## Readable SVG exports and lag explanation

Scalar plots now include readable metric/region labels and plot types, with
replica context on individual exports. Hash-only convergence, autocorrelation
and between-replica distribution-distance names are replaced by descriptive
names with short collision-resistant suffixes. Numerical identifiers remain
unchanged. The sampling review explains lag, time/stride conversion, slow
correlation decay and the distinction from convergence and equilibration.

Two new tests cover sanitized-label collisions and real SVG/report generation,
including figure indexing and HTML lag help. Full suite: 154 tests, the same
three known failures and one skip. Rebuilt the portable archive. Existing
reports can be regenerated without native trajectory analysis; old SVG files
are left intact while report links point to current exports.
