# MD-bud 2.0.0 validation

Observed on 1 November 2026. Linux is the primary deployment target.

| Check | Observed result |
|---|---|
| Ubuntu/WSL, Python 3.14.4, NumPy 2.5.4, Matplotlib 3.11.2 | **198 tests, successful; 1 Windows-only test skipped** |
| Windows, Python 3.12.14 with numerical/plot dependencies | **198 tests, successful; 19 POSIX campaign tests skipped** |
| Linux Bash/Slurm campaign regressions | Passed: script syntax, dependency order, duplicate prevention, engine failure propagation, spool working directories, reset/cancellation/archival recovery and frame parsing, using process doubles |
| Campaign to analysis | Passed: production-only mapping, replica/segment order, saved-frame interval, trim/stride, incomplete stage/frame rejection, checksums, relocated configurations and synthetic end-to-end HTML/numerical reporting |
| Protocol generation | Passed: explicit MC NPT production barostat, restrained density relaxation, final unrestrained step, uniform output spacing and LF shell scripts |
| Production-log QC | Passed: temperature/energy/density parsing, Fortran exponents, exclusion of average blocks, invalid/nonfinite values and engine errors |
| Release build | Linux wheel and source distribution built with setuptools 84.2.1; unified zipapp built with the standard library |
| Installed wheel away from checkout | Passed: version, dependency check, synthetic offline report/template/assets and nonzero failure exit code |
| Portable zipapp away from checkout | Passed: same smoke checks without a source checkout; library dependencies supplied by the environment |
| Generation without numerical packages | Passed via Python -S, both source and portable archive |

The synthetic report smoke test runs from a temporary directory. It verifies
that packages, report templates and viewer license/assets survive packaging.
It is not a molecular calculation. Plotting regressions generate SVG files and
check labels, output identities and replica-boundary warnings; they are not a
full browser/WebGL or visual-layout review.

The final local commands included:

```bash
python -m unittest discover -s tests -v
python -m build --no-isolation
python tools/build_zipapp.py --output dist/md-bud-2.0.0.pyz
# After installing the built wheel into the isolated test environment:
python tools/smoke_release.py --archive dist/md-bud-2.0.0.pyz
```

Local test environments and logs are excluded from Git. The included GitHub
Actions matrix targets Python 3.11/3.12 on Linux and Windows; it has not been
executed on hosted GitHub by this task. Local Python 3.14 results do not replace
the minimum-version CI checks.

## Remaining acceptance work

No real AMBER simulation, native CPPTRAJ calculation, live Slurm submission or
GPU run was performed. Bash/Slurm/AMBER process doubles test orchestration and
failure behavior, not scientific correctness or cluster module compatibility.
The test environment has no native AMBER/CPPTRAJ engines. Resource defaults
and the new pressure-control/equilibration inputs require validation on the
actual cluster before production use.

Recommended site acceptance: run a small prepared periodic protein system
through all stages, check energy/temperature/density and actual frame counts,
compare basic CPPTRAJ results to a reviewed manual calculation, exercise two
replicas and two segments, and test reset/recovery in a disposable campaign.
Then repeat representative nucleic/mixed systems. Review scientific applicability
and remaining work in [SCIENTIFIC_REVIEW.md](SCIENTIFIC_REVIEW.md).

The imported analyser's historical validation record is preserved separately in
[ANALYSIS_VALIDATION_HISTORY.md](ANALYSIS_VALIDATION_HISTORY.md).
