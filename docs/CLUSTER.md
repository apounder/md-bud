# Linux cluster operation

Linux is the deployment target. Use a Python 3.11+ environment on both login
and compute nodes. Slurm management uses `sbatch`, `squeue`, `sacct`, `scancel`
and `srun`; generated scripts use Bash, `realpath` and `flock`. Preparation uses
`pmemd.MPI`, while production defaults to `pmemd.cuda` (`--engine cpu` selects MPI).

## Install once in your own environment

Use your site's documented Python/module environment, then:

```bash
python3 -m venv ~/venvs/md-bud
source ~/venvs/md-bud/bin/activate
python3 -m pip install '.[plots]'
md-bud configure cluster.json --account your-allocation --modules "YOUR MODULE LIST"
md-bud doctor --task cluster
md-bud doctor --task analysis
```

Cluster diagnostics check MPI preparation and CUDA production by default.
Use `doctor --task cluster --engine cpu` when planning CPU-only production.
These checks inspect PATH; they do not submit jobs or test module compatibility.
Run diagnostics after loading the same software stack you will use on compute nodes.

The default module list was inherited from the Fir project:
`StdEnv/2023 gcc/12.3 openmpi/4.1.5 cuda/12.2 amber-pmemd/24.3 ambertools/25.0`.
Availability is site-specific. Use `--modules ""` when inheriting an already loaded
stack, or save the site's modules in your profile. The package neither loads nor
installs modules during generation. Generated jobs load the configured stack.

## One-file deployment

Build once from this checkout:

```bash
python3 tools/build_zipapp.py --output dist/md-bud.pyz
```

Copy `md-bud.pyz` to the cluster. Every command works with
`python3 /path/to/md-bud.pyz` in place of `md-bud`. The archive includes both
engines, HTML assets and license notices. It does not bundle Python, numerical
libraries, AMBER, CPPTRAJ or Slurm. Generation uses only Python; analysis requires
NumPy and optional Matplotlib in the active environment. No network is needed
after dependencies are available.

## Run analysis on a compute node

Use an interactive allocation according to your site's rules, or adapt this
ordinary Slurm script. Replace the account, modules, virtual environment, archive
and campaign paths with your own. Do not run expensive trajectory analyses on a
login node.

```bash
#!/bin/bash
#SBATCH --job-name=md-bud-analysis
#SBATCH --account=your-allocation
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=md-bud-analysis-%j.log
set -euo pipefail
# Load your site's compatible Python and AmberTools modules here.
source "$HOME/venvs/md-bud/bin/activate"
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
md-bud analyze /absolute/path/to/runs/study --output /absolute/path/to/results
```

Submit it with `sbatch analyze.sh` once production finishes. To chain it yourself,
use `--dependency=afterok:ID1:ID2:ID3` with the final production job ID from **every
replica**. MD-bud's normal campaign submission does not automatically schedule
the report job. The separate `combine` job remains optional and is unnecessary
for analysis. Choose RAM and walltime based on system/frame count; pooled matrix
analyses have independent resource guards and can still be expensive.

## Completion, reset and recovery

Campaign submission records job IDs before proceeding, refuses ambiguous states
and uses success dependencies. Reset defaults to preview; applied reset requires
inactive jobs or explicit cancellation and archives affected outputs. A reset
reruns a whole stage from its upstream restart. It does not continue an interrupted
segment in place. Do not edit a generated campaign's copied inputs or mix partial
segments with replacement segments.

`status --offline` inspects local state when scheduler access is unavailable.
`frames` checks expected production-frame counts. `qc` checks readable production
logs. Run both before analysis; stage completion is a file-based execution marker,
not a scientific suitability verdict. Retain Slurm and AMBER logs as well as the
input/analysis manifests.

Copy the complete analysis output folder back to a workstation. `index.html`
opens without a server. `md-bud bundle RESULTS --destination study.zip` makes a
portable results archive. The offline viewer shows a reference structure, not an
entire trajectory animation.
