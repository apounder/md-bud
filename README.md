![md-bud banner](banner.png)

# md-bud

One script to generate and manage AMBER molecular dynamics jobs on Fir.

- Generate minimization, heating, equilibration, and production inputs.
- Run multiple replicas, each with its own files and random seeds.
- Chain production segments with Slurm dependencies.
- Combine trajectories separately for each replica.
- Check job status and trajectory frame counts.
- Reset a stage and its downstream jobs, archiving old results first.

Requires **Linux, Python 3.11+, Bash, and `flock`**. Generation works locally.
Submission needs Slurm; simulation needs AMBER; frame checking needs CPPTRAJ.
No Python packages need installing.

## Start here

```bash
chmod +x md-bud
./md-bud generate -i runs/example -p system.prmtop -s system.inpcrd \
  -A your-allocation -R 3 -n 5 -t 10 -w 0-03:00 -c

./md-bud submit runs/example --dry-run
./md-bud submit runs/example
./md-bud status runs/example
./md-bud frames runs/example
```

This makes three replicas, each with five sequential 10 ns production segments
(50 ns per replica). Each replica follows:

```text
min → heat → eq → md-0001 → md-0002 → … → combine (optional)
```

Replicas run independently. Each job waits for the previous job to succeed.
Generation does not submit unless you add `--submit`. Review generated `.in`
and `.sh` files before running them. Existing campaign directories are never
replaced.

## Most common options

These are the usual options to choose for a run. Short flags and long names do
the same thing; all existing long names still work.

| Short flag | Long name | Meaning | Default |
| --- | --- | --- | --- |
| `-i NAME` | `--output`, `--prefix` | Campaign prefix: the directory holding all run files | Required |
| `-p FILE` | `--topology` | AMBER `.prmtop` file | Required |
| `-s FILE` | `--coordinates` | Starting coordinates or production restart | Required |
| `-A NAME` | `--account` | Slurm allocation account | `SLURM_ACCOUNT`, otherwise required |
| `-n N` | `--segments`, `--chunks` | Sequential production chunks per replica | `1` |
| `-R N` | `--replicas` | Number of independent replicas | `1` |
| `-t NS` | `--ns` | Simulated time per chunk, in nanoseconds | `10` |
| `-w TIME` | `--walltime` | Maximum real runtime for each production job | `0-03:00` |
| `-f MODE` | `--workflow` | `all`, `prep`, or `production` | `all` |
| `-c` | `--combine` | Combine production chunks within each replica | Off |

For example, `-R 3 -n 5 -t 10` requests three replicas with five 10 ns chunks
each: **50 ns per replica**, or 150 ns across all replicas. `-w` controls cluster
walltime separately from `-t`.

`-n`, `-t`, `-w`, `-p`, and `-c` retain their original roles. `-i` now
names the campaign directory instead of prefixing each individual filename:
`-i water` puts predictable filenames inside `water/replica-001/`, etc.
`-R` is new: in `generate` it sets the replica count; in `submit`, `status`,
`reset`, and `frames` it selects a single replica number.

## Choose what to generate

```bash
# Minimization, heating, and equilibration only
./md-bud generate -i runs/prep -p system.prmtop \
  -s system.inpcrd -A your-allocation -f prep

# Production from an existing restart containing velocities and box information
./md-bud generate -i runs/production -p system.prmtop \
  -s equilibrated.rst7 -A your-allocation -f production -R 3 -n 5
```

## Reset and keep old data

```bash
# Preview: replica 2, segment 3 and everything after it
./md-bud reset runs/example -R 2 --from md-3

# Apply: cancel affected active jobs, wait for them to stop, then archive
./md-bud reset runs/example -R 2 --from md-3 --apply --cancel

# Submit the reset stages again
./md-bud submit runs/example -R 2
```

Reset keeps inputs, earlier stages, and other replicas. Old outputs, logs,
completion markers, and submission records go into `archive/TIMESTAMP/`.
`--from` accepts `min`, `heat`, `eq`, `production`, `md-N`, or `combine`.
Omit `--replica` to reset every replica. Nothing is deleted to reset a run.

A reset reruns the selected stage from its original upstream input. It does
**not** continue a partially finished segment from its intermediate checkpoint.
If scheduler state is unknown, reset stops instead of moving potentially active
files. An interrupted submission with an unknown job ID must be reconciled with
Slurm in `md-bud.json` before retrying.

## Commands

| Command | What it does |
| --- | --- |
| `generate` | Create a new campaign and copy its inputs |
| `submit PROJECT` | Submit new stages; keep completed or already submitted jobs |
| `status PROJECT` | Show local completion and Slurm states |
| `reset PROJECT --from STAGE` | Preview a reset; add `--apply` to archive results |
| `frames PROJECT` | Count production frames and compare with expected counts |
| `--help` | Show command help; each command also accepts `--help` |
| `--version` | Show the version |

| Command option | Default when omitted |
| --- | --- |
| `--replica N` | All replicas; available for `submit`, `status`, `reset`, and `frames` |
| `submit --dry-run` | Off: submit jobs; add this flag to print commands instead |
| `status --offline` | Off: query Slurm; add this flag to inspect local files only |
| `reset --from STAGE` | Required: no default starting stage |
| `reset --apply` | Off: preview only; add this flag to archive results |
| `reset --cancel` | Off: refuse reset if affected jobs are active; requires `--apply` to cancel |

`PROJECT` is required for `submit`, `status`, `reset`, and `frames`; it has no
default. `--help` and `--version` display information and exit when requested.

Frame checks never move files or count combined trajectories twice. Load your
AmberTools module before running `frames`. Missing or unexpected frame counts
return exit code `1`; command errors return `2`; success returns `0`.

## Generate options

Every optional setting has either a fixed default, an automatic value, or an
off/unset state. Required inputs have **no default**. Defaults are listed below;
the comparison with the original script follows the option tables.

| Option | Meaning / default |
| --- | --- |
| `--output DIR` | Required: a new campaign directory; no default |
| `--topology FILE` | Required: prepared AMBER topology; no default |
| `--coordinates FILE` | Required: starting coordinates; restart with velocities for production-only; no default |
| `--account NAME` | Defaults to `SLURM_ACCOUNT`; required if that variable is unset |
| `--workflow all\|prep\|production` | Stages to generate; default `all` |
| `--replicas N` | Independent replicas; default `1` |
| `--segments N` | Sequential production segments per replica; default `1` |
| `--ns N` | Nanoseconds per production segment; default `10` |
| `--combine` | Off by default; add a combination job to each replica |
| `--submit` | Off by default; submit immediately after generation |

### Simulation settings

| Option | Meaning / default |
| --- | --- |
| `--temperature K` | Target temperature; `310` K |
| `--ensemble npt\|nvt` | Constant pressure or volume; `npt` |
| `--pressure BAR` | NPT pressure; `1` bar |
| `--dt PS` | Production/equilibration timestep; `0.002` ps, maximum `0.002` |
| `--cutoff ANGSTROMS` | Nonbonded cutoff; `10` |
| `--restraint-mask MASK` | Preparation mask; `'!(:WAT,Na+,Cl-)'` |
| `--restraints FILE` | Optional production distance restraints; off by default |
| `--seed N` | Base random seed; generated and saved if omitted |
| `--ntwx N` | Production trajectory output interval; `25000` steps |
| `--ntpr N` | Production energy output interval; `50000` steps |
| `--ntwr N` | Production restart interval; `50000` steps |

### Cluster settings

| Option | Meaning / default |
| --- | --- |
| `--engine cuda\|cpu` | Production engine; `cuda`; preparation uses MPI CPU |
| `--gpu RESOURCE` | Slurm GPU request; `gpu:1` |
| `--cpus N` | MPI tasks for CPU production; `1` |
| `--prep-tasks N` | MPI tasks for each preparation job; `24` |
| `--mem SIZE` | Total production/analysis memory; `4G` |
| `--prep-mem SIZE` | Total preparation memory; `24G` |
| `--walltime TIME` | Production walltime; `0-03:00` |
| `--prep-walltime TIME` | Walltime for each min/heat/eq job; `0-03:00` |
| `--analysis-walltime TIME` | Combination walltime; `0-01:30` |
| `--partition NAME` | Unset: let Slurm choose its default partition |
| `--modules "LIST"` | Space-separated module list; see below |
| `--mail EMAIL` | Optional notifications; none by default |

Default modules:

```text
StdEnv/2023 gcc/12.3 openmpi/4.1.5 cuda/12.2 amber-pmemd/24.3 ambertools/25.0
```

Use `--modules ""` to inherit an already configured environment. Walltime accepts
`DD-HH:MM`, `HH:MM`, or `HH:MM:SS`. Verify modules, GPU requests, and resource
limits on Fir before your first submission.

### Which defaults came from the original?

The main production defaults were retained: **one segment, 10 ns per segment,
310 K, NPT, 1 bar, a 0.002 ps timestep, a 10 angstrom cutoff, `ntwx=25000`,
`ntpr=50000`, and a three-hour walltime**. CUDA production, 24 MPI tasks for
preparation, the preparation restraint mask, and the production module list
also came from the original. Preparation walltime remains three hours and
combination walltime remains 90 minutes. Combining trajectories and submitting
jobs remain opt-in.

These defaults or behaviors changed:

| Setting | Original | md-bud default |
| --- | --- | --- |
| Production checkpoints (`--ntwr`) | `5000000` steps: every 10 ns | `50000` steps: every 0.1 ns at the default timestep |
| Production memory (`--mem`) | `1G` | `4G` total |
| GPU request (`--gpu`) | `gpu:nvidia_h100_80gb_hbm3_3g.40gb:1` | `gpu:1`; no specific GPU model |
| Preparation memory (`--prep-mem`) | `1000MB` per MPI task | `24G` total; a different allocation request |
| Combination resources | 4 tasks, `2024M` per CPU | 1 task, `4G` total; uses `--mem` |
| Production distance restraints | Enabled by default (`nmropt=1`, `RST.dist`) | Off unless `--restraints FILE` is supplied |
| Random seeds | AMBER chooses seeds using `ig=-1` | An automatically generated base seed is saved; explicit seeds differ by replica and stage |
| Analysis modules | A separate older module list | The same configurable module list as preparation and production |
| Account and notifications | Personal values embedded in scripts | Your account or `SLURM_ACCOUNT`; email notifications off |
| Input/output names | Placeholder topology and prefix values | Explicit input paths and output directory required |
| Workflow | Stages selected individually with flags | `all`: minimization, heating, equilibration, and production |
| Dependencies | Continue after any predecessor termination | Continue only after predecessor success |

`--replicas` is new and defaults to one. Reset, status, frame checks, and their
control flags are also new. The generated workflow is not an exact reproduction
of the old script: preparation is split into three jobs, and production uses
restart files directly instead of imaging them between segments.

To use the original production checkpoint interval, memory, and GPU request:

```bash
./md-bud generate --output runs/original-resources \
  --topology system.prmtop --coordinates system.inpcrd --account your-allocation \
  --ntwr 5000000 --mem 1G --gpu gpu:nvidia_h100_80gb_hbm3_3g.40gb:1
```

This only restores those three settings; it does not enable production distance
restraints. Add `--restraints FILE` if your protocol needs them.

## Files and limitations

A campaign contains `inputs/`, `replica-001/`, further replica directories,
`md-bud.json`, and any reset archives. Inputs are copied under generic filenames
and checked against saved checksums before submission. Moving a campaign is fine;
changing its copied inputs requires a new campaign.

The built-in protocol is for periodic explicit solvent. It uses four minimization
steps, six 20 ps heating steps, and five restrained equilibration steps. It does
not prepare structures, choose force fields, or validate scientific suitability
for membranes, implicit solvent, or special ensembles. Starting from the same
production restart shares initial coordinates/velocities; thermostat seeds differ.
AMBER execution and resource settings still need a real Fir smoke test.

## Publishing and testing

MIT licensed; see [LICENSE](LICENSE). Source files contain no personal account or
email defaults. Your own generated campaigns can contain your account, optional
email, simulation data, and job metadata; keep them out of public repositories.

The release ZIP excludes `.git`, legacy scripts, and simulation data. Existing
Git history still identifies the original author. Use a fresh repository from the
clean release if you need to publish without that history; review the hosting
account and repository URL separately.

Run the offline tests with:

```bash
python3 -m unittest discover -s tests -v
```

References: [Slurm dependencies](https://slurm.schedmd.com/sbatch.html#OPT_dependency)
and [CPPTRAJ frame counting](https://amberhub.chpc.utah.edu/command-line-syntax/).
