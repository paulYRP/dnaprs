# Running dnaprs on Aquarius

## Execution design

`aqua.pbs` starts one small Nextflow controller job. The controller uses the `pbspro`
executor to submit separate compute jobs for each pipeline process. It must be submitted
from the `dnaprs` directory containing `main.nf` and `nextflow.config`.

The current profile uses these confirmed modules:

| Work       | Modules                                      |
| ---------- | -------------------------------------------- |
| Controller | `Java/21.0.8`                                |
| PLINK work | `GCC/13.2.0`, `PLINK/2.0.0-a.6.12`           |
| R work     | `GCC/13.2.0`, `R/4.4.1`                      |
| Report     | `GCC/13.2.0`, `R/4.4.1`, `quarto/1.7.32-x64` |

SBayesRC processes use the fixed R 4.4.1 library supplied through
`DNAPRS_R_LIBS_USER`. The pipeline does not require Docker or Apptainer on Aquarius.

## One-time setup

Install the pinned standalone Nextflow executable in a stable location and cache the
fixed `nf-schema` plugin while network access is available. The launcher defaults to:

```text
$HOME/.local/share/dnaprs/nextflow-25.10.4-dist
$HOME/.nextflow/plugins/nf-schema-2.7.2
```

Cache the required plugin before enabling offline execution:

```bash
$HOME/.local/share/dnaprs/nextflow-25.10.4-dist plugin install nf-schema@2.7.2
```

Confirm the controller can see PBS Pro:

```bash
module purge
module load Java/21.0.8
command -v qsub
command -v qstat
```

Build one R library with `R/4.4.1` and `GCC/13.2.0`. It must contain SBayesRC 0.2.6,
BH 1.84.0-0, data.table, ggplot2, scales, knitr, rmarkdown, glm2, lme4, and openssl.
Keep the package manifest with the run records.

## Directory choice

Copy or clone the complete `dnaprs` directory into the project. From the example project
layout, start here:

```bash
cd /mnt/hpccs01/home/ruizpint/QPASST/dnaprs
```

Manifest paths are resolved from this launch directory. Inputs kept in the parent project
can therefore be written as `../data/...`, or as absolute shared-storage paths.

The work directory must be on shared storage visible to every compute node. Node-local
`/tmp` must not be used for the Nextflow work directory because it prevents reliable
resume after the controller job ends.

## Production submission

```bash
cd /mnt/hpccs01/home/ruizpint/QPASST/dnaprs

export DNAPRS_PARAMS_FILE=runs/ukr/params.yaml
export DNAPRS_R_LIBS_USER=/mnt/hpccs01/home/ruizpint/QPASST/rLibrary/dnaprs_R4.4.1
export DNAPRS_WORK_DIR=/mnt/hpccs01/scratch/ruizpint/dnaprs-work
export DNAPRS_OUTPUT_DIR=/mnt/hpccs01/home/ruizpint/QPASST/results/dnaprs
export DNAPRS_OFFLINE=true

qsub aqua.pbs
```

The full SBayesRC model requests 64 CPUs and 512 GB by default and only one such model is
submitted at a time. Lighter tasks request smaller resources. The configured hard limit
is 64 CPUs and 6 TB; increasing a task above its tested requirement should occur only in
the Aquarius configuration, not inside scientific scripts.

## Logs and monitoring

The PBS launcher creates:

```text
logs/prs/dnaprs_<PBS_JOBID>.log
logs/prs/dnaprs_nextflow_<PBS_JOBID>.log
```

Nextflow writes its execution trace, report, timeline, and DAG under
`<output-dir>/<run_name>/reports/provenance/`. For a failed process, use the work
directory recorded in the trace and inspect
`.command.sh`, `.command.out`, `.command.err`, and `.command.log`.

Monitor the controller and process jobs with:

```bash
qstat -u "$USER"
tail -f logs/prs/dnaprs_<PBS_JOBID>.log
```

A start banner confirms only that the controller passed its software checks. Completion
is confirmed by `dnaprs completed successfully.`, a zero PBS exit status, and the final
result/provenance files.

## Resume after an interruption

Use the same pipeline directory, parameters, inputs, Nextflow work directory, and output
directory:

```bash
export DNAPRS_RESUME=true
qsub aqua.pbs
```

Nextflow will rerun only tasks whose inputs, command, code, or settings changed. Do not
edit completed files in the work directory and do not delete `.nextflow/cache` before a
resume.
