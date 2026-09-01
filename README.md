<p align="center">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="docs/images/nf-core-dnaprs_logo_dark.png">
    <img src="docs/images/nf-core-dnaprs_logo_light.png" alt="dnaprs">
  </picture>
</p>

[![nf-test status](https://github.com/paulYRP/dnaprs/actions/workflows/nf-test.yml/badge.svg)](https://github.com/paulYRP/dnaprs/actions/workflows/nf-test.yml)
[![nf-core linting status](https://github.com/paulYRP/dnaprs/actions/workflows/linting.yml/badge.svg)](https://github.com/paulYRP/dnaprs/actions/workflows/linting.yml)
[![Nextflow](https://img.shields.io/badge/version-%E2%89%A525.10.4-green?style=flat&logo=nextflow&logoColor=white&color=%230DC09D&link=https%3A%2F%2Fnextflow.io)](https://www.nextflow.io/)
[![nf-core template version](https://img.shields.io/badge/nf--core_template-4.1.0-green?style=flat&logo=nfcore&logoColor=white&color=%2324B064&link=https%3A%2F%2Fnf-co.re)](https://github.com/nf-core/tools/releases/tag/4.1.0)
[![nf-test](https://img.shields.io/badge/tested_with-nf--test-337ab7)](https://www.nf-test.com/)
[![MIT licence](https://img.shields.io/badge/license-MIT-111111)](LICENSE)

## Introduction

**dnaprs** is a portable Nextflow DSL2 pipeline for generating and evaluating polygenic
risk scores from raw target genotypes and raw GWAS summary statistics. It is derived
from the nf-core template and follows nf-core practices, but it is maintained at
`paulYRP/dnaprs` and is not an official `nf-core/*` pipeline.

The pipeline discovers PGEN, BED, PED/MAP, BGEN, VCF/BCF, or GenomeStudio FinalReport
inputs; records raw genotype EDA; resolves markers against pinned GRCh37 references;
applies two target-QC checkpoints; and projects participants onto unrelated 1000 Genomes
reference PCA axes. It can impute with Beagle, generate PLINK C+T and/or SBayesRC
scores, apply participant eligibility consistently, and fit one or many declared
phenotype models. Source files are read-only.

![dnaprs workflow](docs/images/dnaprs-workflow.svg)

The workflow publishes one run under `dnaprs/model1/` by default:

- `data/` — generated records, checkpoints, QC, weights, scores, and model tables;
- `figures/` — a convenient copy of every report figure;
- `logs/` — scientific and execution logs;
- `reports/` — the portable HTML website, downloads, and provenance.

## Beginner run

Install Java 17 or newer, Nextflow 25.10.4 or newer, and Docker, Apptainer, or
Singularity. With one raw target dataset in `data/plink/raw/` and raw GWAS files in
`data/gwas/raw/`:

```bash
nextflow run . -profile singularity -resume
```

The defaults are equivalent to:

```bash
nextflow run . \
  -profile singularity \
  --input data/plink/raw \
  --gwas data/gwas/raw \
  --outdir dnaprs \
  --run_name model1 \
  -resume
```

Add one phenotype model directly from the user-provided CSV:

```bash
nextflow run . \
  -profile singularity \
  --input data/plink/raw \
  --gwas data/gwas/raw \
  --phenotype data/pheno/pheno.csv \
  --outcome depression_score \
  --covariates age,sex \
  --model_type gaussian \
  -resume
```

`depression_score`, `age`, and `sex` are examples, not built-in names. If the
participant-ID column cannot be matched uniquely to target IDs, also supply
`--participant_id <column>`. Without phenotype parameters, PRS generation and the
genetic report still complete.

## References

`--reference_mode auto` is the default. The pipeline downloads only the required assets
from a pinned catalogue, verifies size and checksum (or the documented pinned ETag when
the source provides no digest), and reuses each valid asset from
`references/dnaprs/grch37-v1/`. Beagle and unbref3 are handled in the same way.

To use an existing reference collection:

```bash
nextflow run . \
  -profile singularity \
  --references /path/to/rData \
  --reference_mode local \
  -resume
```

To build or verify only the reusable reference cache, set `reference_only: true` in a
small YAML parameter file and run with `-params-file`.

## Advanced and HPC runs

Advanced users declare paths, GWAS column roles, thresholds, selected methods, and
multiple phenotype models in YAML; they never create target/GWAS/reference/model TSV
manifests. See [the complete usage guide](docs/usage.md) and
[`examples/params.yml`](examples/params.yml).

Processes use standard nf-core-style resource labels and pass `task.cpus` and bounded
task memory to capable tools. Nextflow handles cohort, chromosome-imputation, trait,
method, reference, and model tasks concurrently, then validates deterministic gathers.
Scheduler, queue, project, and filesystem settings remain
outside the pipeline, so the same workflow runs locally or through PBS Pro, Slurm, SGE,
LSF, and other Nextflow executors. The workspace-level `test/` directory contains the
minimal QUT PBS/Singularity acceptance launcher.

SBayesRC is a high-memory, long-running method. Use `--methods plink_ct` for a smaller
run; selecting fewer methods changes the scientific work requested, not the executor.

## Report and validation

Open `dnaprs/model1/reports/index.html` after completion. The report preserves every
available plot and provides SVG plus high-resolution PNG, TIFF, and JPEG downloads. Its
Logs page displays execution artifacts in expandable, scrollable panels.
When imputation is enabled, the PLINK page also reports typed-versus-imputed scoring
coverage and agreement between the primary imputed score and direct-genotype sensitivity
score; sensitivity scores are not added to phenotype models.

Run the pipeline-level stub regression test with:

```bash
nf-test test tests/default.nf.test --ci
```

Run the R and report smoke suite from PowerShell with:

```powershell
.\tests\smoke.ps1
```

These are software-contract tests, not biological or clinical validation. Phenotype
associations are research estimates and do not make a PRS a clinical risk prediction.

## Documentation

- [Input, reference, phenotype, and run options](docs/usage.md)
- [Output files and interpretation](docs/output.md)
- [Pinned software environments](containers/README.md)
- [Contributing and testing](docs/CONTRIBUTING.md)
- [Citations](CITATIONS.md)

The code is released under the [MIT licence](LICENSE). External tools, GWAS data, and
reference resources retain their own licences and access conditions.
