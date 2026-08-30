# Running dnaprs

## Input principle

The target genotype and GWAS must have completed their study-specific QC before they
enter dnaprs. The pipeline performs only the compatibility checks needed to construct
scores: file contracts, participant identifiers, genome build, allele columns, finite
GWAS values, required references, score coverage, and participant-level score output.

All paths in a manifest may be absolute or relative to the directory from which
Nextflow is launched. Inputs are read but are never edited in place.

## Parameter file

Use a YAML file for scientific runs so the settings can be reviewed and archived:

```yaml
run_name: ukr_prs
input: manifests/target.tsv
gwas_manifest: manifests/gwas.tsv
reference_manifest: manifests/reference.tsv
manifest_base: .
phenotype_file: data/pheno/pheno.csv
phenotype_models: manifests/phenotype_models.tsv
methods: plink_ct,sbayesrc
genome_build: GRCh37
seed: 20260829
report_enabled: true
```

`run_name` must be unique for a new analysis. Results and technical reports are written
under `<output-dir>/<run_name>/`.

## Target manifest

One row describes each QC-completed scoring cohort.

| Column     | Meaning                                                             |
| ---------- | ------------------------------------------------------------------- |
| `cohort`   | Unique short cohort name                                            |
| `role`     | `target` or `validation`                                            |
| `format`   | `pgen`, `bed`, `vcf`, or `bgen`                                     |
| `genotype` | File, PLINK prefix, or chromosome pattern such as `cohort_chr{chr}` |
| `sample`   | BGEN sample file when required; otherwise blank                     |
| `keep`     | Optional two-column PLINK participant keep file                     |
| `build`    | `GRCh37` or `GRCh38`                                                |
| `ancestry` | Recorded ancestry description                                       |
| `dosage`   | VCF dosage field, normally `DS`                                     |

PGEN inputs require `.pgen`, `.pvar`, and `.psam`; BED inputs require `.bed`, `.bim`,
and `.fam`. SBayesRC scoring expects autosomal chromosome files 1–22 after preparation.

Variant identifiers must be compatible with the selected LD reference and GWAS. dnaprs
preserves existing target identifiers and creates chromosome-position-allele identifiers
only when an input identifier is missing. A cohort whose identifiers were not aligned
during target QC must be aligned before scientific scoring.

## GWAS manifest

One row describes each QC-completed GWAS. dnaprs does not guess source column names.

| Column                                  | Meaning                                                       |
| --------------------------------------- | ------------------------------------------------------------- |
| `trait_id`                              | Unique trait identifier, such as `MDD`                        |
| `prs_name`                              | Stable output name, such as `MDD_PRS`                         |
| `path`                                  | Tab-delimited GWAS file, optionally gzip-compressed           |
| `build`                                 | Genome build matching the target and references               |
| `ancestry`                              | GWAS ancestry description                                     |
| `effect_type`                           | `beta`, `log_or`, or `or`                                     |
| `sample_size`                           | Reported or effective study sample size                       |
| `snp_col`, `chr_col`, `bp_col`          | Source variant identifier and position columns                |
| `effect_allele_col`, `other_allele_col` | Source allele columns                                         |
| `beta_col`, `se_col`, `p_col`           | Source effect, standard error, and P-value columns            |
| `freq_col`                              | Effect-allele-frequency column; required for SBayesRC         |
| `n_col`                                 | Per-variant sample-size column, or blank to use `sample_size` |

`log_or` means that the source effect is already the natural logarithm of the odds ratio.
Use `or` only when the declared source column contains positive odds ratios; dnaprs then
applies the natural logarithm before scoring.

For PLINK C+T, the declared GWAS SNP identifiers must match the independent PLINK LD
reference and the QC-completed target. For SBayesRC, `tidy()` checks the GWAS against its
LD panel, and the resulting scoring identifiers must exist in the target. Variant
coverage files show the retained overlap; the phenotype is never used to repair or tune
that overlap.

## Reference manifest

The selected methods determine the required rows.

| `reference_type` | Used by   | Required content                         |
| ---------------- | --------- | ---------------------------------------- |
| `plink_ld`       | PLINK C+T | Independent PLINK 2 PGEN LD reference    |
| `sbayesrc_ld`    | SBayesRC  | SBayesRC LD directory matching the build |
| `annotation`     | SBayesRC  | Matching functional annotation file      |

Every row also records a unique `reference_id`, build, ancestry, version, and optional
published checksum. Reference data should be installed once on shared storage and
referenced by path rather than copied into every run.

## Phenotype model manifest

Phenotype analysis is optional and begins only after score construction is fixed. Each
row declares one matching phenotype–PRS model.

| Column               | Meaning                                                                   |
| -------------------- | ------------------------------------------------------------------------- |
| `model_id`           | Unique model name                                                         |
| `outcome`            | Outcome column in the phenotype file                                      |
| `prs_name`           | Matching name from the GWAS manifest                                      |
| `family`             | `gaussian`, `binomial`, `poisson`, or another supported R family          |
| `covariates`         | Comma-separated covariate columns                                         |
| `participant_id`     | Participant identifier column                                             |
| `group_id`           | Group or participant column for a mixed model; blank for independent rows |
| `expected_direction` | Prespecified direction for interpretation                                 |
| `primary`            | `TRUE` or `FALSE`                                                         |

Independent Gaussian outcomes use `lm`. Other independent outcomes use `glm2`. A
declared `group_id` uses `lme4` with a random intercept. dnaprs compares the model with
the PRS against the same model without the PRS and reports the PRS coefficient, 95%
confidence interval, P value, and change in fit. These estimates describe association;
they do not validate a score as a clinical predictor.

## Offline terminal assistant

The assistant scans only a chosen workspace and uses numbered terminal selections. It
does not upload files or infer GWAS column meanings.

```bash
bin/dnaprs list /path/to/workspace
bin/dnaprs configure /path/to/workspace
```

Review the created YAML and TSV files. Then check strict syntax and configuration:

```bash
bin/dnaprs check \
  --params-file /path/to/workspace/runs/my_run/params.yaml \
  --profile docker
```

The `check` command is static. The full input preflight runs as the first process of an
actual workflow.

## Local run

Build the fixed image first, then start the workflow:

```bash
docker build --pull -t dnaprs:1.0.0 -f containers/Dockerfile .

nextflow run . \
  -profile docker \
  -params-file runs/my_run/params.yaml \
  --outdir results
```

The `standard` profile may be used when the exact required tools are already on `PATH`.
The `apptainer` profile is available for compatible Linux/HPC systems.

## Resume

Resume only with the same parameter file, inputs, work directory, launch directory, and
pipeline version:

```bash
nextflow run . \
  -profile docker \
  -params-file runs/my_run/params.yaml \
  --outdir results \
  -resume
```

Nextflow reuses a task only when its cached command and inputs still match. Do not remove
the work directory or `.nextflow/cache` while a run may need to resume.

## Test profiles

The committed fixtures are synthetic software tests. Test the R functions and report
from PowerShell:

```powershell
.\tests\smoke.ps1 -OutputDirectory ..\test\synthetic_analysis
```

Use strict stub mode for the complete two-method Nextflow graph:

```bash
nextflow run . -profile test_full -stub-run \
  --outdir ../test/nextflow_results \
  -work-dir ../test/nextflow_work
```

`test` executes the PLINK branch when Linux PLINK and R are available. `test_full`
checks PLINK and SBayesRC process connections. Synthetic and stub tests do not validate
biological performance.
