# Running dnaprs

## Inputs

dnaprs accepts conventional raw-data directories or YAML lists of explicit records.
Both forms are normalised to the same validated, provenance-bearing TSV records. The
pipeline writes `targets.tsv`, `gwas.tsv`, `references.tsv`, and `models.tsv` for
internal hand-off and provenance; these generated files are not launch parameters.

### Target genotypes

`--input` defaults to `data/plink/raw`. The directory must resolve unambiguously to one
coherent PGEN, BED, PED/MAP, BGEN, VCF/BCF, or GenomeStudio FinalReport dataset. PLINK
sets require all companion files; GenomeStudio requires one assay manifest paired with
the FinalReport. This vendor assay manifest is distinct from the generated TSV records.
Discovery never chooses alphabetically between several candidates. For multiple cohorts, use a
YAML `input` list with an `id`, `path`, `format`, and optional companions per cohort.

### GWAS summary statistics

`--gwas` defaults to `data/gwas/raw`. Each supported text or compressed-text GWAS is
inspected independently. Common column names are resolved only when their roles are
unambiguous. An unfamiliar release must be described in a YAML `gwas` list,
including effect scale, study size, build, and source-column mapping.

The workflow validates numeric effects, uncertainty, P values, frequencies, alleles,
coordinates, duplicate IDs, palindromic ambiguity, MAF, optional INFO, and orientation.
It does not recognise a private cohort, trait, or filename.

## Minimal commands

With the conventional folders and automatic references:

```bash
nextflow run . -profile singularity -resume
```

Explicit equivalent:

```bash
nextflow run . \
  -profile singularity \
  --input data/plink/raw \
  --gwas data/gwas/raw \
  --outdir dnaprs \
  --run_name model1 \
  -resume
```

Optional phenotype association:

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

The four phenotype parameters are atomic: omit all four or provide all four.
`--covariates none` requests an unadjusted model. Use `--participant_id` when exactly
one CSV column cannot be matched unambiguously to target sample IDs. Binomial models
accept validated 0/1 data or explicit `--control_value` and `--case_value`; mixed
models also require `--group_column`.

For repeated phenotype records, declare the visit column and retained values together:

```yaml
participant_id: Sample_Name
timepoint_column: Timepoint
timepoint_values: [1]
```

A fixed model requires one value. The pipeline checks agreeing outcome and covariate
values within each participant and requested timepoint, then selects the first source
record. It reports a participant who lacks the requested value. Mixed models may select
several values and retain one agreeing record per participant and value; the declared
covariates still determine whether time is included in the model.

## Configured `params.yml`

Run [`examples/params.yml`](../examples/params.yml) with:

```bash
nextflow run . -profile singularity -params-file examples/params.yml -resume
```

The YAML may use simple paths or explicit record lists. A multi-model declaration is:

Relative input paths are resolved from the directory where `nextflow run` is launched.
Use absolute paths when a file is outside that directory. Files distributed with the
pipeline are resolved separately from the pipeline project directory.

```yaml
phenotype: data/pheno/pheno.csv
participant_id: participant_id
models:
  - id: depression
    outcome: depression_score
    covariates: [age, sex]
    model_type: gaussian
    score_ids: all
    primary: true
  - id: case_control
    outcome: diagnosis
    covariates: [age, sex, batch]
    model_type: binomial
    control_value: control
    case_value: case
    score_ids: [trait1]
```

Use `model_id: depression` to run one declared model. Each selected model-score pair is
an independent task; a one-model run follows the same code path.

## Reference modes

- `auto` (default): use valid cached/local assets and acquire missing pinned assets;
- `local`: require every selected role under `--references` and perform no scientific
  download;
- `download`: acquire the pinned bundle even when another local root is supplied.

A local reference root must contain one BREF3 panel file and one supported genetic map
for each autosome from 1 to 22. Panel names start with `chr<chromosome>` and end in
`.bref3`. Map names use `plink.chr<chromosome>.GRCh37.map` or
`chr<chromosome>.map`. Additional non-autosomal maps are accepted. Missing or duplicate
autosomal files stop validation before imputation.

The dbSNP directory must contain one VCF source and one assembly report. The VCF source
may use a standard VCF extension or a RefSeq name such as `GCF_000001405.25.gz`. When no
index is available, the pipeline creates an indexed copy in the task work directory and
does not modify the source directory.

The persistent cache defaults to `references/dnaprs/grch37-v1/`. It includes dbSNP157,
the GRCh37 FASTA/index, Beagle maps and chromosome panel, 1000 Genomes population and
related-sample metadata, pinned Beagle/unbref3 JARs, and selected SBayesRC resources.
Each asset is checked independently, so a later failure does not restart successful
chromosome downloads.

For a cache-only run, create `references.yml`:

```yaml
reference_only: true
reference_mode: download
reference_dir: references/dnaprs
reference_bundle: grch37-v1
methods: plink_ct,sbayesrc
target_imputation: true
outdir: dnaprs
run_name: references
```

Run `nextflow run . -profile singularity -params-file references.yml -resume`.

## QC and analysis controls

```yaml
genome: GRCh37
methods: plink_ct,sbayesrc
sample_missingness: 0.02
imputation_variant_missingness: 0.10
direct_variant_missingness: 0.01
maf_filter: 0
hwe_filter: 0
ancestry_pcs: 10
ancestry_percentile: 0.99
target_imputation: true
imputation_dr2: 0.80
seed: 20260829
```

The 0.10 checkpoint feeds imputation; the 0.01 checkpoint supports direct-genotype
sensitivity work. MAF and HWE are reported but are not filters by default. Reference
ancestry uses exact typed-variant matches, unrelated 1000 Genomes axes, target
projection, and the empirical European distance percentile. Participant decisions
separate technical eligibility, relatedness, ancestry, score eligibility, and primary
analysis.

Each available autosome is emitted as one independent Beagle task. Nextflow schedules
these tasks concurrently when executor capacity is available and merges them only after
all expected chromosome checks pass. PLINK C+T is also calculated from the stricter
direct-genotype checkpoint; variant coverage and standardised participant scores are
compared with the primary imputed-target score as a sensitivity analysis.

## Portable HPC execution

The pipeline contains task CPU, memory, time, threading, and bounded resource-retry
rules. Do not put scheduler or site details in `params.yml`. Select the institutional
Nextflow executor in the installed environment or an external Nextflow configuration;
Nextflow submits independent cohorts, traits, methods, references, and model-score pairs
as capacity allows.

Keep site-specific module commands, queue/account settings, and scheduler launchers
outside the pipeline repository. A scheduler launcher should request a small Nextflow
driver job; the selected executor then submits the independent pipeline tasks with the
portable process resources defined here. There is no Makefile or site-specific
pipeline configuration.

## Resume

Use `-resume` with the same revision, parameters, inputs, launch directory, and work
directory. Both the Nextflow cache and work directory are required. Published result
copies are never consumed by downstream tasks.
