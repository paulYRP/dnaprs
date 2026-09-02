# Running dnaprs

## Inputs

dnaprs has one workflow and two levels of input. A beginner supplies directories; an
advanced user may replace either directory with a YAML list of explicit records. Both
routes are normalised to the same internal, provenance-bearing records. Users do not
create TSV manifests.

### Target genotypes

`--input` defaults to `data/plink/raw`. The directory must resolve unambiguously to one
coherent PGEN, BED, PED/MAP, BGEN, VCF/BCF, or GenomeStudio FinalReport dataset. PLINK
sets require all companion files; GenomeStudio requires one assay manifest. Discovery
never chooses alphabetically between several candidates. For multiple cohorts, use a
YAML `input` list with an `id`, `path`, `format`, and optional companions per cohort.

### GWAS summary statistics

`--gwas` defaults to `data/gwas/raw`. Each supported text or compressed-text GWAS is
inspected independently. Common column names are resolved only when their roles are
unambiguous. An unfamiliar release must be described in the advanced YAML `gwas` list,
including effect scale, study size, build, and source-column mapping.

The workflow validates numeric effects, uncertainty, P values, frequencies, alleles,
coordinates, duplicate IDs, palindromic ambiguity, MAF, optional INFO, and orientation.
It does not recognise a private cohort, trait, or filename.

## Beginner commands

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

## Advanced `params.yml`

Run [`examples/params.yml`](../examples/params.yml) with:

```bash
nextflow run . -profile singularity -params-file examples/params.yml -resume
```

The YAML may use simple paths or explicit record lists. A multi-model declaration is:

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
