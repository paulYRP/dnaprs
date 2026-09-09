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

`--gwas` defaults to `data/gwas/raw`. The pipeline accepts `.tsv`, `.txt`, `.csv`, and
VCF-style summary statistics as plain text, gzip (`.gz`), or block-gzip (`.bgz`). Each
GWAS is inspected independently. Common column names are resolved only when their roles
are unambiguous. An unfamiliar release must be described in a YAML `gwas` list,
including effect scale, study size, build, and source-column mapping.

Common GWAS QC requires finite effects, uncertainty, P values, frequencies and sample
sizes; autosomal single-base alleles; MAF at least 0.01; and the declared INFO threshold.
It includes INFO when collapsing exact records and excludes conflicting canonical
variant keys. PLINK C+T then removes every palindromic variant and orients effects to
the final target ALT allele. A target-aligned candidate with `P=0` stops clumping;
review the source P-value rather than dropping or replacing it automatically.
SBayesRC receives the common clean GWAS without this PLINK-specific filter.

### SBayesRC scoring genotypes

For imputed targets, the pipeline prepares chromosome-specific SBayesRC PGEN files
from the completed VCFs. It imports `DS` dosages, keeps eligible participants and
preserves existing variant IDs, including rsIDs. Only missing IDs receive the
`chromosome:position:REF:ALT` form. The shared coordinate-based genotypes used by
PLINK C+T remain unchanged.

Preparation runs once per cohort and chromosome and is reused across traits. It
requires unique variant IDs and a complete chromosome 1-22 set for SBayesRC scoring.
The original posterior weights are applied with `1 2 3 header no-mean-imputation`.
Model variants absent from the target are reported; a chromosome with no usable
weight matches stops with an ID and allele-matching error.

No manual conversion or extra parameter is required. This preparation uses completed
imputed VCFs and does not repeat Beagle imputation. Keep the Nextflow cache and work
files when resuming a failed scoring run.

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
record. Analysis values may differ between visits. It reports a participant who lacks
the requested value. Mixed models may select several values and retain one agreeing
record per participant and value; the declared
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
ancestry_pcs: 6
ancestry_percentile: 0.99
target_imputation: true
imputation_dr2: 0.80
seed: 20260829
```

The 0.10 checkpoint feeds imputation; the 0.01 checkpoint supports direct-genotype
sensitivity work. MAF and HWE are reported but are not filters by default. Reference
ancestry uses exact typed-variant matches, excludes the extended MHC, and applies the
declared minor allele frequency, missingness, and LD-pruning rules. It estimates six
unrelated 1000 Genomes reference axes, projects the target, and applies the empirical
European distance percentile. Participant decisions
separate technical eligibility, relatedness, ancestry, score eligibility, and primary
analysis.

Raw marker matching requires the complete assay pair, or its strand complement, to
occur within one complete dbSNP SNP record. Two ALT alleles can therefore identify a
candidate. Records containing non-SNP alleles are excluded as a whole. Exactly one
compatible rsID is required; matching does not establish the final genomic REF allele.

Raw marker preparation counts stored calls before excluding all-missing probes. It
counts Y calls separately from sex-aware missingness. Duplicate rsID groups must share
one manifest-derived assay pair and have completely concordant overlapping calls.
The retained probe has the highest call rate, then an exact source-rsID match, then
the earliest source row. Groups without overlapping calls are excluded.

GRCh37 orientation preserves native genotypes, dosages and sample metadata for QC.
For array inputs, `bcftools +fixref` determines the TOP-to-forward transformation from
marker alleles and reference sequence; PLINK applies the allele recoding to PGEN.
Reference checks remain required before imputation.
An assay/reference conflict stops preparation without replacing an observed allele.
After imputation filtering, every retained participant dosage must be present, finite
and between zero and two. Missing or invalid retained dosages stop scoring.

Required post-QC calculations must complete for retained participants. Missing or
invalid heterozygosity, relatedness or ancestry results stop scoring. Participants
already removed by sample QC remain in the decision table without requiring those
downstream results. A sex check without X data or recorded sex is `NOT_APPLICABLE`,
not a completed calculation.

Existing PLINK PRS columns are preserved as `_NIMP` when no archive exists. If both
columns already exist, the output preserves `_NIMP` and replaces the current score.
The pipeline does not edit the phenotype input. It preserves non-score values and
all timepoints; model selection still uses the declared `timepoint_values`.

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

Target import, marker resolution and GRCh37 finalisation are cached separately.
Reference chromosomes run in up to six groups balanced by source-file size, then merge
once. PLINK alignment reuses one compatibility index per target/reference pair across
traits. The scientific filters and scoring methods remain unchanged.

Default resource limits are 32 CPUs, 144 GB memory and 72 hours per task. Target
preparation has a 30-hour allowance while full-cohort performance is measured.
Compatibility indexing and trait alignment retain 72 GB until the full GWAS trace
supports a lower request.
Institutional configuration may lower these limits to match the scheduler.

Keep site-specific module commands, queue/account settings, and scheduler launchers
outside the pipeline repository. A scheduler launcher should request a small Nextflow
driver job; the selected executor then submits the independent pipeline tasks with the
portable process resources defined here. There is no Makefile or site-specific
pipeline configuration.

## Resume

Use `-resume` with the same revision, parameters, inputs, launch directory, and work
directory. Both the Nextflow cache and work directory are required. Published result
copies are never consumed by downstream tasks.

After a pipeline update, resume from the existing launch and work directories. Nextflow
reruns tasks whose inputs, scripts, parameters or containers changed. Keep the cache
and work files until the resumed run completes; deleting published results is not
required. Reuse of prepared references and compatibility indexes uses this task cache.

## Report resources

`RENDER_REPORT` requests 100 GB initially, then 250 GB and 500 GB for resource-related
retries. Quarto's JavaScript heap allowance is half the allocated task memory, leaving
headroom for other allocations. Test and stub profiles retain small requests.

The selected scheduler queue must support these allocations. Larger requests may wait
longer in the queue. Preserve the work directory and `.nextflow/cache`, and resume the
specific failed session after a report error to reuse unchanged successful analysis tasks.
