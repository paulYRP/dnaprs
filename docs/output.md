# dnaprs outputs

Each run is published below `<outdir>/<run_name>/`; the defaults are
`dnaprs/model1/`. A complete run has exactly four top-level result folders:

```text
model1/
|-- data/
|-- figures/
|-- logs/
`-- reports/
```

## `data/`

`data/` contains scientific tables and reusable checkpoints, organised by stage:

- `inputs/`: internally resolved targets, GWAS, references, models, effective settings,
  checksums, input checks, and the workflow-output index;
- `genotype_eda/<cohort>/`: untouched-input composition, missingness, allele frequency,
  heterozygosity, sex checks, relatedness, and descriptive internal PCA;
- `target_prep/<cohort>/`: target conversion, marker decisions, and checks;
- `target_qc/<cohort>/`: sample/variant decisions, reference ancestry projections,
  distance summary, and integrated `participant_decisions.tsv`;
- `target/prepared/<cohort>/`: separate imputation-ready and direct-genotype PGEN
  checkpoints;
- `target_imputation/<cohort>/`: deterministically gathered chromosome handoff, DR2
  distribution, sample/order, chromosome/key/dosage checks, and SHA-256 receipts;
- `gwas/<trait>/` and `qc/gwas/<trait>/`: harmonised statistics and decisions;
- `reference/`: prepared PLINK and SBayesRC resources used by the run;
- `plink_ct/`, `sbayesrc/`, `scores/`, and `qc/`: weights, method results, combined
  scores, coverage, and agreement;
- `phenotype/`: model declarations, coefficients, fit, permutation, influence, plotting,
  and analysis-ready phenotype/PRS tables;
- `pipeline_info/software_versions.yml`: combined tool versions.

Large public source references remain in `--reference_dir`; the run records their
validated receipt rather than publishing a second source copy.

`participant_decisions.tsv` contains `technical_pass`, `score_eligible`,
`related_flag`, `ancestry_flag`, `ancestry_distance`, `primary_analysis`, and a reason.
Scoring uses score-eligible participants; phenotype models use the integrated primary
analysis decision.

## `figures/`

This is a convenient copy of all report figures. The website retains its own copy under
`reports/figures/` so the complete `reports/` directory remains portable. Every plot is
generated from a published table and, when enabled, is available as vector SVG and
360 dpi PNG, TIFF, and JPEG. No plot is removed when a viewer has several figures.

## `logs/`

Logs are grouped by target, reference, imputation, GWAS, PRS method, and report stage.
The report Logs page displays supported text/HTML artifacts inline in expandable,
collapsible, scrollable panels. Execution trace, timeline, report, and DAG files are
also copied into report provenance when enabled.

## `reports/`

Open `reports/index.html`. The site may include Overview, Genotype EDA, Target PREP,
Target QC, Target Imputation, GWAS QC, PLINK PRS, SBayesRC PRS, Phenotype, and Logs,
depending on completed stages. Only Overview has introductory text. Figure viewers use
Previous, one selector, Next, and a counter; they do not repeat figure names as a second
button row.

Important portable subfolders are:

- `downloads/`: tables, logs, score files, and `dnaprs_report_tables.xlsx`;
- `figures/`: all report figure formats;
- `provenance/`: effective parameters, output/figure indexes, data dictionary, report
  state, software versions, and available Nextflow execution records;
- `assets/` and `site_libs/`: local styling and libraries required offline.

All links are relative. Treat participant-level score and phenotype files as sensitive
research data.

## Score interpretation

`prs_scores_long.tsv` contains one participant, trait, and method per row, including
raw PRS, cohort standardisation values, `prs_z`, variants used, and participant-decision
fields. `prs_scores_wide.tsv` contains analysis-ready score columns plus the same
eligibility fields. `method_concordance.tsv` compares methods when both were selected.
With imputation enabled, `scores/<cohort>/<trait>/sensitivity/` contains the typed-only
PLINK score and participant comparison; the matching QC table records shared,
imputed-only, and direct-only scoring variants plus Pearson and Spearman agreement. The
typed-only score is not included in the combined primary scores or phenotype models.

Within-cohort standardisation expresses effects per cohort standard deviation; it does
not make absolute genetic risk comparable between cohorts. Phenotype coefficients are
research association estimates, not clinical predictions.
