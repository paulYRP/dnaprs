# dnaprs outputs

Each run is isolated under `<output-dir>/<run_name>/`. When reporting is enabled, the
published scientific results are contained in one portable directory:

```text
<run_name>/
└── reports/
    ├── index.html
    ├── inputs.html
    ├── target-gwas.html
    ├── plink.html                # included when PLINK C+T was run
    ├── sbayesrc.html             # included when SBayesRC was run
    ├── prs.html
    ├── phenotype.html            # included when phenotype models were run
    ├── dictionary.html
    ├── logs.html
    ├── assets/
    ├── site_libs/
    ├── downloads/
    ├── figures/
    │   ├── tiff/
    │   ├── png/
    │   └── jpeg/
    └── provenance/
```

Open `reports/index.html`. All report links are relative, so the complete `reports/`
directory can be reviewed without an internet connection. Participant-level score and
phenotype files remain sensitive research data and require the same access controls as
their source data.

## Report pages

- **Overview** summarises the completed run and final score checks.
- **Inputs** contains validated manifests, resolved parameters, and input checksums.
- **Target and GWAS** reports target and GWAS preparation and variant retention.
- **PLINK** is included only when PLINK C+T was run.
- **SBayesRC** is included only when SBayesRC was run.
- **PRS** contains participant distributions, coverage, correlations, and method
  agreement when both methods are available.
- **Phenotype** is included only when phenotype models were declared.
- **Dictionary** defines report files, table columns, and figures.
- **Logs** contains method logs, software versions, the Nextflow task trace, timeline,
  execution report, and pipeline DAG.

Each table, figure, and log is downloadable in the section where it is shown. There is
no separate Downloads page.

## Input and provenance records

`downloads/run/preflight/` contains validated manifests and resolved scientific
settings. `downloads/qc/preflight/preflight_qc.tsv` records the inputs accepted before
large tasks. `downloads/run/preflight/input_checksums.tsv` identifies the source files
with SHA-256 checksums.

`provenance/output_files.tsv` is the file dictionary for the report bundle.
`provenance/data_dictionary.tsv` defines columns in tabular results.
`provenance/figure_manifest.tsv` records each figure, its dimensions, resolution,
description, and available formats.
`provenance/report_state.tsv` records which conditional report pages were included.
`provenance/report_software_versions.tsv` records the software used to build the
website.

## PLINK C+T results

PLINK results are organised under:

```text
downloads/
├── reference/plink_ct/
├── plink_ct/<trait_id>/
├── qc/plink_ct/<trait_id>/
└── logs/plink_ct/
```

The weight QC file reports the harmonised and clumped variant counts and retained
percentage. A low retained count is interpreted as score coverage and is not repaired
by choosing a threshold from the target phenotype.

## SBayesRC results

SBayesRC results are organised under:

```text
downloads/
├── sbayesrc/<trait_id>/summary/
├── sbayesrc/<trait_id>/model/
├── qc/sbayesrc/<trait_id>/
└── logs/sbayesrc/<trait_id>/
```

QC tables report LD-aligned variants, summary-imputed variants, model weights, non-zero
effects, and posterior inclusion probability summaries. The phenotype is not used to
filter or tune SBayesRC weights.

## Participant scores

Method-specific and combined scores are organised under `downloads/scores/`. The main
combined files are:

- `prs_scores_long.tsv`: one participant, trait, and method per row;
- `prs_scores_wide.tsv`: analysis-ready raw and standardised score columns;
- `score_qc.tsv`: participant counts, scoring coverage, and scale checks;
- `method_concordance.tsv`: Pearson and Spearman agreement when both methods exist.

Important long-form fields are:

| Field                      | Meaning                                                    |
| -------------------------- | ---------------------------------------------------------- |
| `raw_prs`                  | Weighted allele or dosage sum on the method's native scale |
| `used_variants`            | Variants contributing to the score                         |
| `cohort_mean`, `cohort_sd` | Values used for within-cohort standardisation              |
| `prs_z`                    | `(raw_prs - cohort_mean) / cohort_sd`                      |

Within-cohort standardisation gives a one-standard-deviation unit for models within the
cohort. It does not compare absolute genetic risk between cohorts.

## Phenotype associations

Phenotype results are under `downloads/phenotype/`:

- `phenotype_associations.tsv` contains the matching PRS coefficient, standard error,
  confidence interval, P value, and change in model fit;
- `phenotype_models_fitted.tsv` records the exact full and reduced formulas;
- `phenotype_plot_data.tsv` contains participant-level fitted values, residuals, and
  adjusted values used by the report.

For independent Gaussian outcomes, incremental fit is the increase in R-squared. For
generalised models it is the recorded pseudo-R-squared change. For mixed models it is
the AIC reduction. These are research association estimates, not clinical risk
predictions.

## Figures

Every quantitative plot is saved in three forms from the same R plot object:

- 300 dpi LZW-compressed TIFF for publication-quality output;
- 300 dpi lossless PNG for browser and presentation use;
- 300 dpi high-quality JPEG with a white background.

The PNG is displayed in the report. TIFF, PNG, and JPEG downloads are available beside
each figure.
