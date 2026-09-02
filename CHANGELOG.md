# Changelog

All notable changes to dnaprs are recorded here.

## 1.0.0dev - 2026-09-01

- Added a strict Nextflow DSL2 workflow with typed parameters and nf-schema
  validation.
- Added directory discovery and explicit YAML records for raw genotype and GWAS
  inputs, with optional phenotype files and parameter overrides.
- Corrected namespaced staging and symlinked PLINK-prefix resolution for
  directory-based target, GWAS, and reference inputs, including paths with spaces,
  and added resolver regression tests.
- Added support for PLINK 1 BED, PLINK 2 PGEN, PED/MAP, BGEN, VCF, and Illumina
  GenomeStudio inputs without dataset-specific names.
- Added validated local, cached, and downloaded reference modes, including pinned
  Beagle, unbref3, dbSNP, genome, genetic-map, ancestry, and SBayesRC assets.
- Added raw-genotype exploration, target conversion, two-stage QC, chromosome-level
  Beagle imputation, chromosome gathering, ancestry projection, and integrated
  participant decisions.
- Added fixed PLINK 2 C+T scoring with an independent LD reference.
- Added the SBayesRC 0.2.6 `tidy`, `impute`, model, and score sequence.
- Added concurrent trait, chromosome, scoring-method, and phenotype-model tasks with
  portable process labels for executor-specific CPU and memory allocation.
- Added combined raw and within-cohort standardised score tables.
- Added optional Gaussian, binomial, and mixed-model phenotype associations
  driven directly by CSV columns and YAML model records.
- Added typed-versus-imputed PLINK sensitivity and scoring-variant coverage checks.
- Added a tabbed Quarto report containing genotype, GWAS, QC, imputation, ancestry,
  PRS, phenotype, provenance, and execution-log results, with downloadable TIFF
  figures and XLSX tables.
- Added standard, Docker, Apptainer, and artificial stub-test profiles.
- Separated process software into pinned analysis, PLINK, SBayesRC, and report
  environments, while keeping production references outside containers.
- Added an executor-independent PBS example that requests HPC CPU and memory once and
  lets Nextflow schedule the pipeline; Apptainer is loaded through Singularity for
  compatibility with the available module environment.
- Kept trace, timeline, report, and DAG replacement independent from scientific
  result replacement so resumed runs refresh their provenance files.
- Organised published results under only `data/`, `figures/`, `logs/`, and `reports/`,
  rooted at `dnaprs/<run_name>` by default.
- Reworked the README, workflow diagram, documentation, and operational test bundle;
  removed unused template and site-specific files.
