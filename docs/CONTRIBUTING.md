# Contributing to dnaprs

## Scope

Version 1 contains PLINK C+T and SBayesRC only. A change should have a clear scientific,
compatibility, or reproducibility purpose. New PRS methods require their own documented
input contract, fixed method settings, tests, outputs, and citations rather than being
inserted into an existing process.

Never commit participant genotypes, phenotypes, identifiers, restricted GWAS files, or
private reference resources. Tests must use artificial or openly licensed data.

## Code structure

- Keep one focused process in each `modules/local/<name>/main.nf` file.
- Keep the authoritative workflow connection in `workflows/dnaprs.nf`; do not duplicate
  the same process or workflow in another file.
- Use stable typed parameters in `main.nf` and keep `nextflow_schema.json` consistent.
- Use the default Nextflow v2 strict parser and avoid deprecated syntax.
- Stage scripts and templates as explicit process inputs.
- Write outputs to the task work directory and publish them through the top-level
  workflow output block. Do not add `publishDir` to processes.
- Use metadata-derived trait, cohort, and method folders so independent outputs cannot
  overwrite each other.
- Do not modify source files in place.

## Scientific rules

- Treat declared source inputs as immutable; record every conversion, correction,
  exclusion, and prepared-reference derivative.
- Keep target phenotypes out of GWAS harmonisation, variant selection, C+T thresholds,
  and SBayesRC model fitting.
- Preserve raw scores and record the exact standardisation population and formula.
- Keep phenotype models prespecified and separate from score construction.
- Record genome build, ancestry, effect scale, reference version, software version,
  variant coverage, and participant count.

## Required checks

Run before proposing a change:

```bash
find bin -maxdepth 1 -name '*.sh' -exec bash -n '{}' \;
nextflow-25.10.4 lint .
python3 -m nf_core pipelines lint --dir .
nf-test test tests/default.nf.test --ci
nextflow-25.10.4 run . -profile test_full -stub-run \
  --outdir /tmp/dnaprs-test-full \
  -work-dir /tmp/dnaprs-test-work
docker build --file containers/analysis/Dockerfile \
  --tag ghcr.io/paulyrp/dnaprs-analysis:1.0.0 .
docker build --file containers/plink2/Dockerfile \
  --tag ghcr.io/paulyrp/dnaprs-plink2:2.0.0-a.6.12 .
docker build --build-arg ANALYSIS_IMAGE=ghcr.io/paulyrp/dnaprs-analysis:1.0.0 \
  --file containers/report/Dockerfile \
  --tag ghcr.io/paulyrp/dnaprs-report:1.0.0 .
nextflow-25.10.4 run . -profile docker,test_full -stub-run \
  --outdir /tmp/dnaprs-docker-test \
  -work-dir /tmp/dnaprs-docker-work
```

On Windows, parse every R script, exercise the phenotype models, create every figure
format, and render the complete report with the artificial fixtures:

```powershell
.\tests\smoke.ps1 -OutputDirectory ..\test\synthetic_analysis
```

A change to PLINK or SBayesRC commands also requires a small real method test and
comparison with the previously accepted outputs using identical inputs.

Check that process and workflow names are unique:

```bash
grep -RHE '^(process|workflow)[[:space:]]+' --include='*.nf' .
```

Review the resolved local and Apptainer configurations before release. Every container
must be built, inspected for tool versions, and published with a stable version. Site
executors, queues, accounts, and storage paths remain in external institutional
configuration files.
