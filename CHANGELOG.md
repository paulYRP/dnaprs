# Changelog

All notable changes to dnaprs are recorded here.

## 1.0.0dev - 2026-09-01

- Corrected participant eligibility totals and added a separate reason breakdown.
- Split QQ and Manhattan figures by GWAS with lightweight previews and independent high-resolution downloads.
- Included table chunks created during rendering in the published report.
- Replaced large inline report tables with on-demand chunks, bounded browser caching, cancellable filters and complete-file downloads.
- Set report memory requests to 100, 250 and 500 GB across three attempts, with a separate Quarto heap allowance and small test-profile requests.
- Added small report-table regression tests and made large GWAS rendering tests opt-in.

- Reduced report memory use with deterministic GWAS display selection, per-trait reads and exact histogram counts, without changing PRS calculations or CSV outputs.
- Added full-data and displayed-point counts to report provenance and integer64 support to report container 1.0.2.
- Reduced image-export memory use by reusing prepared plots and writing vector SVG files with svglite.

- Fixed HTML report preparation by publishing SBayesRC genotype directories separately from report downloads, while preserving `phenoPRS.csv` and existing output paths.
- Added report-input checks for directories, missing files, broken links and unreadable files before copying or calculating checksums.
- Fixed SBayesRC sample selection when original family IDs differ from VCF sample names, preserving original IDs and sample order for scoring.
- Fixed SBayesRC genotype log outputs for compatibility with Nextflow 25.10.4.
