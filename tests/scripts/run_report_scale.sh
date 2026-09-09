#!/usr/bin/env bash
set -euo pipefail
repo_root=$(pwd)
record_dir=${1:?Provide the test-record directory}
mkdir -p "$record_dir"
record_dir=$(cd "$record_dir" && pwd)
test_parent=$(mktemp -d /tmp/dnaprs-report-scale.XXXXXX)
test_dir="$test_parent/report"
Rscript tests/scripts/test_report_scale.R "$test_dir" "${2:-6500000,7000000,7500000}"
cd "$test_dir"
export DNAPRS_REPORT_INPUTS="$test_dir/inputs"
export DNAPRS_OUTPUT_MANIFEST="$test_dir/output_files.tsv"
export QUARTO_VERSION
QUARTO_VERSION=$(quarto --version)
/usr/bin/time -v -o "$record_dir/prepare.resources.txt" \
    Rscript prepare-report.R > "$record_dir/prepare.log" 2>&1
/usr/bin/time -v -o "$record_dir/render.resources.txt" \
    quarto render . > "$record_dir/render.log" 2>&1
Rscript "$repo_root/tests/scripts/check_report_render.R"
cp provenance/gwas_plot_selection.tsv provenance/figure_manifest.tsv "$record_dir/"
cp figures/preview/gwas_qq_*.png figures/preview/gwas_manhattan_*.png "$record_dir/"
