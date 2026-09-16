#!/usr/bin/env bash
set -euo pipefail

pipeline_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
test_dir=$(mktemp -d)
trap 'rm -rf -- "$test_dir"' EXIT
cd "$test_dir"

# Stop at the next validation check, without submitting analysis jobs.
accepted_message='--methods must contain one or more of: plink_ct,sbayesrc'
collision_message='Run output directory contains existing results:'
case_number=0

check_startup() {
    local case_name=$1
    local expected_message=$2
    shift 2
    case_number=$((case_number + 1))
    local exit_code=0
    nextflow -log "$test_dir/case-${case_number}.nextflow.log" run "$pipeline_dir/main.nf" \
        -profile test \
        -c "$pipeline_dir/tests/output_directory.config" \
        -work-dir "$test_dir/work" \
        --outdir "$test_dir/results" \
        --run_name "$case_name" \
        --methods output_guard_test \
        "$@" > "$test_dir/case-${case_number}.txt" 2>&1 || exit_code=$?

    if [[ $exit_code -ne 1 ]] || ! grep -Fq -- "$expected_message" "$test_dir/case-${case_number}.txt"; then
        cat "$test_dir/case-${case_number}.txt"
        printf 'FAIL: %s (exit %s)\n' "$case_name" "$exit_code" >&2
        exit 1
    fi
    printf 'PASS: %s\n' "$case_name"
}

check_startup absent "$accepted_message"

mkdir -p results/empty
check_startup empty "$accepted_message"

mkdir -p results/rounded/reports/provenance
printf 'task_id\thash\n' > results/rounded/reports/provenance/execution_trace.txt
touch -d "@$(date +%s)" results/rounded/reports/provenance/execution_trace.txt
check_startup rounded "$accepted_message"

mkdir -p results/bootstrap/reports/provenance
for filename in execution_trace.txt execution_report.html execution_timeline.html pipeline_dag.html; do
    printf 'Previous startup record\n' > "results/bootstrap/reports/provenance/$filename"
    touch -d '2000-01-01 UTC' "results/bootstrap/reports/provenance/$filename"
done
check_startup bootstrap "$accepted_message"

# Resume uses the same recorded session, not whichever run happened last.
resume_id=$(awk '/Session UUID:/ { print $NF; exit }' "case-${case_number}.nextflow.log")
[[ -n $resume_id ]]
mkdir -p results/bootstrap/data/scores
printf 'Existing scores\n' > results/bootstrap/data/scores/prs_scores_long.tsv
check_startup bootstrap "$collision_message"
grep -Fxq 'Existing scores' results/bootstrap/data/scores/prs_scores_long.tsv
check_startup bootstrap "$accepted_message" --overwrite true
check_startup bootstrap "$accepted_message" -resume "$resume_id"
grep -Fxq 'Existing scores' results/bootstrap/data/scores/prs_scores_long.tsv

# Do not exempt the whole provenance directory or files with matching basenames.
mkdir -p results/other_provenance/reports/provenance
printf 'Existing output\n' > results/other_provenance/reports/provenance/output_files.tsv
check_startup other_provenance "$collision_message"

mkdir -p results/same_basename/data
printf 'Existing output\n' > results/same_basename/data/execution_trace.txt
check_startup same_basename "$collision_message"

printf 'Existing file\n' > results/not_directory
check_startup not_directory 'Run output path exists but is not a directory:'
grep -Fxq 'Existing file' results/not_directory

printf 'Output-directory checks passed (%s cases).\n' "$case_number"
