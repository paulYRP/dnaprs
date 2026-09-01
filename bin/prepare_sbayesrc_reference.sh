#!/usr/bin/env bash
set -euo pipefail

ld_archive="$1"
annotation_archive="$2"
output_dir="$3"

[[ -r "$ld_archive" ]] || { echo "SBayesRC LD archive is unreadable: $ld_archive" >&2; exit 2; }
[[ -r "$annotation_archive" ]] || { echo "SBayesRC annotation archive is unreadable: $annotation_archive" >&2; exit 2; }

mkdir -p extracted_ld extracted_annotation "$output_dir/ld" "$output_dir/annotation"
unzip -q "$ld_archive" -d extracted_ld
unzip -q "$annotation_archive" -d extracted_annotation

first_block=$(find extracted_ld -type f -name 'block*.eigen.bin' -print -quit)
[[ -n "$first_block" ]] || { echo "The SBayesRC LD archive contains no block*.eigen.bin files." >&2; exit 3; }
ld_source=$(dirname "$first_block")
cp -LR "$ld_source"/. "$output_dir/ld/"

annotation_source=$(find extracted_annotation -type f -name 'annot_baseline*.txt' -print -quit)
[[ -n "$annotation_source" ]] || { echo "The annotation archive contains no annot_baseline*.txt file." >&2; exit 3; }
cp "$annotation_source" "$output_dir/annotation/annotation.txt"

block_count=$(find "$output_dir/ld" -type f -name 'block*.eigen.bin' | wc -l | tr -d ' ')
annotation_rows=$(awk 'END { print NR }' "$output_dir/annotation/annotation.txt")
printf 'reference_type\tfiles_or_rows\tstatus\n' > "${output_dir}.summary.tsv"
printf 'sbayesrc_ld\t%s\tPASS\nannotation\t%s\tPASS\n' "$block_count" "$annotation_rows" \
    >> "${output_dir}.summary.tsv"
