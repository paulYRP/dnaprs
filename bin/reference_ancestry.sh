#!/usr/bin/env bash
set -euo pipefail

cohort="$1"
target_dir="$2"
reference_dir="$3"
pcs="$4"
percentile="$5"
threads="$6"
classifier="$7"

target_prefix="${target_dir}/${cohort}"
reference_prefix="${reference_dir}/all_reference"
population="${reference_dir}/population.tsv"
log_file="${cohort}.reference_ancestry.log"
: > "$log_file"

for extension in pgen pvar psam; do
    [[ -s "${target_prefix}.${extension}" ]] || { echo "Missing target ${extension}: ${target_prefix}.${extension}" >&2; exit 2; }
    [[ -s "${reference_prefix}.${extension}" ]] || { echo "Missing ancestry reference ${extension}: ${reference_prefix}.${extension}" >&2; exit 2; }
done
[[ -s "$population" ]] || { echo "Missing reference population metadata: ${population}" >&2; exit 2; }

awk 'BEGIN{FS=OFS="\t"}
    FILENAME==ARGV[1] && !/^#/ { key=$1 FS $2 FS toupper($4) FS toupper($5); target[key]=1; next }
    FILENAME==ARGV[2] && !/^#/ { key=$1 FS $2 FS toupper($4) FS toupper($5); if(key in target) print $3 }
' "${target_prefix}.pvar" "${reference_prefix}.pvar" | sort -u > "${cohort}.ancestry_matched_variants.txt"
matched_variants=$(awk 'NF{n++} END{print n+0}' "${cohort}.ancestry_matched_variants.txt")
[[ "$matched_variants" -ge 2 ]] || { echo "Fewer than two exact target/reference variants are available for ancestry PCA." >&2; exit 3; }

if [[ "$matched_variants" -lt 100 ]]; then
    # Very small validation cohorts cannot support stable LD pruning. Retain all
    # exact common variants and make the decision visible in the log.
    cp "${cohort}.ancestry_matched_variants.txt" "${cohort}.ancestry_prune.prune.in"
    printf 'LD pruning was not applied because fewer than 100 exact variants were available.\n' >> "$log_file"
else
    plink2 --pfile "$reference_prefix" --extract "${cohort}.ancestry_matched_variants.txt" \
        --maf 0.05 --indep-pairwise 200kb 0.2 --threads "$threads" \
        --out "${cohort}.ancestry_prune" >> "$log_file" 2>&1
fi
pruned_variants=$(awk 'NF{n++} END{print n+0}' "${cohort}.ancestry_prune.prune.in")
[[ "$pruned_variants" -ge 2 ]] || { echo "Fewer than two LD-pruned variants remain for ancestry PCA." >&2; exit 3; }

plink2 --pfile "$reference_prefix" --extract "${cohort}.ancestry_prune.prune.in" \
    --freq counts --threads "$threads" \
    --out "${cohort}.reference_pca" >> "$log_file" 2>&1
plink2 --pfile "$reference_prefix" --extract "${cohort}.ancestry_prune.prune.in" \
    --read-freq "${cohort}.reference_pca.acount" \
    --pca allele-wts "$pcs" vcols=chrom,ref,alt \
    --threads "$threads" --out "${cohort}.reference_pca" >> "$log_file" 2>&1
score_end=$((5 + pcs))
for dataset in reference target; do
    if [[ "$dataset" == reference ]]; then prefix="$reference_prefix"; else prefix="$target_prefix"; fi
    plink2 --pfile "$prefix" --extract "${cohort}.ancestry_prune.prune.in" \
        --read-freq "${cohort}.reference_pca.acount" \
        --score "${cohort}.reference_pca.eigenvec.allele" 2 5 header-read no-mean-imputation variance-standardize \
        --score-col-nums "6-${score_end}" --threads "$threads" \
        --out "${cohort}.${dataset}_projection" >> "$log_file" 2>&1
done

Rscript "$classifier" \
    --cohort "$cohort" \
    --reference "${cohort}.reference_projection.sscore" \
    --target "${cohort}.target_projection.sscore" \
    --population "$population" \
    --percentile "$percentile" \
    --matched-variants "$matched_variants" \
    --pruned-variants "$pruned_variants" \
    --reference-output "${cohort}.reference_projection.tsv" \
    --target-output "${cohort}.target_ancestry.tsv" \
    --summary-output "${cohort}.ancestry_summary.tsv"
