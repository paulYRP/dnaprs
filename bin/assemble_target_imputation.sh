#!/usr/bin/env bash
set -euo pipefail

cohort="$1"
role="$2"
build="$3"
ancestry="$4"
threads="$5"
output_dir="${cohort}.imputed"
mkdir -p "$output_dir"

mapfile -t chromosome_dirs < <(find . -maxdepth 1 -type d -name "${cohort}.chr*.imputed" -printf '%f\n' | sort -V)
[[ "${#chromosome_dirs[@]}" -gt 0 ]] || { echo "No completed chromosome-imputation bundles were supplied." >&2; exit 2; }

first_psam=""
prefixes=()
for directory in "${chromosome_dirs[@]}"; do
    mapfile -t pgen < <(find "$directory" -maxdepth 1 -type f -name "${cohort}_chr*.pgen")
    [[ "${#pgen[@]}" -eq 1 ]] || { echo "Imputation bundle ${directory} does not contain exactly one PGEN." >&2; exit 2; }
    prefix="${pgen[0]%.pgen}"
    [[ -s "${prefix}.pvar" && -s "${prefix}.psam" ]] || { echo "Incomplete PGEN bundle: ${prefix}" >&2; exit 2; }
    if [[ -z "$first_psam" ]]; then
        first_psam="${prefix}.psam"
    else
        cmp -s "$first_psam" "${prefix}.psam" || { echo "Chromosome bundles contain different participants or sample order." >&2; exit 3; }
    fi
    cp "${prefix}.pgen" "${prefix}.pvar" "${prefix}.psam" "$output_dir/"
    cp "${prefix}.vcf.gz" "${prefix}.vcf.gz.tbi" "$output_dir/"
    prefixes+=("${output_dir}/$(basename "$prefix")")
done

printf '%s\n' "${prefixes[@]}" > "${output_dir}/${cohort}_merge.txt"
plink2 --pmerge-list "${output_dir}/${cohort}_merge.txt" pfile --make-pgen \
    --threads "$threads" --out "${output_dir}/${cohort}" > "${cohort}.target_imputation.log" 2>&1

combine_tables() {
    local pattern="$1" output="$2"
    mapfile -t tables < <(find . -maxdepth 1 -type f -name "$pattern" -printf '%f\n' | sort -V)
    [[ "${#tables[@]}" -eq "${#chromosome_dirs[@]}" ]] || { echo "Incomplete imputation summary set for ${pattern}." >&2; exit 4; }
    head -n 1 "${tables[0]}" > "$output"
    for table in "${tables[@]}"; do tail -n +2 "$table"; done | sort -t $'\t' -k2,2n >> "$output"
}
combine_tables "${cohort}.chr*.imputation_manifest.tsv" "${cohort}.imputation_manifest.tsv"
combine_tables "${cohort}.chr*.imputation_qc.tsv" "${cohort}.imputation_qc.tsv"
combine_tables "${cohort}.chr*.imputation_dr2.tsv" "${cohort}.imputation_dr2.tsv"

for log in ${cohort}.chr*.target_imputation.log; do
    printf '\n===== %s =====\n' "$log" >> "${cohort}.target_imputation.log"
    cat "$log" >> "${cohort}.target_imputation.log"
done

printf 'cohort\trole\tsource_format\tgenotype\tsample\tkeep\tbuild\tancestry\tdosage\tinput_stage\tassay_manifest\tmarker_map\n' > "${cohort}.imputed_target_manifest.tsv"
printf '%s\t%s\tpgen\t%s\t\t\t%s\t%s\tDS\timputed\t\t\n' "$cohort" "$role" \
    "checkpoints/imputed/${cohort}.imputed/${cohort}.pgen" "$build" "$ancestry" >> "${cohort}.imputed_target_manifest.tsv"
