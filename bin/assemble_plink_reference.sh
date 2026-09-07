#!/usr/bin/env bash
set -euo pipefail

plink2() {
    command plink2 --memory "${PLINK_MEMORY_MB:-1024}" "$@"
}
output_dir="$1"
threads="$2"
genome_build="$3"
expected_chromosomes="$4"
population_panel="$5"
related_samples="$6"
shift 6
mkdir -p "$output_dir/chromosomes"
log_file="$output_dir.prepare.log"
: > "$log_file"
printf 'chromosome\tsource\tsamples\tvariants\tstatus\tsource_sha256\n' > "$output_dir.source_qc.tsv"
for group_dir in "$@"; do
    tail -n +2 "$group_dir.source_qc.tsv" >> "$output_dir.source_qc.tsv"
    cat "$group_dir.prepare.log" >> "$log_file"
done
observed=$(tail -n +2 "$output_dir.source_qc.tsv" | cut -f1 | sort -n | paste -sd, -)
[[ "$observed" == "$expected_chromosomes" ]] || {
    echo "Reference gather expected chromosomes $expected_chromosomes; received $observed. Resume the missing or failed group." >&2
    exit 4
}
prefixes=()
all_prefixes=()
IFS=',' read -ra chromosomes <<< "$expected_chromosomes"
for chromosome in "${chromosomes[@]}"; do
    for population in eur all; do
        sources=()
        for group_dir in "$@"; do
            if [[ -s "$group_dir/chromosomes/${population}_chr$chromosome.pgen" ]]; then
                sources+=("$group_dir/chromosomes/${population}_chr$chromosome")
            fi
        done
        [[ "${#sources[@]}" -eq 1 ]] || { echo "Reference chromosome $chromosome ($population) requires exactly one PGEN dataset." >&2; exit 4; }
        for extension in pgen pvar psam; do
            cp "${sources[0]}.$extension" "$output_dir/chromosomes/${population}_chr$chromosome.$extension"
        done
    done
    prefixes+=("$output_dir/chromosomes/eur_chr$chromosome")
    all_prefixes+=("$output_dir/chromosomes/all_chr$chromosome")
done
[[ "${#prefixes[@]}" -gt 0 ]] || { echo "No autosomal PLINK reference chromosomes were prepared." >&2; exit 4; }
printf '%s\n' "${prefixes[@]}" > reference.merge_list.txt
if [[ "${#prefixes[@]}" -eq 1 ]]; then
    cp "${prefixes[0]}.pgen" "$output_dir/eur_reference.pgen"
    cp "${prefixes[0]}.pvar" "$output_dir/eur_reference.pvar"
    cp "${prefixes[0]}.psam" "$output_dir/eur_reference.psam"
else
    plink2 --pmerge-list reference.merge_list.txt pfile --make-pgen \
        --threads "$threads" --out "$output_dir/eur_reference" >> "$log_file" 2>&1
fi

printf '%s\n' "${all_prefixes[@]}" > reference.all_merge_list.txt
if [[ "${#all_prefixes[@]}" -eq 1 ]]; then
    cp "${all_prefixes[0]}.pgen" "$output_dir/all_reference.pgen"
    cp "${all_prefixes[0]}.pvar" "$output_dir/all_reference.pvar"
    cp "${all_prefixes[0]}.psam" "$output_dir/all_reference.psam"
else
    plink2 --pmerge-list reference.all_merge_list.txt pfile --make-pgen \
        --threads "$threads" --out "$output_dir/all_reference" >> "$log_file" 2>&1
fi
cp "$population_panel" "$output_dir/population.tsv"

printf 'reference_type\tbuild\tancestry\tchromosomes\tsamples\tvariants\tstatus\n' > "${output_dir}.summary.tsv"
printf 'plink_ld\t%s\tEuropean\t%s\t%s\t%s\tPASS\n' \
    "$genome_build" "${#prefixes[@]}" \
    "$(awk 'BEGIN { count=0 } !/^#/ && NF { count++ } END { print count }' "$output_dir/eur_reference.psam")" \
    "$(awk 'BEGIN { count=0 } !/^#/ && NF { count++ } END { print count }' "$output_dir/eur_reference.pvar")" \
    >> "${output_dir}.summary.tsv"
printf 'ancestry_reference\t%s\tMultiple\t%s\t%s\t%s\tPASS\n' \
    "$genome_build" "${#all_prefixes[@]}" \
    "$(awk 'BEGIN { count=0 } !/^#/ && NF { count++ } END { print count }' "$output_dir/all_reference.psam")" \
    "$(awk 'BEGIN { count=0 } !/^#/ && NF { count++ } END { print count }' "$output_dir/all_reference.pvar")" \
    >> "${output_dir}.summary.tsv"

printf 'build\tpopulation\trelated_sample_rule\n%s\tEUR and all populations\tExclude listed related samples\n' "$genome_build" > "$output_dir.identity.tsv"
sha256sum "$population_panel" "$related_samples" > "$output_dir.selection_sources.sha256"
