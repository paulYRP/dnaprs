#!/usr/bin/env bash
set -euo pipefail

plink2() {
    command plink2 --memory "${PLINK_MEMORY_MB:-1024}" "$@"
}

panel="$1"
population_panel="$2"
related_samples="$3"
output_dir="$4"
threads="$5"
genome_build="$6"
expected_chromosomes="$7"

unbref3_jar="${UNBREF3_JAR:-/opt/beagle/unbref3.jar}"
log_file="${output_dir}.prepare.log"
mkdir -p "$output_dir/chromosomes"
: > "$log_file"

[[ "$genome_build" == "GRCh37" ]] || {
    echo "The downloaded 1000 Genomes source bundle is declared for GRCh37, not ${genome_build}." >&2
    exit 2
}
[[ -r "$population_panel" ]] || { echo "Population panel is unreadable: $population_panel" >&2; exit 2; }
[[ -r "$related_samples" ]] || { echo "Related-sample list is unreadable: $related_samples" >&2; exit 2; }

awk -F '\t' '
    NR == 1 {
        for (column = 1; column <= NF; column++) {
            name = tolower($column)
            if (name == "sample") sample_column = column
            if (name == "super_pop") population_column = column
        }
        if (!sample_column || !population_column) exit 2
        next
    }
    $population_column == "EUR" { print $sample_column }
' "$population_panel" | sort -u > eur.samples.all.txt

awk -F '\t' '
    NR == 1 {
        for (column = 1; column <= NF; column++) if (tolower($column) == "sample") sample_column = column
        if (!sample_column) exit 2
        next
    }
    $sample_column != "" { print $sample_column }
' "$population_panel" | sort -u > reference.samples.all.txt

awk -F '[\t ]+' 'NR > 1 && $1 != "" { print $1 }' "$related_samples" | sort -u > related.samples.txt
grep -Fvx -f related.samples.txt eur.samples.all.txt > eur.samples.txt || true
grep -Fvx -f related.samples.txt reference.samples.all.txt > unrelated.samples.txt || true
[[ -s eur.samples.txt ]] || { echo "No unrelated European samples were selected." >&2; exit 3; }
[[ -s unrelated.samples.txt ]] || { echo "No unrelated reference samples were selected." >&2; exit 3; }

mapfile -t panel_files < "$panel"
[[ "${#panel_files[@]}" -gt 0 ]] || { echo "No chromosome reference files were found in $panel." >&2; exit 3; }

prefixes=()
all_prefixes=()
printf 'chromosome\tsource\tsamples\tvariants\tstatus\tsource_sha256\n' > "${output_dir}.source_qc.tsv"
for source_file in "${panel_files[@]}"; do
    base_name=$(basename "$source_file")
    chromosome=$(sed -nE 's/^chr([0-9]+).*/\1/p' <<< "$base_name")
    [[ "$chromosome" =~ ^([1-9]|1[0-9]|2[0-2])$ ]] || continue

    source_vcf="chr${chromosome}.source.vcf.gz"
    if [[ "$source_file" == *.bref3 ]]; then
        [[ -s "$unbref3_jar" ]] || { echo "unbref3 JAR was not found: $unbref3_jar" >&2; exit 2; }
        java -jar "$unbref3_jar" "$source_file" 2>> "$log_file" | bgzip -c > "$source_vcf"
    else
        bcftools view -Oz -o "$source_vcf" "$source_file" 2>> "$log_file"
    fi
    tabix -f -p vcf "$source_vcf"

    filtered_vcf="chr${chromosome}.eur.vcf.gz"
    bcftools view --force-samples --samples-file eur.samples.txt \
        --min-alleles 2 --max-alleles 2 --types snps -Ou "$source_vcf" 2>> "$log_file" |
        bcftools annotate --set-id '%CHROM:%POS:%REF:%FIRST_ALT' -Oz -o "$filtered_vcf" 2>> "$log_file"
    tabix -f -p vcf "$filtered_vcf"

    all_vcf="chr${chromosome}.unrelated.vcf.gz"
    bcftools view --force-samples --samples-file unrelated.samples.txt \
        --min-alleles 2 --max-alleles 2 --types snps -Ou "$source_vcf" 2>> "$log_file" |
        bcftools annotate --set-id '%CHROM:%POS:%REF:%FIRST_ALT' -Oz -o "$all_vcf" 2>> "$log_file"
    tabix -f -p vcf "$all_vcf"

    chromosome_prefix="${output_dir}/chromosomes/eur_chr${chromosome}"
    # PLINK expands these allele placeholders.
    # shellcheck disable=SC2016
    plink2 --vcf "$filtered_vcf" --double-id --set-all-var-ids '@:#:$r:$a' \
        --make-pgen --threads "$threads" --out "$chromosome_prefix" >> "$log_file" 2>&1
    prefixes+=("$chromosome_prefix")
    all_chromosome_prefix="${output_dir}/chromosomes/all_chr${chromosome}"
    # shellcheck disable=SC2016
    plink2 --vcf "$all_vcf" --double-id --set-all-var-ids '@:#:$r:$a' \
        --make-pgen --threads "$threads" --out "$all_chromosome_prefix" >> "$log_file" 2>&1
    all_prefixes+=("$all_chromosome_prefix")

    sample_count=$(awk 'BEGIN { count=0 } !/^#/ && NF { count++ } END { print count }' "${chromosome_prefix}.psam")
    variant_count=$(awk 'BEGIN { count=0 } !/^#/ && NF { count++ } END { print count }' "${chromosome_prefix}.pvar")
    printf '%s\t%s\t%s\t%s\tPASS\t%s\n' "$chromosome" "$base_name" "$sample_count" "$variant_count" "$(sha256sum "$source_file" | cut -d ' ' -f1)" \
        >> "${output_dir}.source_qc.tsv"
done


observed=$(cut -f1 "$output_dir.source_qc.tsv" | tail -n +2 | sort -n | paste -sd, -)
[[ "$expected_chromosomes" == "$observed" ]] || {
    echo "Reference group expected chromosomes $expected_chromosomes; prepared $observed." >&2
    exit 4
}
