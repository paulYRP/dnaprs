#!/usr/bin/env bash
set -euo pipefail

cohort="$1"
format="$2"
genotype="$3"
sample="$4"
keep="$5"
dosage="$6"
threads="$7"

mkdir -p "$cohort"
prefixes=()

import_target() {
    local source_path="$1"
    local output_prefix="$2"
    local chromosome="${3:-}"
    local args=()

    case "$format" in
        pgen)
            source_path="${source_path%.pgen}"
            args=(--pfile "$source_path")
            ;;
        bed)
            source_path="${source_path%.bed}"
            args=(--bfile "$source_path")
            ;;
        vcf)
            args=(--vcf "$source_path" "dosage=${dosage:-DS}" --double-id)
            ;;
        bgen)
            args=(--bgen "$source_path" ref-first)
            if [[ -n "$sample" ]]; then args+=(--sample "$sample"); fi
            ;;
        *)
            echo "Unsupported target format: $format" >&2
            exit 2
            ;;
    esac

    if [[ -n "$chromosome" ]]; then args+=(--chr "$chromosome"); fi
    if [[ -n "$keep" ]]; then args+=(--keep "$keep"); fi
    plink2 "${args[@]}" --set-missing-var-ids '@:#:$r:$a' --make-pgen \
        --threads "$threads" --out "$output_prefix"
}

if [[ "$genotype" == *'{chr}'* || "$genotype" == *'{CHR}'* || "$genotype" == *'{chromosome}'* || "$genotype" == *'#'* ]]; then
    for chromosome in $(seq 1 22); do
        source_path="${genotype//\{chr\}/$chromosome}"
        source_path="${source_path//\{CHR\}/$chromosome}"
        source_path="${source_path//\{chromosome\}/$chromosome}"
        source_path="${source_path//#/$chromosome}"
        if [[ "$format" == "pgen" ]]; then
            source_check="${source_path%.pgen}.pgen"
        elif [[ "$format" == "bed" ]]; then
            source_check="${source_path%.bed}.bed"
        else
            source_check="$source_path"
        fi
        if [[ -e "$source_check" ]]; then
            chromosome_prefix="$cohort/${cohort}_chr${chromosome}"
            import_target "$source_path" "$chromosome_prefix" ""
            prefixes+=("$chromosome_prefix")
        fi
    done
else
    all_prefix="$cohort/${cohort}_all"
    import_target "$genotype" "$all_prefix" ""
    mapfile -t chromosomes < <(awk '!/^#/ {print $1}' "${all_prefix}.pvar" | sed 's/^chr//' | awk '$1 >= 1 && $1 <= 22' | sort -n -u)
    for chromosome in "${chromosomes[@]}"; do
        chromosome_prefix="$cohort/${cohort}_chr${chromosome}"
        plink2 --pfile "$all_prefix" --chr "$chromosome" --make-pgen \
            --threads "$threads" --out "$chromosome_prefix"
        prefixes+=("$chromosome_prefix")
    done
fi

if [[ "${#prefixes[@]}" -eq 0 ]]; then
    echo "No autosomal target genotypes were prepared for $cohort." >&2
    exit 3
fi

printf '%s\n' "${prefixes[@]}" > "$cohort/${cohort}_merge.txt"
if [[ "${#prefixes[@]}" -eq 1 ]]; then
    cp "${prefixes[0]}.pgen" "$cohort/${cohort}.pgen"
    cp "${prefixes[0]}.pvar" "$cohort/${cohort}.pvar"
    cp "${prefixes[0]}.psam" "$cohort/${cohort}.psam"
else
    plink2 --pmerge-list "$cohort/${cohort}_merge.txt" pfile --make-pgen \
        --threads "$threads" --out "$cohort/${cohort}"
fi

sample_count=$(awk 'BEGIN{n=0} !/^#/ && NF>0 {n++} END{print n}' "$cohort/${cohort}.psam")
variant_count=$(awk 'BEGIN{n=0} !/^#/ && NF>0 {n++} END{print n}' "$cohort/${cohort}.pvar")
printf 'cohort\tparticipants\tvariants\tchromosomes\tstatus\n%s\t%s\t%s\t%s\tPASS\n' \
    "$cohort" "$sample_count" "$variant_count" "${#prefixes[@]}" > "$cohort.target_qc.tsv"
