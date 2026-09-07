#!/usr/bin/env bash
set -euo pipefail

plink2() {
    command plink2 --memory "${PLINK_MEMORY_MB:-1024}" "$@"
}

cohort="$1"
format="$2"
genotype="$3"
sample="$4"
keep="$5"
dosage="$6"
threads="$7"
input_stage="${8:-qc_completed}"
assay_manifest="${9:-}"
marker_map="${10:-}"
adapter_script="${11:-}"

mkdir -p "$cohort.imported"
imported_prefixes=()

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
        ped)
            source_path="${source_path%.ped}"
            args=(--pedmap "$source_path")
            ;;
        genomestudio)
            if [[ -z "$assay_manifest" || -z "$adapter_script" ]]; then
                echo "GenomeStudio input requires an assay manifest and target adapter." >&2
                exit 2
            fi
            local converted_prefix="${output_prefix}_genomestudio"
            perl "$adapter_script" finalreport "$source_path" "$assay_manifest" "$converted_prefix"
            args=(--pedmap "$converted_prefix")
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
    # PLINK expands these allele placeholders.
    # shellcheck disable=SC2016
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
            chromosome_prefix="$cohort.imported/${cohort}_chr${chromosome}"
            import_target "$source_path" "$chromosome_prefix" ""
            imported_prefixes+=("$chromosome_prefix")
        fi
    done
else
    all_prefix="$cohort.imported/${cohort}_all"
    import_target "$genotype" "$all_prefix" ""
    imported_prefixes+=("$all_prefix")
fi

if [[ "${#imported_prefixes[@]}" -eq 0 ]]; then
    echo "No target genotypes were prepared for cohort '$cohort'." >&2
    exit 3
fi

printf '%s\n' "${imported_prefixes[@]}" > "$cohort.imported/${cohort}_merge.txt"
if [[ "${#imported_prefixes[@]}" -eq 1 ]]; then
    cp "${imported_prefixes[0]}.pgen" "$cohort.imported/${cohort}.pgen"
    cp "${imported_prefixes[0]}.pvar" "$cohort.imported/${cohort}.pvar"
    cp "${imported_prefixes[0]}.psam" "$cohort.imported/${cohort}.psam"
else
    plink2 --pmerge-list "$cohort.imported/${cohort}_merge.txt" pfile --make-pgen \
        --threads "$threads" --out "$cohort.imported/${cohort}"
fi

if [[ "$input_stage" == "raw" ]]; then
    perl "$adapter_script" annotate-pvar \
        "$cohort.imported/$cohort.pvar" "$assay_manifest" "$marker_map" \
        "$cohort.imported/$cohort.annotated.pvar" "$cohort.initial_marker_decisions.tsv"
    mv "$cohort.imported/$cohort.annotated.pvar" "$cohort.imported/$cohort.pvar"
else
    awk -v stage="$input_stage" 'BEGIN{FS=OFS="\t"; print "source_id","final_id","source_chr","source_pos","final_chr","final_pos","source_ref","source_alt","final_ref","final_alt","decision","reason"}
        /^#/ {next}
        {print $3,$3,$1,$2,$1,$2,$4,$5,$4,$5,"INHERITED","Input entered at " stage " stage"}
    ' "$cohort.imported/$cohort.pvar" > "$cohort.initial_marker_decisions.tsv"
fi
