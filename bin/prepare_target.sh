#!/usr/bin/env bash
set -euo pipefail

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
dbsnp_source="${12:-}"
reference_fasta="${13:-}"

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

if [[ "$input_stage" == "raw" ]]; then
    dbsnp_marker_map=""
    if [[ -n "$dbsnp_source" ]]; then
        if [[ -d "$dbsnp_source" ]]; then
            dbsnp_vcf=$(find "$dbsnp_source" -maxdepth 1 -type f -name 'GCF_*.gz' ! -name '*.md5' -print -quit)
            assembly_report=$(find "$dbsnp_source" -maxdepth 1 -type f -name 'assembly_report.txt' -print -quit)
        else
            dbsnp_vcf="$dbsnp_source"
            assembly_report=""
        fi
        [[ -n "$dbsnp_vcf" && -s "$dbsnp_vcf" ]] || { echo "No dbSNP VCF was found in $dbsnp_source." >&2; exit 4; }
        [[ -n "$assembly_report" && -s "$assembly_report" ]] || { echo "dbSNP assembly_report.txt is required." >&2; exit 4; }
        [[ -s "${dbsnp_vcf}.tbi" || -s "${dbsnp_vcf}.csi" ]] || { echo "The dbSNP VCF requires a tabix/CSI index." >&2; exit 4; }

        awk -F '\t' 'BEGIN { OFS="\t" }
            !/^#/ && $2 == "assembled-molecule" && $3 ~ /^([1-9]|1[0-9]|2[0-2])$/ { print $3, $7 }
        ' "$assembly_report" > chromosome_refseq.tsv
        awk -F '\t' 'BEGIN { OFS="\t" }
            NR == FNR { accession[$1]=$2; next }
            /^##/ || /^#CHROM/ { next }
            {
                chromosome=$1; sub(/^chr/, "", chromosome)
                if (chromosome in accession && $2 ~ /^[0-9]+$/) print accession[chromosome], $2, $2
            }
        ' chromosome_refseq.tsv "$cohort/${cohort}.pvar" | sort -k1,1 -k2,2n -u > dbsnp_regions.tsv

        dbsnp_marker_map="${cohort}.dbsnp_marker_map.tsv"
        printf 'source_id\tnew_id\tchr\tpos\tref\talt\n' > "$dbsnp_marker_map"
        if [[ -s dbsnp_regions.tsv ]]; then
            bcftools view -R dbsnp_regions.tsv -Ou "$dbsnp_vcf" |
                bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\n' > dbsnp_records.tsv
            awk -F '\t' 'BEGIN { OFS="\t" }
                function complement(value, result, offset, base) {
                    result=""
                    for (offset=1; offset<=length(value); offset++) {
                        base=substr(value,offset,1)
                        result=result (base=="A"?"T":base=="T"?"A":base=="C"?"G":base=="G"?"C":base)
                    }
                    return result
                }
                FILENAME == ARGV[1] { chromosome[$2]=$1; next }
                FILENAME == ARGV[2] {
                    key=chromosome[$1] ":" $2
                    if (!(key in db_id) && $3 ~ /^rs[0-9]+$/ && $4 ~ /^[ACGT]+$/ && $5 ~ /^[ACGT]+$/) {
                        db_id[key]=$3; db_ref[key]=toupper($4); db_alt[key]=toupper($5)
                    }
                    next
                }
                /^##/ || /^#CHROM/ { next }
                {
                    chr=$1; sub(/^chr/, "", chr); key=chr ":" $2
                    source_ref=toupper($4); source_alt=toupper($5)
                    if (!(key in db_id) || source_ref !~ /^[ACGT]+$/ || source_alt !~ /^[ACGT]+$/) next
                    direct=(source_ref==db_ref[key] && source_alt==db_alt[key]) || (source_ref==db_alt[key] && source_alt==db_ref[key])
                    flipped=(complement(source_ref)==db_ref[key] && complement(source_alt)==db_alt[key]) ||
                        (complement(source_ref)==db_alt[key] && complement(source_alt)==db_ref[key])
                    if (direct) print $3,db_id[key],chr,$2,source_ref,source_alt
                    else if (flipped) print $3,db_id[key],chr,$2,complement(source_ref),complement(source_alt)
                }
            ' chromosome_refseq.tsv dbsnp_records.tsv "$cohort/${cohort}.pvar" >> "$dbsnp_marker_map"
        fi
    fi

    combined_marker_map="$marker_map"
    if [[ -n "$dbsnp_marker_map" ]]; then
        combined_marker_map="${cohort}.combined_marker_map.tsv"
        if [[ -n "$marker_map" ]]; then
            awk 'NR == 1 || FNR > 1' "$marker_map" "$dbsnp_marker_map" > "$combined_marker_map"
        else
            cp "$dbsnp_marker_map" "$combined_marker_map"
        fi
    fi
    perl "$adapter_script" annotate-pvar \
        "$cohort/${cohort}.pvar" "$assay_manifest" "$combined_marker_map" \
        "$cohort/${cohort}.annotated.pvar" "$cohort.marker_decisions.tsv"
    mv "$cohort/${cohort}.annotated.pvar" "$cohort/${cohort}.pvar"

    reference_prefix="$cohort/${cohort}"
    if [[ -n "$reference_fasta" ]]; then
        plink2 --pfile "$reference_prefix" --ref-from-fa force --fa "$reference_fasta" --make-pgen \
            --threads "$threads" --out "$cohort/${cohort}_referenced"
        reference_prefix="$cohort/${cohort}_referenced"
    fi
    plink2 --pfile "$reference_prefix" --rm-dup exclude-mismatch list --make-pgen \
        --threads "$threads" --out "$cohort/${cohort}_deduplicated"
    mv "$cohort/${cohort}_deduplicated.pgen" "$cohort/${cohort}.pgen"
    mv "$cohort/${cohort}_deduplicated.pvar" "$cohort/${cohort}.pvar"
    mv "$cohort/${cohort}_deduplicated.psam" "$cohort/${cohort}.psam"
else
    awk -v stage="$input_stage" 'BEGIN{FS=OFS="\t"; print "source_id","final_id","source_chr","source_pos","final_chr","final_pos","source_ref","source_alt","final_ref","final_alt","decision","reason"}
        /^##/ || /^#CHROM/ {next}
        {print $3,$3,$1,$2,$1,$2,$4,$5,$4,$5,"INHERITED","Input entered at " stage " stage"}
    ' "$cohort/${cohort}.pvar" > "$cohort.marker_decisions.tsv"
fi

sample_count=$(awk 'BEGIN{n=0} !/^#/ && NF>0 {n++} END{print n}' "$cohort/${cohort}.psam")
variant_count=$(awk 'BEGIN{n=0} !/^#/ && NF>0 {n++} END{print n}' "$cohort/${cohort}.pvar")

printf 'cohort\tinput_stage\tstep\tparticipants\tvariants\tstatus\n%s\t%s\tNormalised to PGEN\t%s\t%s\tPASS\n' \
    "$cohort" "$input_stage" "$sample_count" "$variant_count" > "$cohort.target_prep_summary.tsv"
printf 'cohort\tinput_stage\tparticipants\tvariants\tchromosomes\tstatus\n%s\t%s\t%s\t%s\t%s\tPASS\n' \
    "$cohort" "$input_stage" "$sample_count" "$variant_count" "${#prefixes[@]}" > "$cohort.target_qc.tsv"
