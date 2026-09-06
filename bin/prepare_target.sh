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
marker_resolver="${12:-}"
dbsnp_source="${13:-}"
reference_fasta="${14:-}"

mkdir -p "$cohort"
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
            imported_prefixes+=("$chromosome_prefix")
        fi
    done
else
    all_prefix="$cohort/${cohort}_all"
    import_target "$genotype" "$all_prefix" ""
    imported_prefixes+=("$all_prefix")
fi

if [[ "${#imported_prefixes[@]}" -eq 0 ]]; then
    echo "No target genotypes were prepared for cohort '$cohort'." >&2
    exit 3
fi

printf '%s\n' "${imported_prefixes[@]}" > "$cohort/${cohort}_merge.txt"
if [[ "${#imported_prefixes[@]}" -eq 1 ]]; then
    cp "${imported_prefixes[0]}.pgen" "$cohort/${cohort}.pgen"
    cp "${imported_prefixes[0]}.pvar" "$cohort/${cohort}.pvar"
    cp "${imported_prefixes[0]}.psam" "$cohort/${cohort}.psam"
else
    plink2 --pmerge-list "$cohort/${cohort}_merge.txt" pfile --make-pgen \
        --threads "$threads" --out "$cohort/${cohort}"
fi

if [[ "$input_stage" == "raw" ]]; then
    [[ -n "$dbsnp_source" && -n "$reference_fasta" && -n "$marker_resolver" ]] || {
        echo "Raw input requires dbSNP, GRCh37 FASTA, and the marker-resolution script." >&2
        exit 4
    }
    if [[ -d "$dbsnp_source" ]]; then
        mapfile -t dbsnp_candidates < <(
            find -L "$dbsnp_source" -maxdepth 1 -type f \
                \( -iname '*.vcf.gz' -o -iname '*.vcf.bgz' -o -iname '*.bcf' -o -iname '*.vcf' -o -iname 'GCF_*.gz' \) \
                ! -name '*.md5' -print | sort
        )
        mapfile -t assembly_reports < <(
            find -L "$dbsnp_source" -maxdepth 1 -type f \
                \( -iname 'assembly_report.txt' -o -iname '*assembly_report*.txt' -o -iname '*assembly-report*.txt' \) \
                -print | sort
        )
        [[ "${#dbsnp_candidates[@]}" -eq 1 ]] || {
            echo "Expected exactly one dbSNP VCF in $dbsnp_source; found ${#dbsnp_candidates[@]}." >&2
            exit 4
        }
        [[ "${#assembly_reports[@]}" -eq 1 ]] || {
            echo "Expected exactly one dbSNP assembly report in $dbsnp_source; found ${#assembly_reports[@]}." >&2
            exit 4
        }
        dbsnp_vcf="${dbsnp_candidates[0]}"
        assembly_report="${assembly_reports[0]}"
    else
        dbsnp_vcf="$dbsnp_source"
        assembly_report=""
    fi
    [[ -n "$dbsnp_vcf" && -s "$dbsnp_vcf" ]] || { echo "No dbSNP VCF was found in $dbsnp_source." >&2; exit 4; }
    [[ -n "$assembly_report" && -s "$assembly_report" ]] || { echo "dbSNP assembly_report.txt is required." >&2; exit 4; }
    if [[ ! -s "${dbsnp_vcf}.tbi" && ! -s "${dbsnp_vcf}.csi" ]]; then
        echo "dbSNP has no usable index; creating a task-local bgzip copy." >&2
        bcftools view -Oz -o dbsnp.indexed.vcf.gz "$dbsnp_vcf"
        tabix -f -p vcf dbsnp.indexed.vcf.gz
        dbsnp_vcf="dbsnp.indexed.vcf.gz"
    fi

    # First apply assay-manifest recovery and any explicit marker map while
    # retaining one row per original probe. dbSNP candidate selection follows below.
    perl "$adapter_script" annotate-pvar \
        "$cohort/${cohort}.pvar" "$assay_manifest" "$marker_map" \
        "$cohort/${cohort}.annotated.pvar" "$cohort.initial_marker_decisions.tsv"
    mv "$cohort/${cohort}.annotated.pvar" "$cohort/${cohort}.pvar"

    awk -F '\t' 'BEGIN { OFS="\t" }
        !/^#/ && $2 == "assembled-molecule" && $3 ~ /^([1-9]|1[0-9]|2[0-2]|X|Y|MT)$/ { print $3, $7 }
    ' "$assembly_report" > chromosome_refseq.tsv
    awk -F '\t' 'BEGIN { OFS="\t" }
        NR == FNR { accession[$1]=$2; next }
        /^##/ || /^#CHROM/ { next }
        {
            chromosome=$1; sub(/^chr/, "", chromosome)
            if (chromosome in accession && $2 ~ /^[0-9]+$/) print accession[chromosome], $2, $2
        }
    ' chromosome_refseq.tsv "$cohort/${cohort}.pvar" | sort -k1,1 -k2,2n -u > dbsnp_regions.tsv
    : > dbsnp_records.tsv
    if [[ -s dbsnp_regions.tsv ]]; then
        bcftools view -R dbsnp_regions.tsv -Ou "$dbsnp_vcf" |
            bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\n' > dbsnp_records.tsv
    fi

    plink2 --pfile "$cohort/${cohort}" --export A-transpose \
        --threads "$threads" --out "$cohort/${cohort}.marker_calls"
    Rscript "$marker_resolver" \
        --pvar "$cohort/${cohort}.pvar" \
        --calls "$cohort/${cohort}.marker_calls.traw" \
        --dbsnp-records dbsnp_records.tsv \
        --chromosome-map chromosome_refseq.tsv \
        --initial-decisions "$cohort.initial_marker_decisions.tsv" \
        --output-decisions "$cohort.marker_decisions.tsv" \
        --keep "$cohort.retained_markers.txt" \
        --rename "$cohort.rename_markers.tsv"

    source_count=$(wc -l < "$cohort.retained_markers.txt")
    rename_count=$(wc -l < "$cohort.rename_markers.tsv")
    source_unique=$(sort -u "$cohort.retained_markers.txt" | wc -l)
    final_unique=$(cut -f2 "$cohort.rename_markers.tsv" | sort -u | wc -l)
    [[ "$source_count" -gt 0 && "$source_count" -eq "$rename_count" ]] || {
        echo "Marker resolution for cohort '$cohort' produced different retained and rename counts." >&2
        exit 5
    }
    [[ "$source_count" -eq "$source_unique" ]] || {
        echo "Marker resolution for cohort '$cohort' produced repeated retained source identifiers." >&2
        exit 5
    }
    [[ "$rename_count" -eq "$final_unique" ]] || {
        echo "Marker resolution for cohort '$cohort' produced repeated final identifiers." >&2
        exit 5
    }

    plink2 --pfile "$cohort/${cohort}" \
        --extract "$cohort.retained_markers.txt" \
        --make-pgen --threads "$threads" --out "$cohort/${cohort}_retained"
    plink2 --pfile "$cohort/${cohort}_retained" \
        --update-name "$cohort.rename_markers.tsv" \
        --make-pgen --threads "$threads" --out "$cohort/${cohort}_resolved"

    awk -F '\t' 'NR == 1 {
            for (column = 1; column <= NF; column++) {
                if ($column == "final_id") final_id = column
                if ($column == "decision") decision = column
            }
            next
        }
        $decision ~ /^RETAINED_/ { print $final_id }
    ' "$cohort.marker_decisions.tsv" | sort > expected_final_ids.txt
    awk '!/^#/ { print $3 }' "$cohort/${cohort}_resolved.pvar" | sort > observed_final_ids.txt
    if ! cmp -s expected_final_ids.txt observed_final_ids.txt; then
        echo "Resolved PVAR identifiers for cohort '$cohort' do not match the retained marker decisions." >&2
        exit 5
    fi

    cp "$cohort/${cohort}_resolved.psam" resolved.psam
    plink2 --pfile "$cohort/${cohort}_resolved" --export vcf bgz id-paste=iid \
        --threads "$threads" --out "$cohort/${cohort}_top"
    if [[ -n "$assay_manifest" ]]; then
        bcftools +fixref "$cohort/${cohort}_top.vcf.gz" -Oz -o "$cohort/${cohort}_forward.vcf.gz" -- \
            -f "$reference_fasta" -m top
    else
        cp "$cohort/${cohort}_top.vcf.gz" "$cohort/${cohort}_forward.vcf.gz"
    fi
    bcftools norm -f "$reference_fasta" -c e -Oz \
        -o "$cohort/${cohort}_forward.checked.vcf.gz" "$cohort/${cohort}_forward.vcf.gz"
    tabix -f -p vcf "$cohort/${cohort}_forward.checked.vcf.gz"
    plink2 --vcf "$cohort/${cohort}_forward.checked.vcf.gz" --double-id \
        --make-pgen --threads "$threads" --out "$cohort/${cohort}_referenced"

    awk '!/^#/ { print $2 }' resolved.psam > expected_iids.txt
    awk '!/^#/ { print $2 }' "$cohort/${cohort}_referenced.psam" > observed_iids.txt
    if ! cmp -s expected_iids.txt observed_iids.txt; then
        echo "Participant identifiers changed during GRCh37 orientation for cohort '$cohort'." >&2
        exit 5
    fi
    cp resolved.psam "$cohort/${cohort}_referenced.psam"
    mv "$cohort/${cohort}_referenced.pgen" "$cohort/${cohort}.pgen"
    mv "$cohort/${cohort}_referenced.pvar" "$cohort/${cohort}.pvar"
    mv "$cohort/${cohort}_referenced.psam" "$cohort/${cohort}.psam"
    duplicate_count=$(awk '!/^#/ {count[$3]++} END {n=0; for (id in count) if (count[id] > 1) n++; print n}' "$cohort/${cohort}.pvar")
    [[ "$duplicate_count" == "0" ]] || { echo "Resolved target for cohort '$cohort' contains repeated marker identifiers." >&2; exit 5; }
    invalid_alleles=$(awk '!/^#/ && ($4 !~ /^[ACGT]$/ || $5 !~ /^[ACGT]$/ || $4 == $5) { count++ } END { print count + 0 }' "$cohort/${cohort}.pvar")
    [[ "$invalid_alleles" == "0" ]] || {
        echo "Resolved target for cohort '$cohort' contains $invalid_alleles invalid REF and ALT allele pairs." >&2
        exit 5
    }
else
    awk -v stage="$input_stage" 'BEGIN{FS=OFS="\t"; print "source_id","final_id","source_chr","source_pos","final_chr","final_pos","source_ref","source_alt","final_ref","final_alt","decision","reason"}
        /^##/ || /^#CHROM/ {next}
        {print $3,$3,$1,$2,$1,$2,$4,$5,$4,$5,"INHERITED","Input entered at " stage " stage"}
    ' "$cohort/${cohort}.pvar" > "$cohort.marker_decisions.tsv"
fi

prefixes=()
mapfile -t chromosomes < <(awk '!/^#/ {print $1}' "$cohort/${cohort}.pvar" | sed 's/^chr//' | awk '$1 >= 1 && $1 <= 22' | sort -n -u)
for chromosome in "${chromosomes[@]}"; do
    chromosome_prefix="$cohort/${cohort}_chr${chromosome}"
    plink2 --pfile "$cohort/${cohort}" --chr "$chromosome" --make-pgen \
        --threads "$threads" --out "$chromosome_prefix"
    prefixes+=("$chromosome_prefix")
done
if [[ "${#prefixes[@]}" -eq 0 ]]; then
    echo "Prepared target for cohort '$cohort' contains no autosomal variants." >&2
    exit 3
fi

sample_count=$(awk 'BEGIN{n=0} !/^#/ && NF>0 {n++} END{print n}' "$cohort/${cohort}.psam")
variant_count=$(awk 'BEGIN{n=0} !/^#/ && NF>0 {n++} END{print n}' "$cohort/${cohort}.pvar")

printf 'cohort\tinput_stage\tstep\tparticipants\tvariants\tstatus\n%s\t%s\tNormalised to PGEN\t%s\t%s\tPASS\n' \
    "$cohort" "$input_stage" "$sample_count" "$variant_count" > "$cohort.target_prep_summary.tsv"
printf 'cohort\tinput_stage\tparticipants\tvariants\tchromosomes\tstatus\n%s\t%s\t%s\t%s\t%s\tPASS\n' \
    "$cohort" "$input_stage" "$sample_count" "$variant_count" "${#prefixes[@]}" > "$cohort.target_qc.tsv"
