#!/usr/bin/env bash
set -euo pipefail

plink2() {
    command plink2 --memory "${PLINK_MEMORY_MB:-1024}" "$@"
}
cohort="$1"
imported="$2"
initial="$3"
input_stage="$4"
dbsnp_source="$5"
marker_resolver="$6"
threads="$7"

if [[ "$input_stage" != "raw" ]]; then
    cp "$initial" "$cohort.marker_decisions.tsv"
    : > "$cohort.retained_markers.txt"
    : > "$cohort.rename_markers.tsv"
    : > "$cohort.duplicate_markers.txt"
    exit 0
fi
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
' chromosome_refseq.tsv "$imported/${cohort}.pvar" | sort -k1,1 -k2,2n -u > dbsnp_regions.tsv
: > dbsnp_records.tsv
if [[ -s dbsnp_regions.tsv ]]; then
    bcftools view -R dbsnp_regions.tsv -Ou "$dbsnp_vcf" |
        bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\n' > dbsnp_records.tsv
fi

plink2 --pfile "$imported/$cohort" --missing variant-only vcols=nmissdosage,nobs \
    --threads "$threads" --out "$cohort.missingness"
sample_count=$(awk '!/^#/ && NF {n++} END {print n+0}' "$imported/$cohort.psam")
Rscript "$marker_resolver" --action match --threads "$threads" \
    --pvar "$imported/$cohort.pvar" --initial-decisions "$initial" \
    --missingness "$cohort.missingness.vmiss" --sample-count "$sample_count" \
    --dbsnp-records dbsnp_records.tsv --chromosome-map chromosome_refseq.tsv \
    --candidates "$cohort.marker_candidates.rds" --duplicate-markers "$cohort.duplicate_markers.txt"
if [[ -s "$cohort.duplicate_markers.txt" ]]; then
    plink2 --pfile "$imported/$cohort" --extract "$cohort.duplicate_markers.txt" \
        --export A-transpose --threads "$threads" --out "$cohort.duplicate_calls"
fi
Rscript "$marker_resolver" --action finalise --threads "$threads" \
    --candidates "$cohort.marker_candidates.rds" --calls "$cohort.duplicate_calls.traw" \
    --output-decisions "$cohort.marker_decisions.tsv" \
    --keep "$cohort.retained_markers.txt" --rename "$cohort.rename_markers.tsv"
