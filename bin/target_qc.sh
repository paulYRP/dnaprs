#!/usr/bin/env bash
set -euo pipefail

cohort="$1"
target_dir="$2"
input_stage="$3"
imputation_geno="$4"
direct_geno="$5"
mind="$6"
maf="$7"
hwe="$8"
threads="$9"

source_prefix="${target_dir}/${cohort}"
imputation_dir="${cohort}.imputation_ready"
direct_dir="${cohort}.direct_ready"
imputation_prefix="${imputation_dir}/${cohort}"
direct_prefix="${direct_dir}/${cohort}"
baseline="${cohort}.baseline"
mkdir -p "$imputation_dir" "$direct_dir"

plink2 --pfile "$source_prefix" --missing --freq --hardy midp \
    --threads "$threads" --out "$baseline"

source_samples=$(awk 'BEGIN{n=0} !/^#/ && NF {n++} END{print n}' "${source_prefix}.psam")
source_variants=$(awk 'BEGIN{n=0} !/^#/ && NF {n++} END{print n}' "${source_prefix}.pvar")

is_positive() {
    awk -v value="$1" 'BEGIN { exit !((value + 0) > 0) }'
}

case "$input_stage" in
    raw|corrected)
        qc_args=(--autosome --snps-only just-acgt --mind "$mind" --geno "$imputation_geno")
        if is_positive "$maf"; then qc_args+=(--maf "$maf"); fi
        if is_positive "$hwe"; then qc_args+=(--hwe "$hwe" midp keep-fewhet); fi
        plink2 --pfile "$source_prefix" "${qc_args[@]}" \
            --make-pgen --threads "$threads" --out "$imputation_prefix"
        plink2 --pfile "$imputation_prefix" --geno "$direct_geno" \
            --make-pgen --threads "$threads" --out "$direct_prefix"
        qc_action="FILTERED"
        ;;
    qc_completed|imputed)
        plink2 --pfile "$source_prefix" --autosome --make-pgen \
            --threads "$threads" --out "$imputation_prefix"
        plink2 --pfile "$imputation_prefix" --make-pgen \
            --threads "$threads" --out "$direct_prefix"
        qc_action="INHERITED"
        ;;
    *)
        echo "Unsupported input stage: $input_stage" >&2
        exit 2
        ;;
esac

retained_samples=$(awk 'BEGIN{n=0} !/^#/ && NF {n++} END{print n}' "${imputation_prefix}.psam")
imputation_variants=$(awk 'BEGIN{n=0} !/^#/ && NF {n++} END{print n}' "${imputation_prefix}.pvar")
direct_variants=$(awk 'BEGIN{n=0} !/^#/ && NF {n++} END{print n}' "${direct_prefix}.pvar")
if [[ "$retained_samples" -eq 0 || "$imputation_variants" -eq 0 || "$direct_variants" -eq 0 ]]; then
    echo "Target QC removed every participant or variant for ${cohort}." >&2
    exit 3
fi

awk -v cohort="$cohort" -v stage="$input_stage" -v threshold="$mind" '
BEGIN{FS=OFS="\t"; print "cohort","FID","IID","missingness","decision","reason"}
/^#FID/ {for(i=1;i<=NF;i++) ix[$i]=i; next}
!/^#/ {
    fid=$(ix["#FID"]); iid=$(ix["IID"]); missing=$(ix["F_MISS"])+0
    if(stage=="raw" || stage=="corrected") {
        decision=(missing>threshold ? "FILTER" : "RETAIN")
        reason=(missing>threshold ? "Sample missingness exceeds sample_missingness" : "Sample missingness within threshold")
    } else {decision="INHERITED"; reason="Participant eligibility inherited from declared checkpoint"}
    print cohort,fid,iid,missing,decision,reason
}' "${baseline}.smiss" > "${cohort}.sample_decisions.tsv"

awk -v cohort="$cohort" -v stage="$input_stage" -v impute_geno="$imputation_geno" -v direct_geno="$direct_geno" -v mafmin="$maf" -v hwemin="$hwe" '
BEGIN{FS=OFS="\t"; print "cohort","ID","chromosome","position","ref","alt","missingness","maf","hwe_midp","imputation_decision","direct_decision","decision","reason"}
FILENAME ~ /\.vmiss$/ {
    if(/^#CHROM/){delete ix; for(i=1;i<=NF;i++) ix[$i]=i; next}
    if(!/^#/){missing[$(ix["ID"]) ]=$(ix["F_MISS"])+0}
    next
}
FILENAME ~ /\.afreq$/ {
    if(/^#CHROM/){delete ix; for(i=1;i<=NF;i++) ix[$i]=i; next}
    if(!/^#/){v=$(ix["ALT_FREQS"])+0; frequency[$(ix["ID"])]=(v<=0.5?v:1-v)}
    next
}
FILENAME ~ /\.hardy$/ {
    if(/^#CHROM/){delete ix; for(i=1;i<=NF;i++) ix[$i]=i; next}
    if(!/^#/){
        pcol=("MIDP" in ix ? ix["MIDP"] : ix["P"])
        hardy[$(ix["ID"]) ]=(pcol ? $(pcol)+0 : "NA")
    }
    next
}
/^##/ || /^#CHROM/ {next}
{
    chr=$1; pos=$2; id=$3; ref=toupper($4); alt=toupper($5)
    miss=(id in missing ? missing[id] : "NA")
    frq=(id in frequency ? frequency[id] : "NA")
    hp=(id in hardy ? hardy[id] : "NA")
    if(stage!="raw" && stage!="corrected") {imp="INHERITED"; direct="INHERITED"; reason="Variant eligibility inherited from declared checkpoint"}
    else if(chr !~ /^([1-9]|1[0-9]|2[0-2])$/) {imp="FILTER"; direct="FILTER"; reason="Non-autosomal variant"}
    else if(ref !~ /^[ACGT]$/ || alt !~ /^[ACGT]$/) {imp="FILTER"; direct="FILTER"; reason="Non-biallelic or unresolved SNP alleles"}
    else if(frq!="NA" && mafmin>0 && frq<mafmin) {imp="FILTER"; direct="FILTER"; reason="Minor allele frequency below maf_filter"}
    else if(hp!="NA" && hwemin>0 && hp<hwemin) {imp="FILTER"; direct="FILTER"; reason="Hardy-Weinberg mid-P below hwe_filter"}
    else if(miss!="NA" && miss>impute_geno) {imp="FILTER"; direct="FILTER"; reason="Variant missingness exceeds imputation checkpoint"}
    else if(miss!="NA" && miss>direct_geno) {imp="RETAIN"; direct="FILTER"; reason="Retained for imputation but exceeds direct-scoring checkpoint"}
    else {imp="RETAIN"; direct="RETAIN"; reason="Variant passes both technical checkpoints"}
    print cohort,id,chr,pos,ref,alt,miss,frq,hp,imp,direct,direct,reason
}' "${baseline}.vmiss" "${baseline}.afreq" "${baseline}.hardy" "${source_prefix}.pvar" \
    > "${cohort}.variant_decisions.tsv"

mapfile -t chromosomes < <(awk '!/^#/ {print $1}' "${imputation_prefix}.pvar" | sed 's/^chr//' | awk '$1>=1 && $1<=22' | sort -n -u)
for chromosome in "${chromosomes[@]}"; do
    plink2 --pfile "$imputation_prefix" --chr "$chromosome" --make-pgen \
        --threads "$threads" --out "${imputation_dir}/${cohort}_chr${chromosome}"
    plink2 --pfile "$direct_prefix" --chr "$chromosome" --make-pgen \
        --threads "$threads" --out "${direct_dir}/${cohort}_chr${chromosome}"
done

printf 'cohort\tinput_stage\tqc_action\tsource_participants\tretained_participants\tfiltered_participants\tsource_variants\timputation_variants\tretained_variants\tfiltered_variants\timputation_variant_missingness\tdirect_variant_missingness\tsample_missingness\tmaf_filter\thwe_filter\tchromosomes\tstatus\n' \
    > "${cohort}.target_qc.tsv"
printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\tPASS\n' \
    "$cohort" "$input_stage" "$qc_action" "$source_samples" "$retained_samples" "$((source_samples-retained_samples))" \
    "$source_variants" "$imputation_variants" "$direct_variants" "$((source_variants-direct_variants))" \
    "$imputation_geno" "$direct_geno" "$mind" "$maf" "$hwe" "${#chromosomes[@]}" >> "${cohort}.target_qc.tsv"
