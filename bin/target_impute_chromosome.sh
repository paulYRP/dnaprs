#!/usr/bin/env bash
set -euo pipefail

cohort="$1"
chromosome="$2"
target_pgen="$3"
target_pvar="$4"
target_psam="$5"
panel="$6"
genetic_map="$7"
dr2="$8"
threads="$9"
memory_mb="${10}"
genome_build="${11}"
reference_fasta="${12}"
output_dir="${cohort}.chr${chromosome}.imputed"
log_file="${cohort}.chr${chromosome}.target_imputation.log"
mkdir -p "$output_dir"
: > "$log_file"

beagle_jar="${BEAGLE_JAR:-/opt/beagle/beagle.jar}"
[[ -s "$beagle_jar" ]] || { echo "Beagle JAR was not found: ${beagle_jar}" >&2; exit 2; }
map_root="$genetic_map"
if [[ -f "$genetic_map" && "$genetic_map" == *.zip ]]; then
    mkdir -p genetic_maps
    unzip -q "$genetic_map" -d genetic_maps
    map_root="genetic_maps"
fi

reference_vcf=""
reference_candidates=("${panel}/chr${chromosome}.vcf.gz" "${panel}/panel/chr${chromosome}.vcf.gz" "${panel}/chr${chromosome}.bcf" "${panel}/chr${chromosome}.bref3" "${panel}/panel/chr${chromosome}.bref3")
shopt -s nullglob
reference_candidates+=("${panel}/chr${chromosome}."*.bref3 "${panel}/panel/chr${chromosome}."*.bref3)
shopt -u nullglob
for candidate in "${reference_candidates[@]}"; do
    if [[ -s "$candidate" ]]; then reference_vcf="$candidate"; break; fi
done
[[ -n "$reference_vcf" ]] || { echo "No imputation reference for chromosome ${chromosome} in ${panel}" >&2; exit 3; }

map_arg=()
for candidate in "${map_root}/plink.chr${chromosome}.GRCh37.map" "${map_root}/maps/plink.chr${chromosome}.GRCh37.map" "${map_root}/chr${chromosome}.map"; do
    if [[ -s "$candidate" ]]; then map_arg=("map=${candidate}"); break; fi
done

target_prefix="${cohort}.chr${chromosome}.input"
cp "$target_pgen" "${target_prefix}.pgen"
cp "$target_pvar" "${target_prefix}.pvar"
cp "$target_psam" "${target_prefix}.psam"
plink2 --pfile "$target_prefix" --export vcf bgz id-paste=iid --threads "$threads" --out "${cohort}.chr${chromosome}.typed" >> "$log_file" 2>&1
tabix -f -p vcf "${cohort}.chr${chromosome}.typed.vcf.gz"
bcftools norm --check-ref e --fasta-ref "$reference_fasta" -Ou "${cohort}.chr${chromosome}.typed.vcf.gz" >/dev/null
input_variants=$(bcftools index --nrecords "${cohort}.chr${chromosome}.typed.vcf.gz")

java -Xmx"${memory_mb}"m -jar "$beagle_jar" "gt=${cohort}.chr${chromosome}.typed.vcf.gz" \
    "ref=${reference_vcf}" "out=${cohort}.chr${chromosome}.beagle" "nthreads=${threads}" "${map_arg[@]}" >> "$log_file" 2>&1

bcftools query -l "${cohort}.chr${chromosome}.typed.vcf.gz" > typed.samples
bcftools query -l "${cohort}.chr${chromosome}.beagle.vcf.gz" > imputed.samples
cmp -s typed.samples imputed.samples || { echo "Beagle changed sample identity or order on chromosome ${chromosome}." >&2; exit 5; }
unexpected_chromosome=$(bcftools query -f '%CHROM\n' "${cohort}.chr${chromosome}.beagle.vcf.gz" | sed 's/^chr//' | awk -v expected="$chromosome" '$1 != expected {print; exit}')
[[ -z "$unexpected_chromosome" ]] || { echo "Beagle output contains the wrong chromosome: ${unexpected_chromosome}" >&2; exit 5; }
duplicate_key=$(bcftools query -f '%CHROM:%POS:%REF:%ALT\n' "${cohort}.chr${chromosome}.beagle.vcf.gz" | sort | uniq -d | head -n 1)
[[ -z "$duplicate_key" ]] || { echo "Beagle output contains a duplicate variant key: ${duplicate_key}" >&2; exit 5; }
invalid_dosage=$(bcftools query -f '[\t%DS]\n' "${cohort}.chr${chromosome}.beagle.vcf.gz" | awk -F '[\t,]' '{for(i=1;i<=NF;i++) if($i!="." && $i!="" && ($i+0<0 || $i+0>2)){print $i; exit}}')
[[ -z "$invalid_dosage" ]] || { echo "Beagle output contains a dosage outside [0,2]: ${invalid_dosage}" >&2; exit 5; }

printf 'cohort\tchromosome\tdr2_bin\tvariants\n' > "${cohort}.chr${chromosome}.imputation_dr2.tsv"
bcftools query -i 'INFO/IMP=1' -f '%INFO/DR2\n' "${cohort}.chr${chromosome}.beagle.vcf.gz" |
    awk -v cohort="$cohort" -v chromosome="$chromosome" 'BEGIN{FS=OFS="\t"} {value=$1+0; if(value<0.3) bin="[0,0.3)"; else if(value<0.5) bin="[0.3,0.5)"; else if(value<0.8) bin="[0.5,0.8)"; else if(value<0.9) bin="[0.8,0.9)"; else bin="[0.9,1]"; count[bin]++} END{order[1]="[0,0.3)"; order[2]="[0.3,0.5)"; order[3]="[0.5,0.8)"; order[4]="[0.8,0.9)"; order[5]="[0.9,1]"; for(i=1;i<=5;i++) print cohort,chromosome,order[i],count[order[i]]+0}' \
    >> "${cohort}.chr${chromosome}.imputation_dr2.tsv"

bcftools view -m2 -M2 -v snps -i "INFO/IMP!=1 || INFO/DR2>=${dr2}" -Oz \
    -o "${output_dir}/${cohort}_chr${chromosome}.vcf.gz" "${cohort}.chr${chromosome}.beagle.vcf.gz"
tabix -f -p vcf "${output_dir}/${cohort}_chr${chromosome}.vcf.gz"
bcftools norm --check-ref e --fasta-ref "$reference_fasta" -Ou "${output_dir}/${cohort}_chr${chromosome}.vcf.gz" >/dev/null
retained_variants=$(bcftools index --nrecords "${output_dir}/${cohort}_chr${chromosome}.vcf.gz")
imputed_variants=$(bcftools query -f '%INFO/IMP\n' "${output_dir}/${cohort}_chr${chromosome}.vcf.gz" | awk '$1==1{n++} END{print n+0}')
plink2 --vcf "${output_dir}/${cohort}_chr${chromosome}.vcf.gz" dosage=DS --double-id \
    --set-missing-var-ids '@:#:$r:$a' --make-pgen --threads "$threads" --out "${output_dir}/${cohort}_chr${chromosome}" >> "$log_file" 2>&1

vcf_checksum=$(sha256sum "${output_dir}/${cohort}_chr${chromosome}.vcf.gz" | awk '{print $1}')
index_checksum=$(sha256sum "${output_dir}/${cohort}_chr${chromosome}.vcf.gz.tbi" | awk '{print $1}')
printf 'cohort\tchromosome\tvcf\tindex\tvcf_sha256\tindex_sha256\tbuild\tdr2_threshold\tstatus\n' > "${cohort}.chr${chromosome}.imputation_manifest.tsv"
printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\tPASS\n' "$cohort" "$chromosome" "${cohort}_chr${chromosome}.vcf.gz" "${cohort}_chr${chromosome}.vcf.gz.tbi" "$vcf_checksum" "$index_checksum" "$genome_build" "$dr2" >> "${cohort}.chr${chromosome}.imputation_manifest.tsv"
printf 'cohort\tchromosome\tinput_variants\tretained_variants\timputed_variants\tdr2_threshold\tsample_order\tunique_variant_keys\tdosage_range\treference_alleles\tstatus\n' > "${cohort}.chr${chromosome}.imputation_qc.tsv"
printf '%s\t%s\t%s\t%s\t%s\t%s\tPASS\tPASS\tPASS\tPASS\tPASS\n' "$cohort" "$chromosome" "$input_variants" "$retained_variants" "$imputed_variants" "$dr2" >> "${cohort}.chr${chromosome}.imputation_qc.tsv"
