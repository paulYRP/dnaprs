#!/usr/bin/env bash
set -euo pipefail

cohort="$1"
format="$2"
genotype="$3"
sample="$4"
dosage="$5"
threads="$6"
input_stage="${7:-qc_completed}"
role="${8:-target}"
assay_manifest="${9:-}"
adapter_script="${10:-}"

work_prefix="${cohort}.genotype_eda_work"
raw_prefix="${work_prefix}/${cohort}_raw"
analysis_prefix="${work_prefix}/${cohort}_analysis"
stage_log="${cohort}.genotype_eda.log"
mkdir -p "$work_prefix"
: > "$stage_log"

log_message() {
    printf '[%s] %s\n' "$(date -u +'%Y-%m-%dT%H:%M:%SZ')" "$*" >> "$stage_log"
}

run_plink() {
    local label="$1"
    shift
    log_message "$label: plink2 $*"
    plink2 "$@" >> "$stage_log" 2>&1
}

run_plink1() {
    local label="$1"
    shift
    log_message "$label: plink $*"
    plink "$@" >> "$stage_log" 2>&1
}

import_target() {
    local source_path="$1"
    local output_prefix="$2"
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
                printf 'GenomeStudio input requires an assay manifest and target adapter.\n' >&2
                exit 2
            fi
            local converted_prefix="${output_prefix}_genomestudio"
            perl "$adapter_script" finalreport "$source_path" "$assay_manifest" "$converted_prefix"
            args=(--pedmap "$converted_prefix")
            ;;
        vcf)
            args=(--vcf "$source_path")
            if [[ -n "$dosage" ]]; then args+=(dosage="$dosage"); fi
            args+=(--double-id)
            ;;
        bgen)
            args=(--bgen "$source_path" ref-first)
            if [[ -n "$sample" ]]; then args+=(--sample "$sample"); fi
            ;;
        *)
            printf 'Unsupported target format: %s\n' "$format" >&2
            exit 2
            ;;
    esac

    run_plink "Import raw target" "${args[@]}" --make-pgen --threads "$threads" --out "$output_prefix"
}

pattern_input=false
if [[ "$genotype" == *'{chr}'* || "$genotype" == *'{CHR}'* || "$genotype" == *'{chromosome}'* || "$genotype" == *'#'* ]]; then
    pattern_input=true
fi

if [[ "$pattern_input" == true ]]; then
    prefixes=()
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
            chromosome_prefix="${work_prefix}/${cohort}_source_chr${chromosome}"
            import_target "$source_path" "$chromosome_prefix"
            prefixes+=("$chromosome_prefix")
        fi
    done
    if [[ "${#prefixes[@]}" -eq 0 ]]; then
        printf 'No target files matched the chromosome placeholder for %s.\n' "$cohort" >&2
        exit 3
    fi
    printf '%s\n' "${prefixes[@]}" > "${work_prefix}/${cohort}_merge.txt"
    if [[ "${#prefixes[@]}" -eq 1 ]]; then
        cp "${prefixes[0]}.pgen" "${raw_prefix}.pgen"
        cp "${prefixes[0]}.pvar" "${raw_prefix}.pvar"
        cp "${prefixes[0]}.psam" "${raw_prefix}.psam"
    else
        run_plink "Merge raw chromosome inputs" \
            --pmerge-list "${work_prefix}/${cohort}_merge.txt" pfile \
            --make-pgen --threads "$threads" --out "$raw_prefix"
    fi
else
    import_target "$genotype" "$raw_prefix"
fi

sample_count=$(awk 'BEGIN{n=0} !/^#/ && NF>0 {n++} END{print n}' "${raw_prefix}.psam")
variant_count=$(awk 'BEGIN{n=0} !/^#/ && NF>0 {n++} END{print n}' "${raw_prefix}.pvar")
if [[ "$sample_count" -eq 0 || "$variant_count" -eq 0 ]]; then
    printf 'Target %s contains no participants or variants after format import.\n' "$cohort" >&2
    exit 4
fi

read -r recorded_sex_count phenotype_count duplicated_participants < <(
    awk '
        BEGIN { FS=OFS="\t" }
        NR == 1 {
            for (column = 1; column <= NF; column++) {
                name = $column
                sub(/^#/, "", name)
                column_index[name] = column
                if (name !~ /^(FID|IID|SID|PAT|MAT|SEX)$/) phenotype[column] = 1
            }
            next
        }
        {
            fid = (column_index["FID"] ? $(column_index["FID"]) : "0")
            iid = $(column_index["IID"])
            participant[fid SUBSEP iid]++
            if (column_index["SEX"]) {
                sex = toupper($(column_index["SEX"]))
                if (sex != "" && sex != "." && sex != "NA" && sex != "0") recorded++
            }
            has_phenotype = 0
            for (column in phenotype) {
                value = toupper($(column))
                if (value != "" && value != "." && value != "NA" && value != "-9") has_phenotype = 1
            }
            phenotype_present += has_phenotype
        }
        END {
            duplicates = 0
            for (key in participant) if (participant[key] > 1) duplicates++
            print recorded + 0, phenotype_present + 0, duplicates + 0
        }
    ' "${raw_prefix}.psam"
)

awk -v cohort="$cohort" '
    BEGIN { FS=OFS="\t"; print "cohort", "chromosome", "category", "variants" }
    !/^#/ {
        chromosome=$1
        sub(/^chr/, "", chromosome)
        upper=toupper(chromosome)
        if (upper ~ /^[0-9]+$/ && upper + 0 >= 1 && upper + 0 <= 22) category="autosome"
        else if (upper == "X" || upper == "23") category="sex_chromosome_x"
        else if (upper == "Y" || upper == "24") category="sex_chromosome_y"
        else if (upper == "XY" || upper == "25") category="pseudoautosomal"
        else if (upper == "M" || upper == "MT" || upper == "26") category="mitochondrial"
        else category="unplaced_or_nonstandard"
        count[chromosome SUBSEP category]++
    }
    END {
        for (key in count) {
            split(key, value, SUBSEP)
            print cohort, value[1], value[2], count[key]
        }
    }
' "${raw_prefix}.pvar" | { IFS= read -r header; printf '%s\n' "$header"; sort -t $'\t' -k3,3 -k2,2V; } \
    > "${cohort}.chromosome_counts.tsv"

awk -v cohort="$cohort" '
    BEGIN { FS=OFS="\t"; print "cohort", "chromosome", "bin_start", "bin_end", "variants" }
    !/^#/ {
        chromosome=$1
        sub(/^chr/, "", chromosome)
        if (chromosome ~ /^[0-9]+$/ && chromosome + 0 >= 1 && chromosome + 0 <= 22 && $2 ~ /^[0-9]+$/ && $2 > 0) {
            start=int(($2 - 1) / 1000000) * 1000000 + 1
            count[chromosome SUBSEP start]++
        }
    }
    END {
        for (key in count) {
            split(key, value, SUBSEP)
            print cohort, value[1], value[2], value[2] + 999999, count[key]
        }
    }
' "${raw_prefix}.pvar" | { IFS= read -r header; printf '%s\n' "$header"; sort -t $'\t' -k2,2n -k3,3n; } \
    > "${cohort}.marker_density.tsv"

awk -v cohort="$cohort" '
    BEGIN { FS=OFS="\t"; print "cohort", "identifier_class", "variants", "percent" }
    !/^#/ {
        identifier=$3
        lower=tolower(identifier)
        if (identifier == "" || identifier == "." || identifier == "0") class="missing"
        else if (lower ~ /^rs[0-9]+$/) class="rsid"
        else if (identifier ~ /^(chr)?([0-9]+|X|Y|XY|M|MT)[:_-][0-9]+([:_-].*)?$/) class="coordinate"
        else if (identifier ~ /^[A-Za-z][A-Za-z0-9]*[-_.][A-Za-z0-9_.-]+$/) class="array_probe"
        else class="other"
        count[class]++
        total++
    }
    END {
        order[1]="rsid"; order[2]="array_probe"; order[3]="coordinate"; order[4]="missing"; order[5]="other"
        for (i=1; i<=5; i++) {
            class=order[i]
            print cohort, class, count[class] + 0, (total ? 100 * (count[class] + 0) / total : 0)
        }
    }
' "${raw_prefix}.pvar" > "${cohort}.identifier_classes.tsv"

awk -v cohort="$cohort" '
    BEGIN { FS=OFS="\t"; print "cohort", "allele_state", "variants", "percent" }
    !/^#/ {
        ref=toupper($4); alt=toupper($5)
        ref_missing=(ref == "" || ref == "." || ref == "0")
        alt_missing=(alt == "" || alt == "." || alt == "0")
        if (ref_missing && alt_missing) state="both_alleles_missing"
        else if (ref_missing || alt_missing) state="one_allele_missing"
        else if (index(alt, ",") > 0) state="multiallelic"
        else if (ref ~ /^[ACGT]$/ && alt ~ /^[ACGT]$/) state="acgt_snp"
        else if (length(ref) != 1 || length(alt) != 1 || ref == "-" || alt == "-" || ref ~ /^</ || alt ~ /^</) state="insertion_deletion"
        else state="other"
        count[state]++
        total++
    }
    END {
        order[1]="acgt_snp"; order[2]="one_allele_missing"; order[3]="both_alleles_missing"
        order[4]="insertion_deletion"; order[5]="multiallelic"; order[6]="other"
        for (i=1; i<=6; i++) {
            state=order[i]
            print cohort, state, count[state] + 0, (total ? 100 * (count[state] + 0) / total : 0)
        }
    }
' "${raw_prefix}.pvar" > "${cohort}.allele_states.tsv"

{
    printf 'cohort\tidentifier\toccurrences\n'
    awk 'BEGIN{FS="\t"} !/^#/ && $3 != "" && $3 != "." && $3 != "0" {print $3}' "${raw_prefix}.pvar" |
        sort | uniq -c | awk -v cohort="$cohort" 'BEGIN{OFS="\t"} $1 > 1 {print cohort, $2, $1}'
} > "${cohort}.duplicate_identifiers.tsv"

duplicate_identifier_groups=$(awk 'NR > 1 {n++} END{print n+0}' "${cohort}.duplicate_identifiers.tsv")
duplicate_identifier_variants=$(awk 'NR > 1 {n += $3} END{print n+0}' "${cohort}.duplicate_identifiers.tsv")

autosomal_variants=$(awk -F '\t' 'NR > 1 && $3 == "autosome" {n += $4} END{print n+0}' "${cohort}.chromosome_counts.tsv")
x_variants=$(awk -F '\t' 'NR > 1 && $3 == "sex_chromosome_x" {n += $4} END{print n+0}' "${cohort}.chromosome_counts.tsv")
y_variants=$(awk -F '\t' 'NR > 1 && $3 == "sex_chromosome_y" {n += $4} END{print n+0}' "${cohort}.chromosome_counts.tsv")
mitochondrial_variants=$(awk -F '\t' 'NR > 1 && $3 == "mitochondrial" {n += $4} END{print n+0}' "${cohort}.chromosome_counts.tsv")
unplaced_variants=$(awk -F '\t' 'NR > 1 && $3 == "unplaced_or_nonstandard" {n += $4} END{print n+0}' "${cohort}.chromosome_counts.tsv")
chromosome_count=$(awk -F '\t' 'NR > 1 {seen[$2]=1} END{for (value in seen) n++; print n+0}' "${cohort}.chromosome_counts.tsv")

run_plink "Allele frequency and genotype missingness" \
    --pfile "$raw_prefix" --freq --missing --threads "$threads" --out "$analysis_prefix"
cp "${analysis_prefix}.afreq" "${cohort}.allele_frequency.tsv"
cp "${analysis_prefix}.smiss" "${cohort}.sample_missingness.tsv"
cp "${analysis_prefix}.vmiss" "${cohort}.variant_missingness.tsv"

awk -v cohort="$cohort" '
    BEGIN { FS=OFS="\t"; print "cohort", "missingness_bin", "variants", "percent" }
    NR == 1 { for (i=1; i<=NF; i++) { name=$i; sub(/^#/, "", name); column_index[name]=i }; next }
    {
        value=$(column_index["F_MISS"]) + 0
        if (value == 0) bin="0"
        else if (value <= .01) bin="(0,0.01]"
        else if (value <= .02) bin="(0.01,0.02]"
        else if (value <= .05) bin="(0.02,0.05]"
        else if (value <= .10) bin="(0.05,0.10]"
        else bin="(0.10,1]"
        count[bin]++
        total++
    }
    END {
        order[1]="0"; order[2]="(0,0.01]"; order[3]="(0.01,0.02]"; order[4]="(0.02,0.05]"; order[5]="(0.05,0.10]"; order[6]="(0.10,1]"
        for (i=1; i<=6; i++) {
            bin=order[i]
            print cohort, bin, count[bin] + 0, (total ? 100 * (count[bin] + 0) / total : 0)
        }
    }
' "${analysis_prefix}.vmiss" > "${cohort}.variant_missingness_bins.tsv"

awk -v cohort="$cohort" '
    BEGIN { FS=OFS="\t"; print "cohort", "maf_bin", "variants", "percent" }
    NR == 1 { for (i=1; i<=NF; i++) { name=$i; sub(/^#/, "", name); column_index[name]=i }; next }
    {
        split($(column_index["ALT_FREQS"]), frequency, ",")
        alt=frequency[1]
        if (alt == "" || alt == "." || tolower(alt) == "nan") next
        maf=(alt <= .5 ? alt : 1-alt)
        if (maf == 0) bin="monomorphic"
        else if (maf <= .01) bin="(0,0.01]"
        else if (maf <= .05) bin="(0.01,0.05]"
        else if (maf <= .10) bin="(0.05,0.10]"
        else if (maf <= .20) bin="(0.10,0.20]"
        else bin="(0.20,0.50]"
        count[bin]++
        total++
    }
    END {
        order[1]="monomorphic"; order[2]="(0,0.01]"; order[3]="(0.01,0.05]"; order[4]="(0.05,0.10]"; order[5]="(0.10,0.20]"; order[6]="(0.20,0.50]"
        for (i=1; i<=6; i++) {
            bin=order[i]
            print cohort, bin, count[bin] + 0, (total ? 100 * (count[bin] + 0) / total : 0)
        }
    }
' "${analysis_prefix}.afreq" > "${cohort}.allele_frequency_bins.tsv"

pruned_count=0
if [[ "$sample_count" -ge 2 && "$autosomal_variants" -ge 2 ]]; then
    bad_ld=()
    if [[ "$sample_count" -lt 50 ]]; then bad_ld=(--bad-ld); fi
    if run_plink "LD pruning for heterozygosity, relatedness, and internal PCA" \
        --pfile "$raw_prefix" --maf 0.05 --indep-pairwise 200kb 0.2 \
        "${bad_ld[@]}" --threads "$threads" --out "${analysis_prefix}_prune"; then
        if [[ ! -s "${analysis_prefix}_prune.prune.in" ]]; then
            if run_plink "Retain eligible autosomal markers when no LD pairs are present" \
                --pfile "$raw_prefix" --autosome --maf 0.05 --write-snplist \
                --threads "$threads" --out "${analysis_prefix}_prune_unpaired"; then
                cp "${analysis_prefix}_prune_unpaired.snplist" "${analysis_prefix}_prune.prune.in"
            fi
        fi
        if [[ -s "${analysis_prefix}_prune.prune.in" ]]; then
            pruned_count=$(awk 'NF>0 {n++} END{print n+0}' "${analysis_prefix}_prune.prune.in")
        fi
    fi
fi

heterozygosity_status="NOT_RUN"
heterozygosity_reason="Fewer than two LD-pruned autosomal variants were available."
printf 'cohort\tFID\tIID\tobserved_homozygotes\texpected_homozygotes\tobservations\tinbreeding_coefficient\theterozygosity_rate\tmissingness\theterozygosity_z\tstatus\n' \
    > "${cohort}.heterozygosity.tsv"
if [[ "$pruned_count" -ge 2 ]]; then
    if run_plink "Autosomal heterozygosity" \
        --pfile "$raw_prefix" --extract "${analysis_prefix}_prune.prune.in" --het \
        --threads "$threads" --out "${analysis_prefix}_heterozygosity"; then
        awk -v cohort="$cohort" '
            BEGIN { FS=OFS="\t" }
            NR == FNR {
                if (FNR == 1) {
                    for (i=1; i<=NF; i++) { name=$i; sub(/^#/, "", name); missing_index[name]=i }
                } else {
                    fid=(missing_index["FID"] ? $(missing_index["FID"]) : "0")
                    iid=$(missing_index["IID"])
                    missing[fid SUBSEP iid]=$(missing_index["F_MISS"])
                }
                next
            }
            FNR == 1 {
                for (i=1; i<=NF; i++) { name=$i; sub(/^#/, "", name); het_index[name]=i }
                next
            }
            {
                n++
                fid_value[n]=(het_index["FID"] ? $(het_index["FID"]) : "0")
                iid_value[n]=$(het_index["IID"])
                observed[n]=$(het_index["O(HOM)"])
                expected[n]=$(het_index["E(HOM)"])
                observations[n]=$(het_index["OBS_CT"])
                coefficient[n]=$(het_index["F"])
                rate[n]=(observations[n] > 0 ? (observations[n] - observed[n]) / observations[n] : "")
                if (rate[n] != "") { sum += rate[n]; finite++ }
            }
            END {
                print "cohort", "FID", "IID", "observed_homozygotes", "expected_homozygotes", "observations", "inbreeding_coefficient", "heterozygosity_rate", "missingness", "heterozygosity_z", "status"
                mean=(finite ? sum / finite : 0)
                for (i=1; i<=n; i++) if (rate[i] != "") square += (rate[i] - mean)^2
                sd=(finite > 1 ? sqrt(square / (finite - 1)) : 0)
                for (i=1; i<=n; i++) {
                    key=fid_value[i] SUBSEP iid_value[i]
                    z=(rate[i] != "" && sd > 0 ? (rate[i] - mean) / sd : 0)
                    status=(sd > 0 && (z < -3 || z > 3) ? "REVIEW" : "PASS")
                    print cohort, fid_value[i], iid_value[i], observed[i], expected[i], observations[i], coefficient[i], rate[i], missing[key], z, status
                }
            }
        ' "${analysis_prefix}.smiss" "${analysis_prefix}_heterozygosity.het" > "${cohort}.heterozygosity.tsv"
        heterozygosity_status="PASS"
        heterozygosity_reason="Autosomal heterozygosity was calculated from the LD-pruned marker set."
    else
        heterozygosity_reason="PLINK could not calculate autosomal heterozygosity; inspect the stage log."
    fi
fi

printf 'cohort\tFID\tIID\trecorded_sex\tgenetic_sex\tx_inbreeding_coefficient\tstatus\treason\n' > "${cohort}.sex_check.tsv"
sex_status="NOT_RUN"
sex_reason="Recorded sex or X-chromosome variants were unavailable."
if [[ "$x_variants" -gt 0 && "$recorded_sex_count" -gt 0 ]]; then
    if run_plink "Reported sex check" \
        --pfile "$raw_prefix" --read-freq "${analysis_prefix}.afreq" \
        --check-sex max-female-xf=0.2 min-male-xf=0.8 \
        --threads "$threads" --out "${analysis_prefix}_sex"; then
        awk -v cohort="$cohort" '
            BEGIN { FS=OFS="\t" }
            NR == 1 {
                for (i=1; i<=NF; i++) { name=$i; sub(/^#/, "", name); column_index[name]=i }
                print "cohort", "FID", "IID", "recorded_sex", "genetic_sex", "x_inbreeding_coefficient", "status", "reason"
                next
            }
            {
                fid=(column_index["FID"] ? $(column_index["FID"]) : "0")
                recorded=(column_index["PEDSEX"] ? $(column_index["PEDSEX"]) : "NA")
                inferred=(column_index["SNPSEX"] ? $(column_index["SNPSEX"]) : "NA")
                coefficient=(column_index["F"] ? $(column_index["F"]) : (column_index["XF"] ? $(column_index["XF"]) : ""))
                status=(column_index["STATUS"] ? $(column_index["STATUS"]) : "NA")
                normal_status=(status == "OK" ? "PASS" : "REVIEW")
                reason=(normal_status == "PASS" ? "Recorded and genetic sex agree." : "Recorded or genetic sex is missing or discordant.")
                print cohort, fid, $(column_index["IID"]), recorded, inferred, coefficient, normal_status, reason
            }
        ' "${analysis_prefix}_sex.sexcheck" > "${cohort}.sex_check.tsv"
        sex_problem_count=$(awk -F '\t' 'NR > 1 && $7 != "PASS" {n++} END{print n+0}' "${cohort}.sex_check.tsv")
        sex_status=$([[ "$sex_problem_count" -gt 0 ]] && printf 'REVIEW' || printf 'PASS')
        sex_reason="X-chromosome coefficients were compared with recorded sex."
    else
        sex_reason="PLINK could not complete the sex check; inspect the stage log."
    fi
fi

printf 'cohort\tFID1\tIID1\tFID2\tIID2\tvariants\tpi_hat\tz0\tz1\tz2\trelationship_category\treview_status\n' \
    > "${cohort}.relatedness.tsv"
printf 'cohort\tpi_hat_bin\tpairs\tpercent\n' > "${cohort}.relatedness_bins.tsv"
printf 'cohort\tFID\tIID\tPC1\tPC2\tPC3\tPC4\tPC5\tPC6\tPC7\tPC8\tPC9\tPC10\n' \
    > "${cohort}.internal_pca.tsv"
printf 'cohort\tcomponent\teigenvalue\tpercent_of_reported_eigenvalues\n' \
    > "${cohort}.pca_eigenvalues.tsv"

relatedness_status="NOT_RUN"
relatedness_reason="At least two participants and two informative autosomal variants are required."
pca_status="NOT_RUN"
pca_reason="$relatedness_reason"
if [[ "$pruned_count" -ge 2 ]]; then
    if run_plink "Prepare PLINK 1 relatedness input" \
        --pfile "$raw_prefix" --extract "${analysis_prefix}_prune.prune.in" \
        --make-bed --threads "$threads" --out "${analysis_prefix}_relatedness_input" \
        && run_plink1 "Pairwise PLINK 1 identity by descent" \
        --bfile "${analysis_prefix}_relatedness_input" --genome full --threads "$threads" \
        --out "${analysis_prefix}_relatedness"; then
        awk -v cohort="$cohort" -v variants="$pruned_count" '
            BEGIN { OFS="\t" }
            NR == 1 {
                for (i=1; i<=NF; i++) { name=$i; sub(/^#/, "", name); column_index[name]=i }
                print "cohort", "FID1", "IID1", "FID2", "IID2", "variants", "pi_hat", "z0", "z1", "z2", "relationship_category", "review_status"
                next
            }
            {
                pi_hat=$(column_index["PI_HAT"])
                if (tolower(pi_hat) == "nan" || pi_hat == "") category="unresolved"
                else if (pi_hat >= .9) category="duplicate_or_monozygotic"
                else if (pi_hat >= .375) category="first_degree"
                else if (pi_hat >= .1875) category="second_degree"
                else if (pi_hat >= .0884) category="third_degree"
                else category="unrelated"
                status=(pi_hat >= .1875 || category == "unresolved" ? "REVIEW" : "PASS")
                print cohort, $(column_index["FID1"]), $(column_index["IID1"]), $(column_index["FID2"]), $(column_index["IID2"]), variants, pi_hat, $(column_index["Z0"]), $(column_index["Z1"]), $(column_index["Z2"]), category, status
            }
        ' "${analysis_prefix}_relatedness.genome" > "${cohort}.relatedness.tsv"
        awk -v cohort="$cohort" '
            BEGIN { FS=OFS="\t"; print "cohort", "pi_hat_bin", "pairs", "percent" }
            NR > 1 {
                value=$7
                if (tolower(value) == "nan" || value == "") bin="unresolved"
                else if (value < 0) bin="negative"
                else if (value < .0884) bin="[0,third_degree)"
                else if (value < .1875) bin="[third_degree,second_degree)"
                else if (value < .375) bin="[second_degree,first_degree)"
                else if (value < .9) bin="[first_degree,duplicate)"
                else bin="duplicate_or_monozygotic"
                count[bin]++
                total++
            }
            END {
                order[1]="negative"; order[2]="[0,third_degree)"; order[3]="[third_degree,second_degree)"; order[4]="[second_degree,first_degree)"; order[5]="[first_degree,duplicate)"; order[6]="duplicate_or_monozygotic"; order[7]="unresolved"
                for (i=1; i<=7; i++) { bin=order[i]; print cohort, bin, count[bin]+0, (total ? 100*(count[bin]+0)/total : 0) }
            }
        ' "${cohort}.relatedness.tsv" > "${cohort}.relatedness_bins.tsv"
        related_problem_count=$(awk -F '\t' 'NR > 1 && $12 != "PASS" {n++} END{print n+0}' "${cohort}.relatedness.tsv")
        relatedness_status=$([[ "$related_problem_count" -gt 0 ]] && printf 'REVIEW' || printf 'PASS')
        relatedness_reason="PLINK 1 identity by descent was calculated from the LD-pruned exploratory marker set; PI_HAT >= 0.1875 was flagged."
    else
        relatedness_reason="PLINK could not calculate identity by descent; inspect the stage log."
    fi

    pc_count=10
    if [[ "$pc_count" -gt $((sample_count - 1)) ]]; then pc_count=$((sample_count - 1)); fi
    if [[ "$pc_count" -gt "$pruned_count" ]]; then pc_count="$pruned_count"; fi
    if [[ "$pc_count" -ge 2 ]] && run_plink "Internal target PCA" \
        --pfile "$raw_prefix" --extract "${analysis_prefix}_prune.prune.in" \
        --pca "$pc_count" --threads "$threads" --out "${analysis_prefix}_pca"; then
        awk -v cohort="$cohort" '
            BEGIN { FS=OFS="\t" }
            NR == 1 {
                for (i=1; i<=NF; i++) { name=$i; sub(/^#/, "", name); column_index[name]=i }
                print "cohort", "FID", "IID", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10"
                next
            }
            {
                printf "%s\t%s\t%s", cohort, (column_index["FID"] ? $(column_index["FID"]) : "0"), $(column_index["IID"])
                for (component=1; component<=10; component++) {
                    name="PC" component
                    printf "\t%s", (column_index[name] ? $(column_index[name]) : "")
                }
                printf "\n"
            }
        ' "${analysis_prefix}_pca.eigenvec" > "${cohort}.internal_pca.tsv"
        awk -v cohort="$cohort" '
            NF > 0 { n++; value[n]=$1; total += $1 }
            END {
                print "cohort\tcomponent\teigenvalue\tpercent_of_reported_eigenvalues"
                for (i=1; i<=n; i++) print cohort "\tPC" i "\t" value[i] "\t" (total ? 100*value[i]/total : 0)
            }
        ' "${analysis_prefix}_pca.eigenval" > "${cohort}.pca_eigenvalues.tsv"
        pca_status="PASS"
        pca_reason="Internal PCs were calculated from the LD-pruned target marker set."
    else
        pca_reason="PLINK could not calculate at least two internal PCs; inspect marker count and the stage log."
    fi
else
    relatedness_reason="LD pruning retained fewer than two informative variants; relatedness was not estimated."
    pca_reason="LD pruning retained fewer than two informative variants; internal PCA was not estimated."
fi

identifier_status=$([[ "$duplicate_identifier_groups" -gt 0 ]] && printf 'REVIEW' || printf 'PASS')
identifier_reason=$([[ "$duplicate_identifier_groups" -gt 0 ]] && printf '%s duplicated identifier group(s) require review.' "$duplicate_identifier_groups" || printf 'No duplicated non-missing identifiers were found.')

{
    printf 'cohort\tcheck\tstatus\tvalue\treason\n'
    printf '%s\tformat_import\tPASS\t%s participants; %s variants\tThe supplied target was imported without changing the source files.\n' "$cohort" "$sample_count" "$variant_count"
    printf '%s\tvariant_identifiers\t%s\t%s duplicate groups\t%s\n' "$cohort" "$identifier_status" "$duplicate_identifier_groups" "$identifier_reason"
    printf '%s\theterozygosity\t%s\t%s participants\t%s\n' "$cohort" "$heterozygosity_status" "$sample_count" "$heterozygosity_reason"
    printf '%s\treported_sex\t%s\t%s recorded; %s X variants\t%s\n' "$cohort" "$sex_status" "$recorded_sex_count" "$x_variants" "$sex_reason"
    printf '%s\trelatedness\t%s\t%s pruned variants\t%s\n' "$cohort" "$relatedness_status" "$pruned_count" "$relatedness_reason"
    printf '%s\tinternal_pca\t%s\t%s pruned variants\t%s\n' "$cohort" "$pca_status" "$pruned_count" "$pca_reason"
} > "${cohort}.genotype_eda_checks.tsv"

review_count=$(awk -F '\t' 'NR > 1 && $3 == "REVIEW" {n++} END{print n+0}' "${cohort}.genotype_eda_checks.tsv")
overall_status=$([[ "$review_count" -gt 0 ]] && printf 'REVIEW' || printf 'PASS')

printf 'cohort\trole\tinput_stage\tformat\tparticipants\tvariants\tchromosomes\tautosomal_variants\tx_variants\ty_variants\tmitochondrial_variants\tunplaced_or_nonstandard_variants\trecorded_sex_participants\tphenotype_participants\tduplicated_participant_identifiers\tduplicated_variant_identifier_groups\tduplicated_variant_identifier_records\treview_items\tstatus\n' \
    > "${cohort}.genotype_eda_summary.tsv"
printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$cohort" "$role" "$input_stage" "$format" "$sample_count" "$variant_count" "$chromosome_count" \
    "$autosomal_variants" "$x_variants" "$y_variants" "$mitochondrial_variants" "$unplaced_variants" \
    "$recorded_sex_count" "$phenotype_count" "$duplicated_participants" "$duplicate_identifier_groups" \
    "$duplicate_identifier_variants" "$review_count" "$overall_status" >> "${cohort}.genotype_eda_summary.tsv"

log_message "Genotype EDA completed with status ${overall_status}."
