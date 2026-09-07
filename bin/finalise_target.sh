#!/usr/bin/env bash
set -euo pipefail

plink2() {
    command plink2 --memory "${PLINK_MEMORY_MB:-1024}" "$@"
}
cohort="$1"
imported="$2"
input_stage="$3"
assay_manifest="$4"
reference_fasta="$5"
threads="$6"

mkdir -p "$cohort"
for extension in pgen pvar psam; do
    cp "$imported/$cohort.$extension" "$cohort/$cohort.$extension"
done
if [[ "$input_stage" == "raw" ]]; then
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
    export_prefix="$cohort/${cohort}_resolved"
    export_args=(--export vcf bgz id-paste=iid)
    if [[ -z "$assay_manifest" ]]; then
        export_args+=(vcf-dosage=DS)
        awk -F '\t' 'BEGIN { OFS="\t" }
            NR == 1 { for (i=1; i<=NF; i++) column[$i]=i; next }
            $(column["decision"]) ~ /^RETAINED_/ {
                print $(column["final_id"]), $(column["final_ref"]), $(column["final_alt"])
            }
        ' "$cohort.marker_decisions.tsv" > reference_alleles.tsv
        # Complement both allele labels together without changing their PGEN indices.
        # PLINK then changes REF order and recodes genotypes and dosages together.
        awk -F '\t' 'BEGIN { OFS="\t"; complement["A"]="T"; complement["T"]="A"; complement["C"]="G"; complement["G"]="C" }
            NR == FNR { ref[$1]=$2; alt[$1]=$3; next }
            /^#/ { print; next }
            {
                direct=($4 == ref[$3] && $5 == alt[$3]) || ($5 == ref[$3] && $4 == alt[$3])
                reverse=(complement[$4] == ref[$3] && complement[$5] == alt[$3]) || (complement[$5] == ref[$3] && complement[$4] == alt[$3])
                if (!direct && !reverse) {
                    printf "Marker %s has no compatible resolved allele orientation.\n", $3 > "/dev/stderr"
                    exit 5
                }
                if (!direct) { $4=complement[$4]; $5=complement[$5] }
                print
            }
        ' reference_alleles.tsv "$export_prefix.pvar" > strand.pvar
        mv strand.pvar "$export_prefix.pvar"
        plink2 --pfile "$export_prefix" --ref-allele force reference_alleles.tsv 2 1 \
            --make-pgen --threads "$threads" --out "$cohort/${cohort}_oriented"
        export_prefix="$cohort/${cohort}_oriented"
    fi
    plink2 --pfile "$export_prefix" "${export_args[@]}" \
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
    plink2 --vcf "$cohort/${cohort}_forward.checked.vcf.gz" dosage=DS --double-id \
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
    awk -F '\t' -v cohort="$cohort" '
        NR == FNR {
            if (FNR == 1) { for (i=1; i<=NF; i++) column[$i]=i; next }
            if ($(column["decision"]) ~ /^RETAINED_/) {
                id=$(column["final_id"])
                expected[id]=$(column["final_chr"]) FS $(column["final_pos"]) FS $(column["final_ref"]) FS $(column["final_alt"])
                retained++
            }
            next
        }
        /^#/ { next }
        {
            observed++
            chromosome=$1; sub(/^chr/, "", chromosome)
            if (!($3 in expected) || expected[$3] != chromosome FS $2 FS $4 FS $5) {
                printf "Cohort %s marker %s disagrees with the retained GRCh37 coordinate or REF/ALT pair.\n", cohort, $3 > "/dev/stderr"
                failed=1
            }
        }
        END {
            if (observed != retained) {
                printf "Cohort %s has %d final markers; expected %d retained decisions.\n", cohort, observed, retained > "/dev/stderr"
                failed=1
            }
            exit failed
        }
    ' "$cohort.marker_decisions.tsv" "$cohort/$cohort.pvar"
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
