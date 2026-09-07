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
decisions="${7:-$cohort.marker_decisions.tsv}"
decision_input="$cohort.marker_decisions.input.tsv"
cp "$decisions" "$decision_input"

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
    ' "$decision_input" | sort > expected_final_ids.txt
    awk '!/^#/ { print $3 }' "$cohort/${cohort}_resolved.pvar" | sort > observed_final_ids.txt
    if ! cmp -s expected_final_ids.txt observed_final_ids.txt; then
        echo "Resolved PVAR identifiers for cohort '$cohort' do not match the retained marker decisions." >&2
        exit 5
    fi

    cp "$cohort/${cohort}_resolved.psam" resolved.psam
    export_prefix="$cohort/${cohort}_resolved"
    if [[ -z "$assay_manifest" ]]; then
        awk -F '\t' 'BEGIN { OFS="\t" }
            NR == 1 { for (i=1; i<=NF; i++) column[$i]=i; next }
            $(column["decision"]) ~ /^RETAINED_/ {
                print $(column["final_id"]), $(column["assay_ref_a"]), $(column["assay_ref_b"])
            }
        ' "$decision_input" > reference_alleles.tsv
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
    # A task-local marker VCF tracks where the original REF allele moves. The 0/0
    # index is not a participant genotype. Keep it outside the published PGEN directory.
    awk 'BEGIN { FS=OFS="\t"; print "##fileformat=VCFv4.2" }
        !/^#/ { chromosome[$1]=1 }
        END {
            for (chr in chromosome) print "##contig=<ID=" chr ">"
            print "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Original REF allele index\">"
            print "#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","ALLELE_INDEX"
        }
    ' "$export_prefix.pvar" > "$cohort.allele_orientation.vcf"
    awk 'BEGIN { FS=OFS="\t" } !/^#/ { print $1,$2,$3,$4,$5,".","PASS",".","GT","0/0" }' \
        "$export_prefix.pvar" >> "$cohort.allele_orientation.vcf"
    if [[ -n "$assay_manifest" ]]; then
        bcftools +fixref "$cohort.allele_orientation.vcf" -Oz -o "$cohort.forward.vcf.gz" -- \
            -f "$reference_fasta" -m top
    else
        bcftools view -Oz -o "$cohort.forward.vcf.gz" "$cohort.allele_orientation.vcf"
    fi
    if ! bcftools norm -f "$reference_fasta" -c e -Oz \
        -o "$cohort.forward.checked.vcf.gz" "$cohort.forward.vcf.gz" 2> reference_check.log; then
        cat reference_check.log >&2
        awk -F '\t' '
            NR == FNR {
                if (FNR == 1) { for(i=1;i<=NF;i++) column[$i]=i; next }
                if ($(column["decision"]) ~ /^RETAINED_/) {
                    key=$(column["final_chr"]) ":" $(column["final_pos"])
                    context[key]="Source marker " $(column["source_id"]) ", rsID " $(column["final_id"]) ", coordinate " key ", assay " $(column["assay_ref_a"]) "/" $(column["assay_ref_b"])
                }
                next
            }
            { count=split($0, token, /[[:space:]]+/); for(i=1;i<=count;i++) if(token[i] in context) print context[token[i]] }
        ' "$decision_input" reference_check.log >&2
        echo "Cohort '$cohort' has an assay/reference conflict. Review its source alleles, genome build and FASTA; candidate rsID matching does not validate genomic REF/ALT." >&2
        exit 5
    fi
    cat reference_check.log >&2
    tabix -f -p vcf "$cohort.forward.checked.vcf.gz"
    if [[ -n "$assay_manifest" ]]; then
        # Use the reference-checked TOP transformation, including ambiguous SNPs.
        # Apply its strand change to allele labels without changing PGEN indices.
        bcftools query -f '%ID\t%REF\t%ALT[\t%GT]\n' \
            "$cohort.forward.checked.vcf.gz" > reference_alleles.tsv
        awk -F '\t' 'BEGIN { OFS="\t"; complement["A"]="T"; complement["T"]="A"; complement["C"]="G"; complement["G"]="C" }
            NR == FNR {
                if ($4 != "0/0" && $4 != "1/1") {
                    printf "Marker %s has no completed TOP-to-forward transformation.\n", $1 > "/dev/stderr"
                    exit 5
                }
                ref[$1]=$2; alt[$1]=$3; original_ref[$1]=($4 == "0/0" ? $2 : $3); next
            }
            /^#/ { print; next }
            {
                if ($4 != original_ref[$3]) { $4=complement[$4]; $5=complement[$5] }
                if ($4 != original_ref[$3]) {
                    printf "Marker %s has an invalid original REF transformation.\n", $3 > "/dev/stderr"
                    exit 5
                }
                if (!(($4 == ref[$3] && $5 == alt[$3]) || ($5 == ref[$3] && $4 == alt[$3]))) {
                    printf "Marker %s disagrees with its reference-checked TOP transformation.\n", $3 > "/dev/stderr"
                    exit 5
                }
                print
            }
        ' reference_alleles.tsv "$export_prefix.pvar" > strand.pvar
        mv strand.pvar "$export_prefix.pvar"
        # PLINK recodes hard calls and dosages together when changing REF order.
        plink2 --pfile "$export_prefix" --ref-allele force reference_alleles.tsv 2 1 \
            --make-pgen --threads "$threads" --out "$cohort/${cohort}_oriented"
        export_prefix="$cohort/${cohort}_oriented"
    fi
    # Keep native genotypes for QC. A VCF reimport can lose raw dosages or sex metadata.
    if ! cmp -s resolved.psam "$export_prefix.psam"; then
        echo "Sample metadata changed during GRCh37 orientation for cohort '$cohort'." >&2
        exit 5
    fi
    for extension in pgen pvar psam; do
        cp "$export_prefix.$extension" "$cohort/${cohort}.$extension"
    done
    duplicate_count=$(awk '!/^#/ {count[$3]++} END {n=0; for (id in count) if (count[id] > 1) n++; print n}' "$cohort/${cohort}.pvar")
    [[ "$duplicate_count" == "0" ]] || { echo "Resolved target for cohort '$cohort' contains repeated marker identifiers." >&2; exit 5; }
    invalid_alleles=$(awk '!/^#/ && ($4 !~ /^[ACGT]$/ || $5 !~ /^[ACGT]$/ || $4 == $5) { count++ } END { print count + 0 }' "$cohort/${cohort}.pvar")
    [[ "$invalid_alleles" == "0" ]] || {
        echo "Resolved target for cohort '$cohort' contains $invalid_alleles invalid REF and ALT allele pairs." >&2
        exit 5
    }
    # Write a new decision table; never edit a staged upstream symlink.
    awk -F '\t' -v cohort="$cohort" 'BEGIN { OFS="\t" }
        NR == FNR {
            if (!/^#/) { chr[$3]=$1; pos[$3]=$2; ref[$3]=$4; alt[$3]=$5 }
            next
        }
        FNR == 1 { for(i=1;i<=NF;i++) column[$i]=i; print; next }
        {
            if ($(column["decision"]) ~ /^RETAINED_/) {
                id=$(column["final_id"])
                if (!(id in chr) || chr[id] != $(column["final_chr"]) || pos[id] != $(column["final_pos"])) {
                    printf "Cohort %s marker %s disagrees with its resolved coordinate.\n", cohort, id > "/dev/stderr"
                    failed=1
                }
                $(column["final_chr"])=chr[id]; $(column["final_pos"])=pos[id]
                $(column["final_ref"])=ref[id]; $(column["final_alt"])=alt[id]
            }
            print
        }
        END { exit failed }
    ' "$cohort/$cohort.pvar" "$decision_input" > "$cohort.marker_decisions.new.tsv"
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
    ' "$cohort.marker_decisions.new.tsv" "$cohort/$cohort.pvar"
else
    cp "$decision_input" "$cohort.marker_decisions.new.tsv"
fi
mv "$cohort.marker_decisions.new.tsv" "$cohort.marker_decisions.tsv"

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
