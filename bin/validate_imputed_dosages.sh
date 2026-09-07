#!/usr/bin/env bash
set -euo pipefail

vcf="$1"
cohort="$2"
chromosome="$3"
samples="${cohort}.chr${chromosome}.dosage_samples.txt"
bcftools query -l "$vcf" > "$samples"
[[ -s "$samples" ]] || { echo "Cohort '$cohort', chromosome $chromosome: retained VCF has no participants." >&2; exit 5; }
if ! bcftools view -h "$vcf" | awk '/^##FORMAT=<ID=DS,/ { found=1 } END { exit !found }'; then
    echo "Cohort '$cohort', chromosome $chromosome: retained VCF requires a FORMAT/DS header." >&2
    exit 5
fi

# Stream the retained records. Drain the input on failure so diagnostics are not
# replaced by an upstream broken-pipe error.
bcftools query -f '%CHROM\t%POS\t%ID[\t%DS]\n' "$vcf" |
    awk -F '\t' -v cohort="$cohort" -v chromosome="$chromosome" '
        FILENAME == ARGV[1] { sample[++samples]=$0; next }
        {
            records++
            if (NF != samples + 3) {
                if (errors++ < 5) printf "Cohort %s, variant %s:%s (%s): expected %d dosages, found %d.\n", cohort, $1, $2, $3, samples, NF-3 > "/dev/stderr"
            }
            for (i=4; i<=NF; i++) {
                value=$i
                if (value !~ /^[+-]?([0-9]+([.][0-9]*)?|[.][0-9]+)([eE][+-]?[0-9]+)?$/ || value+0 < 0 || value+0 > 2) {
                    if (errors++ < 5) printf "Cohort %s, participant %s, variant %s:%s (%s): DS=%s; a finite dosage in [0,2] is required.\n", cohort, sample[i-3], $1, $2, $3, value > "/dev/stderr"
                }
            }
        }
        END {
            if (!records) {
                printf "Cohort %s, chromosome %s: no retained dosage records.\n", cohort, chromosome > "/dev/stderr"
                errors++
            }
            exit (errors > 0)
        }
    ' "$samples" -
