#!/usr/bin/env bash
set -euo pipefail

repo_root=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
test_root=$(mktemp -d)
cd "$test_root"
printf '%s\n' \
    '##fileformat=VCFv4.2' \
    '##contig=<ID=1,length=1000>' \
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">' \
    '##FORMAT=<ID=DS,Number=1,Type=Float,Description="ALT dosage">' > dosage.vcf
printf '%b\n' \
    '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2' \
    '1\t100\tdosage1\tG\tA\t.\tPASS\t.\tGT:DS\t./.:1.25\t./.:0.75' \
    '1\t200\tdosage2\tT\tC\t.\tPASS\t.\tGT:DS\t./.:0.25\t./.:1.75' >> dosage.vcf

bash "$repo_root/bin/import_target.sh" DOSAGE vcf dosage.vcf '' '' DS 1 raw '' '' "$repo_root/bin/target_adapter.pl"
bash "$repo_root/bin/resolve_target_markers.sh" DOSAGE DOSAGE.imported DOSAGE.initial_marker_decisions.tsv raw \
    "$repo_root/tests/data/reference/dbsnp_source" "$repo_root/bin/resolve_target_markers.R" 1
bash "$repo_root/bin/finalise_target.sh" DOSAGE DOSAGE.imported raw '' \
    "$repo_root/tests/data/reference/human_g1k_v37.fasta" 1
plink2 --pfile DOSAGE/DOSAGE --export A-transpose --threads 1 --memory 1024 --out observed
awk 'NR==2 { if ($5 != "A" || $7 != 1.25 || $8 != 0.75) exit 1 }
     NR==3 { if ($5 != "C" || $7 != 0.25 || $8 != 1.75) exit 1 }
     END { if (NR != 3) exit 1 }' observed.traw
test ! -s DOSAGE.duplicate_markers.txt
printf 'Marker resolution preserves dosages and corrects reversed reference alleles.\n'
