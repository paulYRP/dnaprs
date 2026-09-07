#!/usr/bin/env bash
set -euo pipefail

repo_root=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
test_root=$(mktemp -d)
cd "$test_root"
printf '%s\n' '##fileformat=VCFv4.2' '##contig=<ID=1,length=1000>' \
    '##INFO=<ID=IMP,Number=0,Type=Flag,Description="Imputed">' \
    '##INFO=<ID=DR2,Number=1,Type=Float,Description="Dosage quality">' \
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">' \
    '##FORMAT=<ID=DS,Number=1,Type=Float,Description="ALT dosage">' > header.vcf
printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n' >> header.vcf
cp header.vcf good.vcf
printf '1\t100\trs100\tA\tG\t.\tPASS\tDR2=1\tGT:DS\t0/0:0\t0/1:1.25\n' >> good.vcf
bash "$repo_root/bin/validate_imputed_dosages.sh" good.vcf TEST 1
for dosage in . -0.1 2.1 nan inf 1,0 malformed; do
    cp header.vcf invalid.vcf
    printf '1\t100\trs100\tA\tG\t.\tPASS\tDR2=1\tGT:DS\t0/0:0\t0/1:%s\n' "$dosage" >> invalid.vcf
    if bash "$repo_root/bin/validate_imputed_dosages.sh" invalid.vcf TEST 1 > check.log 2>&1; then
        printf 'Invalid retained DS=%s was accepted.\n' "$dosage" >&2
        exit 1
    fi
    if [[ "$dosage" == . ]]; then
        grep -q 'participant S2, variant 1:100 (rs100): DS=.' check.log
    fi
done
sed '/^##FORMAT=<ID=DS,/d; s/GT:DS/GT/g; s/0\/0:0/0\/0/g; s/0\/1:1.25/0\/1/g' good.vcf > no_ds.vcf
if bash "$repo_root/bin/validate_imputed_dosages.sh" no_ds.vcf TEST 1 > check.log 2>&1; then exit 1; fi
grep -q 'requires a FORMAT/DS header' check.log
cp good.vcf low_quality.vcf
printf '1\t200\trs200\tC\tT\t.\tPASS\tIMP;DR2=0.2\tGT:DS\t0/0:.\t0/1:.\n' >> low_quality.vcf
bcftools view -m2 -M2 -v snps -i 'INFO/IMP!=1 || INFO/DR2>=0.8' -Ov -o retained.vcf low_quality.vcf
bash "$repo_root/bin/validate_imputed_dosages.sh" retained.vcf TEST 1
printf 'Retained dosage completeness, validity and post-filter validation passed.\n'
