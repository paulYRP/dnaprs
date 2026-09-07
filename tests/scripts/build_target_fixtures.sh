#!/usr/bin/env bash
set -euo pipefail
repo_root=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
test_root=$(mktemp -d)
cd "$test_root"
plink2 --vcf "$repo_root/tests/data/target/test_target.vcf" --double-id --make-pgen \
    --threads 1 --memory 1024 --out test_target
plink2 --pfile test_target --make-bed --threads 1 --memory 1024 --out target_bed
plink2 --pfile test_target --export ped --threads 1 --memory 1024 --out target_ped
plink2 --pfile test_target --export bgen-1.2 bits=16 --threads 1 --memory 1024 --out target_bgen
cp test_target.pgen test_target.pvar test_target.psam "$repo_root/tests/data/target/"
cp target_bed.bed target_bed.bim target_bed.fam target_ped.ped target_ped.map \
    target_bgen.bgen target_bgen.sample "$repo_root/tests/data/target/formats/"
cp target_ped.ped "$repo_root/tests/data/resolver inputs/ukr/ukr.ped"
cp target_ped.map "$repo_root/tests/data/resolver inputs/ukr/ukr.map"
awk 'BEGIN {
    OFS=","
    print "[Header]"
    print "Descriptor File Name,target_bed_manifest.csv"
    print "[Assay]"
    print "IlmnID,Name,IlmnStrand,SNP,GenomeBuild,Chr,MapInfo,RefStrand"
}
{print "probe" NR,$2,"TOP","[" $6 "/" $5 "]",37,$1,$4,"+"}
' target_bed.bim > "$repo_root/tests/data/target/formats/target_bed_manifest.csv"
