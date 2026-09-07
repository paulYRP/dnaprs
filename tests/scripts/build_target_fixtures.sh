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
