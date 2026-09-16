#!/usr/bin/env bash
set -euo pipefail
repo_root=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
task_dir=$(mktemp -d)
cd "$task_dir"
plink2 --vcf "$repo_root/tests/data/target/relatedness_target.vcf" --make-pgen --out source >/dev/null
cp source.pvar original.pvar
awk 'BEGIN {FS=OFS="\t"} /^#/ {print; next}
    {row++; if(row==1) {$1=0; $2=0; $5=$4}
    if(row==2 || row==3) $3="duplicate"
    if(row==4) $3="."
    print}' original.pvar > source.pvar
sha256sum source.pgen source.pvar source.psam > before.sha256
bash "$repo_root/bin/genotype_eda.sh" TEST pgen source '' DS 1 raw target '' ''
sha256sum -c before.sha256
awk -F '\t' 'NR>1 && $2==1 {if($8!="excluded") exit 1; found=1} END {exit !found}' TEST.diagnostic_markers.tsv
awk -F '\t' 'NR>1 && $2=="heterozygosity" {if($3!="PASS" && $3!="REVIEW") exit 1; found=1} END {exit !found}' TEST.genotype_eda_checks.tsv
awk -F '\t' 'NR>1 && $2=="reported_sex" {if($3!="NOT_RUN") exit 1; found=1} END {exit !found}' TEST.genotype_eda_checks.tsv
if grep -q 'Duplicate allele code' TEST.genotype_eda.log; then
    echo 'A downstream diagnostic still read the invalid raw allele record.' >&2
    exit 1
fi
real_plink=$(command -v plink2)
mkdir failing_bin
cat > failing_bin/plink2 <<'WRAPPER'
#!/usr/bin/env bash
for arg in "$@"; do
    if [[ "$arg" == --het || "$arg" == --pca ]]; then exit 9; fi
done
exec "$REAL_PLINK" "$@"
WRAPPER
chmod +x failing_bin/plink2
REAL_PLINK="$real_plink" PATH="$task_dir/failing_bin:$PATH" \
    bash "$repo_root/bin/genotype_eda.sh" FAILED pgen source '' DS 1 raw target '' ''
awk -F '\t' 'NR>1 && ($2=="heterozygosity" || $2=="internal_pca") {if($3!="FAIL") bad=1; found++} END {exit bad || found!=2}' FAILED.genotype_eda_checks.tsv
awk -F '\t' 'NR==1 {for(i=1;i<=NF;i++) ix[$i]=i; next}
    {if($(ix["status"])!="FAIL" || $(ix["fail_items"])<2 || $(ix["not_run_items"])<1) exit 1}' FAILED.genotype_eda_summary.tsv
printf 'Diagnostic subsets preserve source records; attempted failures and skipped checks are distinct.\n'

# Invalid X records must not reach sex-check frequencies.
awk 'BEGIN {FS=OFS="\t"} /^#/ {print; next}
    {row++; if(row<=3) {$1="X"; $2=4000000+row*1000000}
    if(row==1) $5=$4; print}' original.pvar > source.pvar
awk 'BEGIN {FS=OFS="\t"} NR==1 {
    for(i=1;i<=NF;i++) if($i=="SEX") sex=i
    if(!sex) {sex=NF+1; $sex="SEX"}; print; next}
    {$sex=(NR<=4 ? 1 : 2); print}' source.psam > sex.psam
mv sex.psam source.psam
bash "$repo_root/bin/genotype_eda.sh" XTEST pgen source '' DS 1 raw target '' ''
awk -F '\t' 'NR>1 && $2=="x_frequencies" {if($3!="PASS") bad=1; found=1}
    END {exit bad || !found}' XTEST.genotype_eda_checks.tsv
awk -F '\t' 'NR>1 && $2=="reported_sex" {if($3!="PASS" && $3!="REVIEW") bad=1; found=1}
    END {exit bad || !found}' XTEST.genotype_eda_checks.tsv

# A failed prerequisite is FAIL; dependent diagnostics are NOT_RUN.
cat > failing_bin/plink2 <<'WRAPPER'
#!/usr/bin/env bash
for arg in "$@"; do
    if [[ "$arg" == --indep-pairwise ]]; then exit 9; fi
done
exec "$REAL_PLINK" "$@"
WRAPPER
REAL_PLINK="$real_plink" PATH="$task_dir/failing_bin:$PATH" \
    bash "$repo_root/bin/genotype_eda.sh" PRUNING pgen source '' DS 1 raw target '' ''
awk -F '\t' 'NR>1 && $2=="diagnostic_pruning" {if($3!="FAIL") bad=1; found++}
    NR>1 && ($2=="heterozygosity" || $2=="relatedness" || $2=="internal_pca") {if($3!="NOT_RUN") bad=1; found++}
    END {exit bad || found!=4}' PRUNING.genotype_eda_checks.tsv
printf 'Valid X subsets and failed diagnostic prerequisites passed.\n'
