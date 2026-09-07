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
mkdir upstream
mv DOSAGE.marker_decisions.tsv upstream/
ln -s upstream/DOSAGE.marker_decisions.tsv DOSAGE.marker_decisions.tsv
sha256sum upstream/DOSAGE.marker_decisions.tsv DOSAGE.imported/DOSAGE.* > upstream.sha256
bash "$repo_root/bin/finalise_target.sh" DOSAGE DOSAGE.imported raw '' \
    "$repo_root/tests/data/reference/human_g1k_v37.fasta" 1
sha256sum --check --strict upstream.sha256
test ! -L DOSAGE.marker_decisions.tsv
Rscript -e '
    library(data.table)
    decisions <- fread("DOSAGE.marker_decisions.tsv")[startsWith(decision, "RETAINED_")]
    variants <- fread("DOSAGE/DOSAGE.pvar")
    variants <- variants[match(decisions$final_id, ID)]
    stopifnot(identical(decisions$final_ref, variants$REF), identical(decisions$final_alt, variants$ALT))
'
plink2 --pfile DOSAGE/DOSAGE --export A-transpose --threads 1 --memory 1024 --out observed
awk 'NR==2 { if ($5 != "A" || $7 != 1.25 || $8 != 0.75) exit 1 }
     NR==3 { if ($5 != "C" || $7 != 0.25 || $8 != 1.75) exit 1 }
     END { if (NR != 3) exit 1 }' observed.traw
test ! -s DOSAGE.duplicate_markers.txt
printf 'Marker resolution preserves dosages and corrects reversed reference alleles.\n'

mkdir manifest_dosages
cp -r DOSAGE.imported manifest_dosages/
cp DOSAGE.retained_markers.txt DOSAGE.rename_markers.tsv DOSAGE.marker_decisions.tsv manifest_dosages/
(
    cd manifest_dosages
    bash "$repo_root/bin/finalise_target.sh" DOSAGE DOSAGE.imported raw top \
        "$repo_root/tests/data/reference/human_g1k_v37.fasta" 1
    plink2 --pfile DOSAGE/DOSAGE --export Av --threads 1 --memory 1024 --out observed
    cmp ../observed.traw observed.traw
)
printf 'TOP-to-forward correction preserves raw dosages.\n'

# Retained calls on sex chromosomes must survive preparation before biological QC.
printf '%s\n' \
    'F M 0 0 1 -9 A G A G A G G G' \
    'F W 0 0 2 -9 G G A G A A A G' \
    'F U 0 0 0 -9 A A A A G G 0 0' > chromosomes.ped
printf '%b\n' '1\tauto\t0\t100' '23\tx\t0\t100' '24\ty\t0\t100' '26\tmt\t0\t100' > chromosomes.map
mkdir -p CHROM.imported
plink2 --pedmap chromosomes --make-pgen --threads 1 --memory 1024 --out CHROM.imported/CHROM
plink2 --pfile CHROM.imported/CHROM --export Av --threads 1 --memory 1024 --out before
plink2 --pfile CHROM.imported/CHROM --missing variant-only vcols=nmissdosage,nobs \
    --threads 1 --memory 1024 --out before
Rscript -e '
    library(data.table)
    calls <- fread("before.traw"); missing <- fread("before.vmiss")
    count <- rowSums(!is.na(as.matrix(calls[, 7:ncol(calls), with=FALSE])))
    fast <- missing$OBS_CT - missing$MISSING_DOSAGE_CT
    stopifnot(all(count[calls$CHR != "Y"] == fast[calls$CHR != "Y"]), count[calls$CHR == "Y"] > fast[calls$CHR == "Y"])
'
awk '!/^#/ {print $3}' CHROM.imported/CHROM.pvar > CHROM.retained_markers.txt
awk 'BEGIN {OFS="\t"} !/^#/ {print $3,$3}' CHROM.imported/CHROM.pvar > CHROM.rename_markers.tsv
printf 'final_id\tfinal_chr\tfinal_pos\tfinal_ref\tfinal_alt\tdecision\n' > CHROM.marker_decisions.tsv
awk 'BEGIN {OFS="\t"} !/^#/ {print $3,$1,$2,"A","G","RETAINED_UNIQUE"}' CHROM.imported/CHROM.pvar >> CHROM.marker_decisions.tsv
for chromosome in 1 X Y MT; do
    printf '>%s\n' "$chromosome"
    printf '%01000d\n' 0 | tr '0' 'A'
done > chromosomes.fasta
awk 'BEGIN {OFS="\t"} /^>/ {name=substr($0,2); offset+=length($0)+1; next}
     {print name,length($0),offset,length($0),length($0)+1; offset+=length($0)+1}' chromosomes.fasta > chromosomes.fasta.fai
bash "$repo_root/bin/finalise_target.sh" CHROM CHROM.imported raw top chromosomes.fasta 1
plink2 --pfile CHROM/CHROM --export Av --threads 1 --memory 1024 --out after
Rscript -e '
    library(data.table)
    before <- fread("before.traw"); after <- fread("after.traw")
    after <- after[match(before$SNP, SNP)]
    stopifnot(identical(names(before), names(after)), identical(before$SNP, after$SNP))
    first <- as.matrix(before[, 7:ncol(before), with=FALSE])
    second <- as.matrix(after[, 7:ncol(after), with=FALSE])
    reverse <- which(before$COUNTED != after$COUNTED)
    second[reverse,] <- 2 - second[reverse,,drop=FALSE]
    stopifnot(identical(is.na(first), is.na(second)), isTRUE(all.equal(first, second, check.attributes=FALSE)))
'
cmp CHROM.imported/CHROM.psam CHROM/CHROM.psam
printf 'Autosomal, X, Y and mitochondrial stored calls and sample metadata are preserved.\n'

# An ambiguous A/T SNP can require a strand change despite an unchanged REF/ALT pair.
mkdir ambiguous
cd ambiguous
printf '%s\n' 'F A 0 0 2 -9 A A' 'F B 0 0 1 -9 A T' 'F C 0 0 0 -9 T T' > input.ped
printf '1\tambiguous\t0\t100\n' > input.map
mkdir AMBIG.imported
plink2 --pedmap input --make-pgen --threads 1 --memory 1024 --out AMBIG.imported/AMBIG
printf 'ambiguous\n' > AMBIG.retained_markers.txt
printf 'ambiguous\tambiguous\n' > AMBIG.rename_markers.tsv
printf 'final_id\tfinal_chr\tfinal_pos\tfinal_ref\tfinal_alt\tdecision\nambiguous\t1\t100\tA\tT\tRETAINED_UNIQUE\n' > AMBIG.marker_decisions.tsv
printf '>1\n' > reference.fasta
awk 'BEGIN {for(i=1;i<=1000;i++) printf (i == 99 ? "C" : "A"); printf "\n"}' >> reference.fasta
printf '1\t1000\t3\t1000\t1001\n' > reference.fasta.fai
bash "$repo_root/bin/finalise_target.sh" AMBIG AMBIG.imported raw top reference.fasta 1
plink2 --pfile AMBIG/AMBIG --export Av --threads 1 --memory 1024 --out observed
Rscript -e '
    library(data.table)
    calls <- fread("observed.traw")
    value <- as.numeric(calls[1, 7:9, with=FALSE])
    if (calls$COUNTED == "A") value <- 2 - value
    stopifnot(identical(value, c(2, 1, 0)))
'
printf 'Ambiguous TOP alleles follow the reference-context strand transformation.\n'

# A unique ALT/ALT candidate must not introduce the unobserved genomic REF base.
mkdir ../reference_conflict
cd ../reference_conflict
printf '%s\n' 'F A 0 0 2 -9 C C' 'F B 0 0 1 -9 C G' > conflict.ped
printf '1\tconflict\t0\t100\n' > conflict.map
mkdir panel
cp "$repo_root/tests/data/reference/dbsnp_source/assembly_report.txt" panel/
printf '%s\n' '##fileformat=VCFv4.2' '##contig=<ID=NC_000001.10,length=1000>' > panel/dbsnp.vcf
printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nNC_000001.10\t100\trs100\tA\tC,G\t.\tPASS\t.\n' >> panel/dbsnp.vcf
bash "$repo_root/bin/import_target.sh" CONFLICT ped conflict.ped '' '' DS 1 raw '' '' "$repo_root/bin/target_adapter.pl"
bash "$repo_root/bin/resolve_target_markers.sh" CONFLICT CONFLICT.imported CONFLICT.initial_marker_decisions.tsv raw panel "$repo_root/bin/resolve_target_markers.R" 1
grep -q '^conflict' CONFLICT.retained_markers.txt
if bash "$repo_root/bin/finalise_target.sh" CONFLICT CONFLICT.imported raw top \
    "$repo_root/tests/data/reference/human_g1k_v37.fasta" 1 > conflict.log 2>&1; then
    echo 'An ALT/ALT assay with an incompatible genomic reference was accepted.' >&2
    exit 1
fi
grep -q 'assay/reference conflict' conflict.log
grep -q 'Source marker conflict, rsID rs100' conflict.log
printf 'ALT/ALT matching preserves the fail-fast reference policy.\n'
