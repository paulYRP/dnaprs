# Raw-format test fixtures

These files are deterministic PLINK 2 exports of `../test_target.vcf`. They are
small test inputs, not biological reference data.

The six-participant fixture includes four additional markers for valid LD-pruned
participant QC. The two GWAS scoring markers and their genotypes are unchanged.
Regenerate the PLINK and BGEN files with `tests/scripts/build_target_fixtures.sh`
inside the existing imputation container.

- `target_bed.{bed,bim,fam}`: PLINK 1 binary trio.
- `target_ped.{ped,map}`: PED/MAP text pair.
- `target_bgen.{bgen,sample}`: BGEN 1.2 with 16-bit probabilities.

The BGEN fixture uses 16-bit probabilities because the tested PLINK alpha 6.12
build cannot read its own 8-bit export for this very small dataset.
