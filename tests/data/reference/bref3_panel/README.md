# BREF3 conversion fixtures

`chr1.bref3` and `chr2.bref3` contain the 24 synthetic samples and two phased SNPs
from [the VCF fixture](../source_panel/chr1.vcf). Chromosome 2 repeats these records
with chromosome labels changed to 2. These files contain no study data.

To regenerate them, run from the repository root in the imputation container:

```bash
docker run --rm --volume "$PWD:/work/dnaprs" --workdir /work/dnaprs \
    ghcr.io/paulyrp/dnaprs-imputation:1.1.0 \
    bash tests/scripts/generate_bref3_fixtures.sh
```

The generator downloads the official Beagle `bref3.27Feb25.75f.jar`, verifies its
SHA-256 checksum and converts the VCF records. The converter is a fixture-generation
tool only. Tests use the existing container's unbref3 JAR and do not download it.

The conversion tests check sample order, variant positions, identifiers, alleles
and phased genotypes against the synthetic VCF. They also check sample filtering,
concurrent conversions, Java diagnostic separation and invalid-output rejection.
