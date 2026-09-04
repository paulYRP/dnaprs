#!/usr/bin/env bash
set -euo pipefail

repo_root=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
resolver="$repo_root/bin/resolve_inputs.R"
validator="$repo_root/bin/validate_manifests.R"
imputation_script="$repo_root/bin/target_impute_chromosome.sh"
empty="$repo_root/assets/empty_input"
test_root=$(mktemp -d)
resolved_test_root=$(realpath "$test_root")
case "$resolved_test_root" in
    /tmp/tmp.*) ;;
    *) printf 'ERROR: unexpected temporary test path: %s\n' "$resolved_test_root" >&2; exit 1 ;;
esac
trap 'rm -rf -- "$resolved_test_root"' EXIT

reference_root="$test_root/references"
mkdir -p "$reference_root/dbsnp" "$reference_root/map" "$reference_root/panel"

printf 'dbSNP fixture\n' > "$reference_root/dbsnp/GCF_000001405.25.gz"
printf 'assembly fixture\n' > "$reference_root/dbsnp/assembly_report.txt"
printf '>1\nA\n' > "$reference_root/human_g1k_v37.fasta"
printf '1\t2\t3\t2\t3\n' > "$reference_root/human_g1k_v37.fasta.fai"
printf 'sample\tpopulation\nTEST\tEUR\n' > "$reference_root/population.tsv"
printf 'TEST\n' > "$reference_root/related.txt"
printf 'jar\n' > "$reference_root/beagle.test.jar"
printf 'jar\n' > "$reference_root/unbref3.test.jar"
cp "$repo_root/tests/data/reference/sbayesrc_source.zip" "$reference_root/sbayesrc_source.zip"
cp "$repo_root/tests/data/reference/annotation_source.zip" "$reference_root/annotation_source.zip"

for chromosome in $(seq 1 22); do
    printf 'position\trate\n' > "$reference_root/map/plink.chr${chromosome}.GRCh37.map"
    printf 'BREF3 fixture\n' > "$reference_root/panel/chr${chromosome}.synthetic.bref3"
done
printf 'position\trate\n' > "$reference_root/map/chrX.map"
printf 'position\trate\n' > "$reference_root/map/chrY.map"
printf 'position\trate\n' > "$reference_root/map/chrMT.map"

encode_json_string() {
    printf '"%s"' "$1" | base64 | tr -d '\n'
}

run_resolver() {
    local workdir="$1"
    mkdir -p "$workdir"
    (
        cd "$workdir"
        Rscript "$resolver" \
            --input-spec 'W10=' \
            --gwas-spec 'W10=' \
            --reference-spec "$(encode_json_string "$reference_root")" \
            --models-spec 'W10=' \
            --genome-build 'GRCh37' \
            --methods 'plink_ct' \
            --target-imputation 'true' \
            --stop-after 'prs' \
            --reference-only 'true' \
            --reference-mode 'local' \
            --reference-bundle 'test' \
            --phenotype '' \
            --outcome '' \
            --covariates '' \
            --model-type '' \
            --participant-id '' \
            --timepoint-column '' \
            --timepoint-values-spec 'W10=' \
            --group-column '' \
            --control-value '' \
            --case-value '' \
            --beagle-jar '' \
            --unbref3-jar '' \
            --target-staged "$empty" \
            --gwas-staged "$empty" \
            --reference-staged "$reference_root"
    )
}

run_validator() {
    local workdir="$1"
    local resolver_workdir="$2"
    mkdir -p "$workdir"
    (
        cd "$workdir"
        Rscript "$validator" \
            --run-name 'reference-layout' \
            --target-manifest "$resolver_workdir/targets.tsv" \
            --gwas-manifest "$resolver_workdir/gwas.tsv" \
            --reference-manifest "$resolver_workdir/references.tsv" \
            --phenotype-file "$empty" \
            --phenotype-models "$resolver_workdir/models.tsv" \
            --target-assets "$empty" \
            --gwas-assets "$empty" \
            --reference-assets "$reference_root" \
            --methods 'plink_ct' \
            --genome-build 'GRCh37' \
            --seed '1' \
            --report-enabled 'false' \
            --target-imputation 'true' \
            --launch-dir "$repo_root" \
            --reference-base "$reference_root" \
            --run-plan "$resolver_workdir/run_plan.tsv"
    )
}

# A complete autosomal set remains valid when three non-autosomal maps are present.
resolver_valid="$test_root/resolver-valid"
run_resolver "$resolver_valid"
grep -Fq $'GENETIC_MAP_GRCh37\tgenetic_map' "$resolver_valid/references.tsv"
grep -Fq $'IMPUTATION_PANEL_GRCh37\timputation_panel' "$resolver_valid/references.tsv"
run_validator "$test_root/validator-valid" "$resolver_valid"

# Missing autosomal maps fail during discovery and manifest validation.
rm "$reference_root/map/plink.chr22.GRCh37.map"
if run_resolver "$test_root/resolver-missing" 2> "$test_root/resolver-missing.err"; then
    printf 'ERROR: reference discovery accepted a missing chromosome 22 map.\n' >&2
    exit 1
fi
grep -Fq 'missing chromosome(s): 22' "$test_root/resolver-missing.err"
if run_validator "$test_root/validator-missing" "$resolver_valid" 2> "$test_root/validator-missing.err"; then
    printf 'ERROR: manifest validation accepted a missing chromosome 22 map.\n' >&2
    exit 1
fi
grep -Fq "missing genetic-map file(s) for chromosome(s): 22" "$test_root/validator-missing.err"

# Two supported names for one autosome fail instead of selecting one silently.
printf 'position\trate\n' > "$reference_root/map/plink.chr22.GRCh37.map"
printf 'position\trate\n' > "$reference_root/map/chr1.map"
if run_resolver "$test_root/resolver-duplicate" 2> "$test_root/resolver-duplicate.err"; then
    printf 'ERROR: reference discovery accepted duplicate chromosome 1 maps.\n' >&2
    exit 1
fi
grep -Fq 'duplicate chromosome(s): 1' "$test_root/resolver-duplicate.err"
if run_validator "$test_root/validator-duplicate" "$resolver_valid" 2> "$test_root/validator-duplicate.err"; then
    printf 'ERROR: manifest validation accepted duplicate chromosome 1 maps.\n' >&2
    exit 1
fi
grep -Fq "more than one genetic-map file for chromosome(s): 1" "$test_root/validator-duplicate.err"

# Chromosome imputation reports the same missing and ambiguous map conditions.
runtime_root="$test_root/runtime"
mkdir -p "$runtime_root/panel" "$runtime_root/map-missing" "$runtime_root/map-duplicate"
printf 'BREF3 fixture\n' > "$runtime_root/panel/chr1.bref3"
printf 'target\n' > "$runtime_root/target.pgen"
printf 'target\n' > "$runtime_root/target.pvar"
printf 'target\n' > "$runtime_root/target.psam"
printf 'jar\n' > "$runtime_root/beagle.jar"
printf '>1\nA\n' > "$runtime_root/reference.fasta"

if (
    cd "$runtime_root"
    BEAGLE_JAR="$runtime_root/beagle.jar" bash "$imputation_script" \
        TEST 1 "$runtime_root/target.pgen" "$runtime_root/target.pvar" "$runtime_root/target.psam" \
        "$runtime_root/panel" "$runtime_root/map-missing" 0.8 1 1024 GRCh37 "$runtime_root/reference.fasta"
) 2> "$test_root/runtime-missing.err"; then
    printf 'ERROR: chromosome imputation accepted a missing chromosome 1 map.\n' >&2
    exit 1
fi
grep -Fq 'has no supported map for chromosome 1' "$test_root/runtime-missing.err"

printf 'position\trate\n' > "$runtime_root/map-duplicate/plink.chr1.GRCh37.map"
printf 'position\trate\n' > "$runtime_root/map-duplicate/chr1.map"
if (
    cd "$runtime_root"
    BEAGLE_JAR="$runtime_root/beagle.jar" bash "$imputation_script" \
        TEST 1 "$runtime_root/target.pgen" "$runtime_root/target.pvar" "$runtime_root/target.psam" \
        "$runtime_root/panel" "$runtime_root/map-duplicate" 0.8 1 1024 GRCh37 "$runtime_root/reference.fasta"
) 2> "$test_root/runtime-duplicate.err"; then
    printf 'ERROR: chromosome imputation accepted duplicate chromosome 1 maps.\n' >&2
    exit 1
fi
grep -Fq 'has 2 supported maps for chromosome 1' "$test_root/runtime-duplicate.err"
