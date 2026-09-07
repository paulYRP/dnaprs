#!/usr/bin/env bash
set -Eeuo pipefail

repo_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
prepare_script="${1:-$repo_dir/bin/prepare_plink_reference.sh}"
fixture_dir="$repo_dir/tests/data/reference"
test_dir=$(mktemp -d)
# Synthetic test records must remain readable outside a root-run container.
chmod a+rx "$test_dir"
export UNBREF3_JAR="${UNBREF3_JAR:-/opt/beagle/unbref3.jar}"
export DNAPRS_TEST_REAL_JAVA
DNAPRS_TEST_REAL_JAVA=$(command -v java)
export DNAPRS_TEST_JAVA_HELPER="$repo_dir/tests/scripts/helpers/java_perf_collision.sh"
java() { bash "$DNAPRS_TEST_JAVA_HELPER" "$@"; }
export -f java

sha256sum "$fixture_dir"/bref3_panel/*.bref3 "$fixture_dir/source_panel/chr1.vcf" \
    "$fixture_dir/population.tsv" "$fixture_dir/related.txt" "$UNBREF3_JAR" > "$test_dir/inputs.sha256"
# Exercise both unrelated population subsets without changing the shared fixtures.
awk 'BEGIN { FS=OFS="\t" } $1 == "REF03" { $3="AFR" } { print }' \
    "$fixture_dir/population.tsv" > "$test_dir/population.tsv"
printf 'sample\nREF02\n' > "$test_dir/related.txt"

run_conversion() {
    local label="$1" mode="$2" source_file="$3" chromosome="$4"
    mkdir -p "$test_dir/$label"
    (
        cd "$test_dir/$label"
        printf '%s\n' "$source_file" > sources.txt
        export DNAPRS_TEST_JAVA_MODE="$mode"
        bash "$prepare_script" sources.txt "$test_dir/population.tsv" "$test_dir/related.txt" \
            "$label" 1 GRCh37 "$chromosome" > task.log 2>&1
    )
}

check_conversion() {
    local label="$1" chromosome="$2"
    local task_dir="$test_dir/$label"
    local converted="$task_dir/chr${chromosome}.source.vcf.gz"
    bcftools view --header-only "$converted" > /dev/null
    tabix --list-chroms "$converted" | diff - <(printf '%s\n' "$chromosome")
    bcftools query --list-samples "$converted" | diff - <(bcftools query --list-samples "$fixture_dir/source_panel/chr1.vcf")
    bcftools query --format '%CHROM\t%POS\t%ID\t%REF\t%ALT[\t%GT]\n' "$converted" |
        diff - <(bcftools query --format '%CHROM\t%POS\t%ID\t%REF\t%ALT[\t%GT]\n' \
            "$fixture_dir/source_panel/chr1.vcf" | sed "s/^1\t/$chromosome\t/")
    [[ $(bcftools query --list-samples "$task_dir/chr${chromosome}.eur.vcf.gz" | wc -l) == 22 ]]
    [[ $(bcftools query --list-samples "$task_dir/chr${chromosome}.unrelated.vcf.gz" | wc -l) == 23 ]]
    if bcftools query --list-samples "$task_dir/chr${chromosome}.eur.vcf.gz" | grep -Eq '^REF0[23]$'; then
        echo 'European reference retained an excluded sample.' >&2
        exit 1
    fi
    if bcftools query --list-samples "$task_dir/chr${chromosome}.unrelated.vcf.gz" | grep -qx REF02; then
        echo 'Unrelated reference retained a related sample.' >&2
        exit 1
    fi
    bcftools query --list-samples "$task_dir/chr${chromosome}.unrelated.vcf.gz" | grep -qx REF03
    [[ $(awk 'NR > 1 { print $3, $4, $5 }' "$task_dir/$label.source_qc.tsv") == '22 2 PASS' ]]
    for population in eur all; do
        [[ -s "$task_dir/$label/chromosomes/${population}_chr${chromosome}.pgen" ]]
        awk '!/^#/ { print $2, $3, $4, $5 }' "$task_dir/$label/chromosomes/${population}_chr${chromosome}.pvar" |
            diff - <(printf '100 %s:100:A:G A G\n200 %s:200:C:T C T\n' "$chromosome" "$chromosome")
    done
    if zcat "$converted" | grep -Fq '[warning]'; then
        echo "Java warning entered the converted VCF: $converted" >&2
        exit 1
    fi
}

trap 'echo "Reference conversion test failed; records: $test_dir" >&2' ERR
run_conversion plain plain "$fixture_dir/bref3_panel/chr1.bref3" 1
check_conversion plain 1

# The helper retains an OS file lock across exec, triggering the actual JVM warning.
run_conversion warning warning "$fixture_dir/bref3_panel/chr1.bref3" 1
check_conversion warning 1
grep -F '[warning][perf,memops]' "$test_dir/warning/warning.prepare.log"
grep -Fq 'locked by another process' "$test_dir/warning/warning.prepare.log"

run_conversion concurrent1 locked "$fixture_dir/bref3_panel/chr1.bref3" 1 &
first_pid=$!
run_conversion concurrent2 locked "$fixture_dir/bref3_panel/chr2.bref3" 2 &
second_pid=$!
first_status=0
second_status=0
wait "$first_pid" || first_status=$?
wait "$second_pid" || second_status=$?
[[ "$first_status" == 0 && "$second_status" == 0 ]]
check_conversion concurrent1 1
check_conversion concurrent2 2
if grep -Fq 'locked by another process' "$test_dir"/concurrent*/concurrent*.prepare.log; then
    echo 'Shared Java performance counters were not disabled.' >&2
    exit 1
fi

run_conversion vcf plain "$fixture_dir/source_panel/chr1.vcf" 1
check_conversion vcf 1

mkdir -p "$test_dir/invalid"
printf 'invalid BREF3\n' > "$test_dir/invalid/chr1.bref3"
for failure in invalid malformed; do
    source_file="$fixture_dir/bref3_panel/chr1.bref3"
    mode=malformed
    if [[ "$failure" == invalid ]]; then
        source_file="$test_dir/invalid/chr1.bref3"
        mode=plain
    fi
    if run_conversion "$failure" "$mode" "$source_file" 1; then
        echo "Expected $failure conversion to fail." >&2
        exit 1
    fi
    grep -F "Reference group '$failure', chromosome 1" "$test_dir/$failure/task.log"
    grep -Fq "$source_file" "$test_dir/$failure/task.log"
    grep -Fq "$failure.prepare.log" "$test_dir/$failure/task.log"
    [[ ! -e "$test_dir/$failure/chr1.source.vcf.gz.tbi" ]]
    [[ $(wc -l < "$test_dir/$failure/$failure.source_qc.tsv") == 1 ]]
done
sha256sum --check --strict "$test_dir/inputs.sha256"
printf 'Reference conversion tests passed. Records: %s\n' "$test_dir"
