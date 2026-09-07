#!/usr/bin/env bash
set -euo pipefail

repo_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
fixture_dir="$repo_dir/tests/data/reference/bref3_panel"
test_dir=$(mktemp -d)
converter="$test_dir/bref3.27Feb25.75f.jar"

curl --fail --location --retry 3 \
    https://faculty.washington.edu/browning/beagle/bref3.27Feb25.75f.jar \
    --output "$converter"
printf '%s  %s\n' 6166426f63b2c1cfed9e2cda2f9ba59ef3b3c19fdfc988a9ec6036358457e3f2 "$converter" | sha256sum --check --strict

mkdir -p "$fixture_dir"
for chromosome in 1 2; do
    sed -E "s/ID=1,/ID=$chromosome,/;s/^1[[:space:]]/$chromosome\t/" \
        "$repo_dir/tests/data/reference/source_panel/chr1.vcf" > "$test_dir/chr${chromosome}.vcf"
    java -XX:+PerfDisableSharedMem -Xlog:disable -Xlog:all=warning:stderr \
        -jar "$converter" "$test_dir/chr${chromosome}.vcf" > "$fixture_dir/chr${chromosome}.bref3"
done
sha256sum "$fixture_dir"/*.bref3
