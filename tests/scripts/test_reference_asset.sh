#!/usr/bin/env bash
set -euo pipefail

repo_root=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
asset_script="$repo_root/bin/reference_asset.sh"
fixture="$repo_root/assets/reference/human_g1k_v37.fasta.fai"
empty="$repo_root/assets/empty_input"
test_root=$(mktemp -d)
trap 'rm -rf "$test_root"' EXIT

checksum=$(sha256sum "$fixture" | awk '{print $1}')
size=$(stat -c '%s' "$fixture")
cache="$test_root/cache/reference.fai"

run_asset() {
    local workdir="$1"
    local asset_id="$2"
    local expected_checksum="$3"
    local cached="$4"
    mkdir -p "$workdir"
    (
        cd "$workdir"
        bash "$asset_script" "$asset_id" 'asset://fixture' sha256 "$expected_checksum" \
            "$size" '' "$cached" "$cache" "$fixture" "$asset_id"
    )
}

# A missing cache is populated only after checksum and size validation.
run_asset "$test_root/first" reference_fai "$checksum" "$empty"
test -s "$cache"
test "$(sha256sum "$cache" | awk '{print $1}')" = "$checksum"

# A valid cache is reused without replacing its content.
before=$(stat -c '%Y:%s' "$cache")
run_asset "$test_root/reuse" reference_fai "$checksum" "$cache"
after=$(stat -c '%Y:%s' "$cache")
test "$before" = "$after"

# An invalid cached target is replaced atomically by the verified fixture.
printf 'corrupt\n' > "$cache"
run_asset "$test_root/replace" reference_fai "$checksum" "$empty"
test "$(sha256sum "$cache" | awk '{print $1}')" = "$checksum"

# A wrong checksum is rejected and cannot create a shared cache entry.
rm -f "$cache"
if run_asset "$test_root/reject" reference_fai "$(printf '0%.0s' {1..64})" "$empty"; then
    echo 'ERROR: corrupt reference fixture was accepted.' >&2
    exit 1
fi
test ! -e "$cache"

# Concurrent writers converge on the same verified cache object.
run_asset "$test_root/concurrent-a" reference_a "$checksum" "$empty" &
pid_a=$!
run_asset "$test_root/concurrent-b" reference_b "$checksum" "$empty" &
pid_b=$!
wait "$pid_a"
wait "$pid_b"
test "$(sha256sum "$cache" | awk '{print $1}')" = "$checksum"
test -L "$test_root/concurrent-a/reference_a"
test -L "$test_root/concurrent-b/reference_b"
