#!/usr/bin/env bash
set -euo pipefail

repo_root=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
asset_script="$repo_root/bin/reference_asset.sh"
fixture="$repo_root/assets/reference/human_g1k_v37.fasta.fai"
empty="$repo_root/assets/empty_input"
test_root=$(mktemp -d)
server_pid=''
cleanup() {
    if [[ -n "$server_pid" ]]; then
        kill "$server_pid" 2>/dev/null || true
        wait "$server_pid" 2>/dev/null || true
    fi
    rm -rf "$test_root"
}
trap cleanup EXIT

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

# Use real curl retries against a small local HTTP server, without remote downloads.
python3 "$repo_root/tests/scripts/helpers/reference_http_server.py" "$fixture" \
    > "$test_root/server.url" 2> "$test_root/server.log" &
server_pid=$!
for _attempt in $(seq 1 50); do
    [[ -s "$test_root/server.url" ]] && break
    sleep 0.1
done
if [[ ! -s "$test_root/server.url" ]]; then
    cat "$test_root/server.log" >&2
    echo 'HTTP reference fixture did not start.' >&2
    exit 1
fi
read -r base_url < "$test_root/server.url"
export NO_PROXY="${NO_PROXY:+$NO_PROXY,}127.0.0.1"
export no_proxy="$NO_PROXY"

run_http_asset() {
    local mode="$1"
    local expected_exit="$2"
    local expected_message="${3:-}"
    local workdir="$test_root/http-$mode"
    local http_cache="$test_root/cache/$mode.fai"
    local exit_code=0
    mkdir -p "$workdir"
    (
        cd "$workdir"
        bash "$asset_script" reference_fai "$base_url/$mode" sha256 "$checksum" \
            "$size" 'reference-fixture-v1' "$empty" "$http_cache" "$fixture" reference_fai
    ) > "$workdir/task.log" 2>&1 || exit_code=$?

    if [[ "$exit_code" != "$expected_exit" ]]; then
        cat "$workdir/task.log" >&2
        printf 'Unexpected exit status for %s: %s (expected %s).\n' "$mode" "$exit_code" "$expected_exit" >&2
        exit 1
    fi
    if [[ "$expected_exit" == 0 ]]; then
        cmp "$fixture" "$http_cache"
        cmp "$fixture" "$workdir/reference_fai"
        test -L "$workdir/reference_fai"
    else
        grep -Fq "$expected_message" "$workdir/task.log"
        test ! -e "$http_cache"
        test ! -e "$workdir/reference_fai"
    fi
    printf 'PASS: HTTP reference %s\n' "$mode"
}

run_http_asset valid 0
run_http_asset retry_once 0
run_http_asset retry_four 0
grep -Fq 'retrying all errors' "$test_root/http-retry_once/task.log"
grep -Fq 'retrying all errors' "$test_root/http-retry_four/task.log"
run_http_asset conflicting 11 'ETag mismatch'
run_http_asset incorrect 11 'ETag mismatch'
run_http_asset missing 11 'ETag mismatch'
run_http_asset bad_checksum 12 'SHA-256 or size validation failed'
run_http_asset bad_size 12 'SHA-256 or size validation failed'

# A verified HTTP download can be reused with an unavailable source URL.
mkdir -p "$test_root/http-reuse"
(
    cd "$test_root/http-reuse"
    bash "$asset_script" reference_fai 'http://127.0.0.1:1/unavailable' sha256 "$checksum" \
        "$size" 'reference-fixture-v1' "$test_root/cache/valid.fai" \
        "$test_root/cache/valid.fai" "$fixture" reference_fai
)
cmp "$fixture" "$test_root/http-reuse/reference_fai"
printf 'Verified reference caching and HTTP download tests passed.\n'
