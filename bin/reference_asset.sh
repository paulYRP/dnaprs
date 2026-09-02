#!/usr/bin/env bash
set -euo pipefail

asset_id="$1"
url="$2"
algorithm="$3"
checksum="$4"
expected_size="$5"
expected_etag="$6"
cached="$7"
cache_target="$8"
embedded="$9"
output="${10}"

verify_asset() {
    local value="$1"
    [[ -s "$value" ]] || return 1
    [[ -n "$checksum" && "$algorithm" == "sha256" ]] || return 1
    if [[ -n "$expected_size" ]]; then
        [[ "$(stat -c '%s' "$value")" == "$expected_size" ]] || return 1
    fi
    echo "$checksum  $value" | sha256sum --check --strict --status
}

link_verified() {
    local value="$1"
    rm -f "$output"
    ln -s "$(readlink -f "$value")" "$output"
}

if [[ "$asset_id" == "cache_complete" ]]; then
    printf 'No download required.\n' > "$output"
    exit 0
fi

if [[ "$(basename "$cached")" != "empty_input" ]] && verify_asset "$cached"; then
    link_verified "$cached"
    exit 0
fi

[[ -n "$cache_target" ]] || { echo "Cache target is empty for ${asset_id}." >&2; exit 10; }
mkdir -p "$(dirname "$cache_target")"

download="${output}.download"
rm -f "$download"
if [[ "$url" == asset://* ]]; then
    cp "$embedded" "$download"
else
    etag_file="${output}.etag"
    curl --fail --location --retry 10 --retry-all-errors --connect-timeout 60 \
        --continue-at - --etag-save "$etag_file" --output "$download" "$url"
    if [[ -n "$expected_etag" ]]; then
        observed_etag=$(tr -d '"\r\n' < "$etag_file")
        [[ "$observed_etag" == "$expected_etag" ]] || {
            echo "ETag mismatch for ${asset_id}: expected ${expected_etag}, observed ${observed_etag}" >&2
            exit 11
        }
    fi
fi

if ! verify_asset "$download"; then
    echo "SHA-256 or size validation failed for ${asset_id}." >&2
    exit 12
fi

# Only a fully verified file can enter the shared cache. A per-asset directory lock
# prevents two concurrent runs from publishing partial or competing copies.
lock_dir="${cache_target}.lock"
lock_acquired=false
for _attempt in $(seq 1 300); do
    if mkdir "$lock_dir" 2>/dev/null; then
        lock_acquired=true
        break
    fi
    if verify_asset "$cache_target"; then
        rm -f "$download"
        link_verified "$cache_target"
        exit 0
    fi
    sleep 1
done
[[ "$lock_acquired" == "true" ]] || { echo "Timed out waiting for reference cache lock: ${cache_target}" >&2; exit 13; }
trap 'rmdir "$lock_dir" 2>/dev/null || true' EXIT

if ! verify_asset "$cache_target"; then
    cache_part="${cache_target}.part.$$"
    rm -f "$cache_part"
    cp "$download" "$cache_part"
    verify_asset "$cache_part" || { echo "Cache staging validation failed for ${asset_id}." >&2; exit 14; }
    mv -f "$cache_part" "$cache_target"
fi
rm -f "$download"
verify_asset "$cache_target" || { echo "Published cache validation failed for ${asset_id}." >&2; exit 15; }
link_verified "$cache_target"
