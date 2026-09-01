#!/usr/bin/env bash
set -euo pipefail

asset_id="$1"
url="$2"
algorithm="$3"
checksum="$4"
expected_size="$5"
expected_etag="$6"
cached="$7"
embedded="$8"
output="$9"

verify_asset() {
    local value="$1"
    [[ -s "$value" ]] || return 1
    if [[ -n "$expected_size" ]]; then
        [[ "$(stat -c '%s' "$value")" == "$expected_size" ]] || return 1
    fi
    if [[ -n "$checksum" ]]; then
        case "$algorithm" in
            sha256) echo "$checksum  $value" | sha256sum --check --strict --status ;;
            md5) echo "$checksum  $value" | md5sum --check --strict --status ;;
            *) echo "Unsupported checksum algorithm for ${asset_id}: ${algorithm}" >&2; return 1 ;;
        esac
    fi
}

if [[ "$asset_id" == "cache_complete" ]]; then
    printf 'No download required.\n' > "$output"
    exit 0
fi

if [[ "$(basename "$cached")" != "empty_input" ]] && verify_asset "$cached"; then
    ln -s "$(readlink -f "$cached")" "$output"
    exit 0
fi

if [[ "$url" == asset://* ]]; then
    cp "$embedded" "$output"
else
    part="${output}.part"
    etag_file="${output}.etag"
    curl --fail --location --retry 10 --retry-all-errors --connect-timeout 60 \
        --continue-at - --etag-save "$etag_file" --output "$part" "$url"
    if [[ -n "$expected_etag" ]]; then
        observed_etag=$(tr -d '"\r\n' < "$etag_file")
        [[ "$observed_etag" == "$expected_etag" ]] || {
            echo "ETag mismatch for ${asset_id}: expected ${expected_etag}, observed ${observed_etag}" >&2
            exit 11
        }
    fi
    mv "$part" "$output"
fi

if ! verify_asset "$output"; then
    echo "Size or checksum validation failed for ${asset_id}." >&2
    exit 12
fi
