#!/usr/bin/env bash
set -euo pipefail

repo_root=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
failed=0

while IFS= read -r module_dir; do
    module=$(basename "$module_dir")
    main="$module_dir/main.nf"
    meta="$module_dir/meta.yml"

    for required in "$main" "$meta"; do
        if [[ ! -s "$required" ]]; then
            printf 'ERROR: %s is missing or empty.\n' "${required#"$repo_root/"}" >&2
            failed=1
        fi
    done
    [[ -s "$main" && -s "$meta" ]] || continue

    label_count=$(grep -Ec "^[[:space:]]+label 'process_(single|low|medium|high)'$" "$main" || true)
    if [[ "$label_count" -ne 1 ]]; then
        printf 'ERROR: %s must use exactly one bundled process label.\n' "$module" >&2
        failed=1
    fi
    if ! grep -Eq '^[[:space:]]+stub:' "$main"; then
        printf 'ERROR: %s has no stub block.\n' "$module" >&2
        failed=1
    fi
    if grep -Eq 'versions[.]yml.*topic:[[:space:]]*versions' "$main"; then
        printf 'ERROR: %s still publishes a legacy versions.yml topic record.\n' "$module" >&2
        failed=1
    fi
    case "$module" in
        collect_versions|public_figures|render_report) ;;
        *)
            if ! grep -Eq 'tuple val[(]"[$][{]task[.]process[}]"[)], val[(].*[)], (eval|val)[(].*[)], emit:[[:space:]]*versions_[a-z0-9_]+,[[:space:]]*topic:[[:space:]]*versions' "$main"; then
                printf 'ERROR: %s does not publish nf-core version tuples to the versions topic.\n' "$module" >&2
                failed=1
            fi
            if ! grep -Eq '^topics:' "$meta" || ! grep -Eq '^[[:space:]]+versions:' "$meta"; then
                printf 'ERROR: %s metadata does not document the versions topic.\n' "$module" >&2
                failed=1
            fi
            ;;
    esac
    for key in authors maintainers; do
        if ! grep -Eq "^${key}:" "$meta"; then
            printf 'ERROR: %s metadata has no %s entry.\n' "$module" "$key" >&2
            failed=1
        fi
    done
done < <(find "$repo_root/modules/local" -mindepth 1 -maxdepth 1 -type d | sort)

if find "$repo_root/modules/local" -mindepth 1 -maxdepth 1 -type d -empty | grep -q .; then
    printf 'ERROR: empty local module directories are not allowed.\n' >&2
    failed=1
fi

if grep -RHE '@sha256:' --include='main.nf' "$repo_root/modules/local"; then
    printf 'ERROR: process container references must not use tag@digest syntax.\n' >&2
    failed=1
fi

exit "$failed"
