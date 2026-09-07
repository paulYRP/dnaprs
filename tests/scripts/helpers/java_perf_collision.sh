#!/usr/bin/env bash
set -euo pipefail

# Run inside a disposable test container. Lock only this process's new perf file.
case "$DNAPRS_TEST_JAVA_MODE" in
    malformed)
        printf 'not a VCF\n'
        exit 0
        ;;
    plain)
        exec "$DNAPRS_TEST_REAL_JAVA" "$@"
        ;;
    locked|warning)
        perf_dir="/tmp/hsperfdata_$(id -un)"
        mkdir -p "$perf_dir"
        set -o noclobber
        exec 9> "$perf_dir/$BASHPID"
        flock --exclusive --nonblock 9
        ;;
    *)
        echo "Unknown test Java mode: $DNAPRS_TEST_JAVA_MODE" >&2
        exit 2
        ;;
esac

arguments=()
for argument in "$@"; do
    # Force a real runtime warning while retaining the production logging options.
    if [[ "$DNAPRS_TEST_JAVA_MODE" == warning && "$argument" == -jar ]]; then
        arguments+=(-XX:-PerfDisableSharedMem)
    fi
    arguments+=("$argument")
done
exec "$DNAPRS_TEST_REAL_JAVA" "${arguments[@]}"
