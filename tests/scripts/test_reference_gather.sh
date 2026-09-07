#!/usr/bin/env bash
set -euo pipefail

repo_root=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
test_root=$(mktemp -d)
cd "$test_root"
mkdir -p group/chromosomes
printf 'chromosome\tsource\tsamples\tvariants\tstatus\tsource_sha256\n1\tchr1.vcf\t2\t1\tPASS\tfixture\n' > group.source_qc.tsv
printf 'fixture\n' > group.prepare.log

if bash "$repo_root/bin/assemble_plink_reference.sh" reference 1 GRCh37 1,2 population.tsv related.txt group 2> missing.err; then
    printf 'ERROR: reference merge accepted a missing chromosome.\n' >&2
    exit 1
fi
grep -Fq 'Reference gather expected chromosomes 1,2; received 1.' missing.err

if bash "$repo_root/bin/assemble_plink_reference.sh" reference_duplicate 1 GRCh37 1 population.tsv related.txt group group 2> duplicate.err; then
    printf 'ERROR: reference merge accepted a repeated chromosome.\n' >&2
    exit 1
fi
grep -Fq 'Reference gather expected chromosomes 1; received 1,1.' duplicate.err
printf 'Reference gather rejects missing and repeated chromosomes before PLINK.\n'
