process PLINK_REFERENCE_FREQ {
    tag "${reference.reference_id}"
    label 'process_high'

    container 'ghcr.io/paulyrp/dnaprs-plink2:2.0.0-a.6.12'

    input:
    tuple val(reference), path(reference_files)

    output:
    tuple val(reference), path('reference.afreq'), path('reference.log'), emit: frequency
    tuple val("${task.process}"), val('plink2'), eval("command -v plink2 >/dev/null && plink2 --version 2>&1 | head -n 1 | cut -d ' ' -f 2 | sed 's/^v//' || printf stub"), emit: versions_plink2, topic: versions
    script:
    """
    plink2 \
        --pfile '${reference.path}' \
        --freq \
        --threads '${task.cpus}' \
        --out reference
    """

    stub:
    """
    printf '#CHROM\tID\tREF\tALT\tALT_FREQS\tOBS_CT\n1\t1:100:A:G\tA\tG\t0.25\t4\n' > reference.afreq
    printf 'PLINK reference frequency stub\n' > reference.log
    """
}
