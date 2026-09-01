process PLINK_REFERENCE_FREQ {
    tag "${reference.reference_id}"
    label 'process_high'
    label 'process_plink'

    container 'ghcr.io/paulyrp/dnaprs-plink2:2.0.0-a.6.12'

    input:
    tuple val(reference), path(reference_files)

    output:
    tuple val(reference), path('reference.afreq'), path('reference.log'), emit: frequency
    path 'versions.yml', emit: versions

    script:
    """
    plink2 \
        --pfile '${reference.path}' \
        --freq \
        --threads '${task.cpus}' \
        --out reference

    cat > versions.yml <<-VERSIONS
    "${task.process}:${reference.reference_id}":
        plink2: \$(plink2 --version 2>&1 | head -n 1 | cut -d ' ' -f 2 | sed 's/^v//')
    VERSIONS
    """

    stub:
    """
    printf '#CHROM\tID\tREF\tALT\tALT_FREQS\tOBS_CT\n1\t1:100:A:G\tA\tG\t0.25\t4\n' > reference.afreq
    printf 'PLINK reference frequency stub\n' > reference.log
    printf '"${task.process}:${reference.reference_id}":\n  plink2: stub\n' > versions.yml
    """
}
