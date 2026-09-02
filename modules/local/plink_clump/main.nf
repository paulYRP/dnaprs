process PLINK_CLUMP {
    tag "${meta.trait_id}"
    label 'process_high'

    container 'ghcr.io/paulyrp/dnaprs-plink2:2.0.0-a.6.12'

    input:
    tuple val(meta), path(cojo), path(clump_input), path(harmonisation_qc), val(reference), path(reference_files)

    output:
    tuple val(meta), path(cojo), path("${meta.trait_id}.clumps"), path(harmonisation_qc), path("${meta.trait_id}.clump.log"), emit: clumped
    path 'versions.yml', emit: versions

    script:
    """
    plink2 \
        --pfile '${reference.path}' \
        --clump '${clump_input}' \
        --clump-id-field ID \
        --clump-p-field P \
        --clump-a1-field A1 \
        --clump-force-a1 \
        --clump-p1 1 \
        --clump-p2 1 \
        --clump-r2 0.10 \
        --clump-kb 250 \
        --threads '${task.cpus}' \
        --out '${meta.trait_id}.clump'

    mv '${meta.trait_id}.clump.clumps' '${meta.trait_id}.clumps'

    cat > versions.yml <<-VERSIONS
    "${task.process}:${meta.trait_id}":
        plink2: \$(plink2 --version 2>&1 | head -n 1 | cut -d ' ' -f 2 | sed 's/^v//')
    VERSIONS
    """

    stub:
    """
    printf '#CHROM\tPOS\tID\tREF\tALT\tPROVISIONAL_REF?\tA1\tF\tP\tTOTAL\tNONMAJOR\n1\t100\t1:100:A:G\tA\tG\tN\tG\t1\t0.001\t1\t1\n' > ${meta.trait_id}.clumps
    printf 'PLINK clumping stub\n' > ${meta.trait_id}.clump.log
    printf '"${task.process}:${meta.trait_id}":\n  plink2: stub\n' > versions.yml
    """
}
