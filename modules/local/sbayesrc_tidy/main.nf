process SBAYESRC_TIDY {
    tag "${meta.trait_id}"
    label 'process_high'
    label 'process_r'

    input:
    tuple val(meta), path(cojo), path(clump_input), path(harmonisation_qc), val(ld_reference)
    path sbayesrc_script

    output:
    tuple val(meta), path("${meta.trait_id}.tidy.ma"), path("${meta.trait_id}.tidy_qc.tsv"), path(harmonisation_qc), emit: tidy
    tuple val(meta), path("${meta.trait_id}.sbayesrc.tidy.log"), emit: logs
    path 'versions.yml', emit: versions

    script:
    """
    export OMP_NUM_THREADS='${task.cpus}'
    Rscript ${sbayesrc_script} \
        --action tidy \
        --input '${cojo}' \
        --ld-dir '${ld_reference.path}' \
        --trait-id '${meta.trait_id}'

    cat > versions.yml <<-VERSIONS
    "${task.process}:${meta.trait_id}":
        SBayesRC: \$(Rscript -e "cat(as.character(packageVersion('SBayesRC')))" 2>/dev/null)
        R: \$(Rscript -e 'cat(as.character(getRversion()))')
    VERSIONS
    """

    stub:
    """
    cp ${cojo} ${meta.trait_id}.tidy.ma
    printf 'trait_id\tinput_variants\tld_aligned_variants\tretained_percent\treview_below_70_percent\n${meta.trait_id}\t1\t1\t100\tFALSE\n' > ${meta.trait_id}.tidy_qc.tsv
    printf 'SBayesRC tidy stub\n' > ${meta.trait_id}.sbayesrc.tidy.log
    printf '"${task.process}:${meta.trait_id}":\n  SBayesRC: stub\n  R: stub\n' > versions.yml
    """
}
