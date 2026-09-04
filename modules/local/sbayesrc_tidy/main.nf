process SBAYESRC_TIDY {
    tag "${meta.trait_id}"
    label 'process_high'

    container 'docker.io/zhiliz/sbayesrc:0.2.6'

    input:
    tuple val(meta), path(cojo), path(clump_input), path(harmonisation_qc), val(ld_reference), path(ld_reference_files)
    path sbayesrc_script

    output:
    tuple val(meta), path("${meta.trait_id}.tidy.ma"), path("${meta.trait_id}.tidy_qc.tsv"), path(harmonisation_qc), emit: tidy
    tuple val(meta), path("${meta.trait_id}.sbayesrc.tidy.log"), emit: logs
    tuple val("${task.process}"), val('SBayesRC'), eval("Rscript -e 'cat(as.character(packageVersion(\"SBayesRC\")))' 2>/dev/null || printf stub"), emit: versions_sbayesrc, topic: versions
    tuple val("${task.process}"), val('R'), eval("Rscript -e 'cat(as.character(getRversion()))' 2>/dev/null || printf stub"), emit: versions_r, topic: versions
    script:
    """
    export OMP_NUM_THREADS='${task.cpus}'
    Rscript ${sbayesrc_script} \
        --action tidy \
        --input '${cojo}' \
        --ld-dir '${ld_reference.path}' \
        --trait-id '${meta.trait_id}'
    """

    stub:
    """
    cp ${cojo} ${meta.trait_id}.tidy.ma
    printf 'trait_id\tinput_variants\tld_aligned_variants\tretained_percent\treview_below_70_percent\n${meta.trait_id}\t1\t1\t100\tFALSE\n' > ${meta.trait_id}.tidy_qc.tsv
    printf 'SBayesRC tidy stub\n' > ${meta.trait_id}.sbayesrc.tidy.log
    """
}
