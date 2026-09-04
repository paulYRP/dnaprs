process SBAYESRC_IMPUTE {
    tag "${meta.trait_id}"
    label 'process_high'

    container 'docker.io/zhiliz/sbayesrc:0.2.6'

    input:
    tuple val(meta), path(tidy), path(tidy_qc), path(harmonisation_qc), val(ld_reference), path(ld_reference_files)
    path sbayesrc_script

    output:
    tuple val(meta), path("${meta.trait_id}.imputed.ma"), path("${meta.trait_id}.impute_qc.tsv"), path(tidy_qc), path(harmonisation_qc), emit: imputed
    tuple val(meta), path("${meta.trait_id}.sbayesrc.impute.log"), emit: logs
    tuple val("${task.process}"), val('SBayesRC'), eval("Rscript -e 'cat(as.character(packageVersion(\"SBayesRC\")))' 2>/dev/null || printf stub"), emit: versions_sbayesrc, topic: versions
    tuple val("${task.process}"), val('R'), eval("Rscript -e 'cat(as.character(getRversion()))' 2>/dev/null || printf stub"), emit: versions_r, topic: versions
    script:
    """
    export OMP_NUM_THREADS='${task.cpus}'
    Rscript ${sbayesrc_script} \
        --action impute \
        --input '${tidy}' \
        --ld-dir '${ld_reference.path}' \
        --trait-id '${meta.trait_id}'
    """

    stub:
    """
    cp ${tidy} ${meta.trait_id}.imputed.ma
    printf 'trait_id\tld_aligned_variants\tsummary_imputed_variants\tstatus\n${meta.trait_id}\t1\t1\tPASS\n' > ${meta.trait_id}.impute_qc.tsv
    printf 'SBayesRC impute stub\n' > ${meta.trait_id}.sbayesrc.impute.log
    """
}
