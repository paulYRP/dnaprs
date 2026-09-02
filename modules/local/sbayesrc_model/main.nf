process SBAYESRC_MODEL {
    tag "${meta.trait_id}"
    label 'process_high'

    container 'docker.io/zhiliz/sbayesrc:0.2.6@sha256:5a6139e6c3ab471799c059fe7032870f22bacc07dfadf14fc0517314db5ea4bb'

    input:
    tuple val(meta), path(imputed), path(impute_qc), path(tidy_qc), path(harmonisation_qc), val(ld_reference), path(ld_reference_files), val(annotation_reference), path(annotation_reference_files)
    path sbayesrc_script
    val seed

    output:
    tuple val(meta), path("${meta.trait_id}.sbayesrc.txt"), path("${meta.trait_id}.sbayesrc.par"), path("${meta.trait_id}.sbayesrc.model_qc.tsv"), path(impute_qc), path(tidy_qc), path(harmonisation_qc), emit: model
    tuple val(meta), path("${meta.trait_id}.sbayesrc.model.log"), emit: logs
    tuple val("${task.process}"), val('SBayesRC'), eval("Rscript -e 'cat(as.character(packageVersion(\"SBayesRC\")))' 2>/dev/null || printf stub"), emit: versions_sbayesrc, topic: versions
    tuple val("${task.process}"), val('R'), eval("Rscript -e 'cat(as.character(getRversion()))' 2>/dev/null || printf stub"), emit: versions_r, topic: versions
    script:
    """
    export OMP_NUM_THREADS='${task.cpus}'
    Rscript ${sbayesrc_script} \
        --action model \
        --input '${imputed}' \
        --ld-dir '${ld_reference.path}' \
        --annotation '${annotation_reference.path}' \
        --trait-id '${meta.trait_id}' \
        --seed '${seed}'
    """

    stub:
    """
    printf 'SNP\tA1\tBETA\tPIP\n1:100:A:G\tG\t0.10\t0.20\n' > ${meta.trait_id}.sbayesrc.txt
    printf 'parameter\tvalue\nheritability\t0.10\n' > ${meta.trait_id}.sbayesrc.par
    printf 'trait_id\tweights\tnonzero_weights\tmaximum_pip\tstatus\n${meta.trait_id}\t1\t1\t0.20\tPASS\n' > ${meta.trait_id}.sbayesrc.model_qc.tsv
    printf 'SBayesRC model stub\n' > ${meta.trait_id}.sbayesrc.model.log
    """
}
