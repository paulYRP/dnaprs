process SBAYESRC_IMPUTE {
    tag "${meta.trait_id}"
    label 'process_high'
    label 'process_r'

    container 'docker.io/zhiliz/sbayesrc:0.2.6@sha256:5a6139e6c3ab471799c059fe7032870f22bacc07dfadf14fc0517314db5ea4bb'

    input:
    tuple val(meta), path(tidy), path(tidy_qc), path(harmonisation_qc), val(ld_reference), path(ld_reference_files)
    path sbayesrc_script

    output:
    tuple val(meta), path("${meta.trait_id}.imputed.ma"), path("${meta.trait_id}.impute_qc.tsv"), path(tidy_qc), path(harmonisation_qc), emit: imputed
    tuple val(meta), path("${meta.trait_id}.sbayesrc.impute.log"), emit: logs
    path 'versions.yml', emit: versions

    script:
    """
    export OMP_NUM_THREADS='${task.cpus}'
    Rscript ${sbayesrc_script} \
        --action impute \
        --input '${tidy}' \
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
    cp ${tidy} ${meta.trait_id}.imputed.ma
    printf 'trait_id\tld_aligned_variants\tsummary_imputed_variants\tstatus\n${meta.trait_id}\t1\t1\tPASS\n' > ${meta.trait_id}.impute_qc.tsv
    printf 'SBayesRC impute stub\n' > ${meta.trait_id}.sbayesrc.impute.log
    printf '"${task.process}:${meta.trait_id}":\n  SBayesRC: stub\n  R: stub\n' > versions.yml
    """
}
