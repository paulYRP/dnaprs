process BUILD_CT_WEIGHTS {
    tag "${meta.trait_id}"
    label 'process_single'

    container 'ghcr.io/paulyrp/dnaprs-analysis:1.0.0'

    input:
    tuple val(meta), path(cojo), path(clumps), path(harmonisation_qc), path(clump_log)
    path weight_script

    output:
    tuple val(meta), path("${meta.trait_id}.plink_ct.weights.tsv"), path("${meta.trait_id}.plink_ct.weight_qc.tsv"), path(harmonisation_qc), path(clump_log), emit: weights
    tuple val("${task.process}"), val('R'), eval("Rscript -e 'cat(as.character(getRversion()))' 2>/dev/null || printf stub"), emit: versions_r, topic: versions
    tuple val("${task.process}"), val('data.table'), eval("Rscript -e 'cat(as.character(packageVersion(\"data.table\")))' 2>/dev/null || printf stub"), emit: versions_data_table, topic: versions
    script:
    """
    Rscript ${weight_script} \
        --cojo '${cojo}' \
        --clumps '${clumps}' \
        --trait-id '${meta.trait_id}' \
        --prs-name '${meta.prs_name}'
    """

    stub:
    """
    printf 'SNP\tA1\tBETA\n1:100:A:G\tG\t0.10\n' > ${meta.trait_id}.plink_ct.weights.tsv
    printf 'trait_id\tprs_name\tharmonised_variants\tclumped_variants\tretained_percent\tstatus\n${meta.trait_id}\t${meta.prs_name}\t1\t1\t100\tPASS\n' > ${meta.trait_id}.plink_ct.weight_qc.tsv
    """
}
