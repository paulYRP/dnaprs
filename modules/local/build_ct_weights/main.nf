process BUILD_CT_WEIGHTS {
    tag "${meta.trait_id}"
    label 'process_low'
    label 'process_r'

    container 'ghcr.io/paulyrp/dnaprs-analysis:1.0.0'

    input:
    tuple val(meta), path(cojo), path(clumps), path(harmonisation_qc), path(clump_log)
    path weight_script

    output:
    tuple val(meta), path("${meta.trait_id}.plink_ct.weights.tsv"), path("${meta.trait_id}.plink_ct.weight_qc.tsv"), path(harmonisation_qc), path(clump_log), emit: weights
    path 'versions.yml', emit: versions

    script:
    """
    Rscript ${weight_script} \
        --cojo '${cojo}' \
        --clumps '${clumps}' \
        --trait-id '${meta.trait_id}' \
        --prs-name '${meta.prs_name}'

    cat > versions.yml <<-VERSIONS
    "${task.process}:${meta.trait_id}":
        R: \$(Rscript -e 'cat(as.character(getRversion()))')
        data.table: \$(Rscript -e "cat(as.character(packageVersion('data.table')))" 2>/dev/null)
    VERSIONS
    """

    stub:
    """
    printf 'SNP\tA1\tBETA\n1:100:A:G\tG\t0.10\n' > ${meta.trait_id}.plink_ct.weights.tsv
    printf 'trait_id\tprs_name\tharmonised_variants\tclumped_variants\tretained_percent\tstatus\n${meta.trait_id}\t${meta.prs_name}\t1\t1\t100\tPASS\n' > ${meta.trait_id}.plink_ct.weight_qc.tsv
    printf '"${task.process}:${meta.trait_id}":\n  R: stub\n  data.table: stub\n' > versions.yml
    """
}
