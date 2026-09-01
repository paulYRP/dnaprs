process SUMMARISE_GENERATION_QC {
    tag 'variant flow'
    label 'process_low'
    label 'process_r'

    container 'ghcr.io/paulyrp/dnaprs-analysis:1.0.0'

    input:
    path generation_qc_files
    path score_qc
    path generation_qc_script

    output:
    path 'variant_flow.tsv', emit: variant_flow
    path 'versions.yml', emit: versions

    script:
    qcARG = generation_qc_files.collect { qc_file -> qc_file.toString() }.join(',')
    """
    Rscript ${generation_qc_script} \
        --qc-files '${qcARG}' \
        --score-qc '${score_qc}'

    cat > versions.yml <<-VERSIONS
    "${task.process}":
        R: \$(Rscript -e 'cat(as.character(getRversion()))')
        data.table: \$(Rscript -e "cat(as.character(packageVersion('data.table')))" 2>/dev/null)
    VERSIONS
    """

    stub:
    """
    printf 'cohort\trole\ttrait_id\tprs_name\tmethod\tstage\tstage_order\tvariant_count\tpercent_of_source\tpercent_of_previous\nTEST\ttarget\tMDD\tMDD_PRS\tplink_ct\tSource GWAS\t1\t1\t100\t100\nTEST\ttarget\tMDD\tMDD_PRS\tplink_ct\tScored\t4\t1\t100\t100\n' > variant_flow.tsv
    printf '"${task.process}":\n  R: stub\n  data.table: stub\n' > versions.yml
    """
}
