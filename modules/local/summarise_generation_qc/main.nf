process SUMMARISE_GENERATION_QC {
    tag 'variant flow'
    label 'process_single'

    container 'ghcr.io/paulyrp/dnaprs-analysis:1.0.0'

    input:
    path generation_qc_files
    path score_qc
    path generation_qc_script

    output:
    path 'variant_flow.tsv', emit: variant_flow
    tuple val("${task.process}"), val('R'), eval("Rscript -e 'cat(as.character(getRversion()))' 2>/dev/null || printf stub"), emit: versions_r, topic: versions
    tuple val("${task.process}"), val('data.table'), eval("Rscript -e 'cat(as.character(packageVersion(\"data.table\")))' 2>/dev/null || printf stub"), emit: versions_data_table, topic: versions
    script:
    qcARG = generation_qc_files.collect { qc_file -> qc_file.toString() }.join(',')
    """
    Rscript ${generation_qc_script} \
        --qc-files '${qcARG}' \
        --score-qc '${score_qc}'
    """

    stub:
    """
    printf 'cohort\trole\ttrait_id\tprs_name\tmethod\tstage\tstage_order\tvariant_count\tpercent_of_source\tpercent_of_previous\nTEST\ttarget\tMDD\tMDD_PRS\tplink_ct\tSource GWAS\t1\t1\t100\t100\nTEST\ttarget\tMDD\tMDD_PRS\tplink_ct\tScored\t4\t1\t100\t100\n' > variant_flow.tsv
    """
}
