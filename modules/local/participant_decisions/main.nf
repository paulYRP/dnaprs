process PARTICIPANT_DECISIONS {
    tag "${meta.cohort}"
    label 'process_single'

    container 'ghcr.io/paulyrp/dnaprs-analysis:1.0.1'

    input:
    tuple val(meta), path(sample_decisions), path(eda_tables), path(target_ancestry)
    path decision_script
    val sample_missingness
    val heterozygosity_z_threshold

    output:
    tuple val(meta), path("${meta.cohort}.participant_decisions.tsv"), path("${meta.cohort}.score_eligible.keep"), emit: decisions
    tuple val("${task.process}"), val('R'), eval("Rscript -e 'cat(as.character(getRversion()))' 2>/dev/null || printf stub"), emit: versions_r, topic: versions
    tuple val("${task.process}"), val('data.table'), eval("Rscript -e 'cat(as.character(packageVersion(\"data.table\")))' 2>/dev/null || printf stub"), emit: versions_data_table, topic: versions
    script:
    relatedness = eda_tables.find { table -> table.name.endsWith('.relatedness.tsv') }
    heterozygosity = eda_tables.find { table -> table.name.endsWith('.heterozygosity.tsv') }
    sex_check = eda_tables.find { table -> table.name == "${meta.cohort}.sex_check.tsv" }
    qc_checks = eda_tables.find { table -> table.name.endsWith('.genotype_eda_checks.tsv') }
    if (!relatedness) error "The genotype EDA output for ${meta.cohort} has no relatedness table."
    if (!heterozygosity) error "The genotype EDA output for ${meta.cohort} has no heterozygosity table."
    if (!sex_check) error "The genotype EDA output for ${meta.cohort} has no sex-check table."
    if (!qc_checks) error "The sample QC output for ${meta.cohort} has no diagnostic-status table."
    """
    Rscript ${decision_script} \
        --sample-decisions '${sample_decisions}' \
        --relatedness '${relatedness}' \
        --heterozygosity '${heterozygosity}' \
        --sex-check '${sex_check}' \
        --qc-checks '${qc_checks}' \
        --ancestry '${target_ancestry}' \
        --sample-missingness '${sample_missingness}' \
        --heterozygosity-z-threshold '${heterozygosity_z_threshold}' \
        --output '${meta.cohort}.participant_decisions.tsv' \
        --keep '${meta.cohort}.score_eligible.keep'
    """

    stub:
    """
    printf 'cohort\tFID\tIID\tmissingness\tretained_after_qc\tsample_missingness_flag\tsample_missingness_pass\tsample_missingness_threshold\theterozygosity_z\theterozygosity_flag\theterozygosity_pass\theterozygosity_z_threshold\tsex_check_pass\tsex_check_status\ttechnical_pass\tscore_eligible\trelated_flag\tancestry_flag\tancestry_distance\tprimary_analysis\tqc_status\treason\n${meta.cohort}\tTEST01\tTEST01\t0\tTRUE\tFALSE\tTRUE\t${sample_missingness}\t0\tFALSE\tTRUE\t${heterozygosity_z_threshold}\tTRUE\tPASS\tTRUE\tTRUE\tFALSE\tPASS\t0.5\tTRUE\tPASS\tEligible\n${meta.cohort}\tTEST02\tTEST02\t0\tTRUE\tFALSE\tTRUE\t${sample_missingness}\t0\tFALSE\tTRUE\t${heterozygosity_z_threshold}\tTRUE\tPASS\tTRUE\tTRUE\tFALSE\tPASS\t0.6\tTRUE\tPASS\tEligible\n' > ${meta.cohort}.participant_decisions.tsv
    printf 'TEST01\tTEST01\nTEST02\tTEST02\n' > ${meta.cohort}.score_eligible.keep
    """
}
