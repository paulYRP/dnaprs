process PARTICIPANT_DECISIONS {
    tag "${meta.cohort}"
    label 'process_low'
    label 'process_r'

    container 'ghcr.io/paulyrp/dnaprs-analysis:1.0.0'

    input:
    tuple val(meta), path(sample_decisions), path(eda_tables), path(target_ancestry)
    path decision_script

    output:
    tuple val(meta), path("${meta.cohort}.participant_decisions.tsv"), path("${meta.cohort}.score_eligible.keep"), emit: decisions
    path 'versions.yml', emit: versions

    script:
    relatedness = eda_tables.find { it.name == "${meta.cohort}.relatedness.tsv" }
    if (!relatedness) error "The genotype EDA output for ${meta.cohort} has no relatedness table."
    """
    Rscript ${decision_script} \
        --sample-decisions '${sample_decisions}' \
        --relatedness '${relatedness}' \
        --ancestry '${target_ancestry}' \
        --output '${meta.cohort}.participant_decisions.tsv' \
        --keep '${meta.cohort}.score_eligible.keep'

    cat > versions.yml <<-VERSIONS
    "${task.process}:${meta.cohort}":
        R: \$(Rscript -e 'cat(as.character(getRversion()))')
        data.table: \$(Rscript -e "cat(as.character(packageVersion('data.table')))" 2>/dev/null)
    VERSIONS
    """

    stub:
    """
    printf 'cohort\tFID\tIID\tmissingness\ttechnical_pass\tscore_eligible\trelated_flag\tancestry_flag\tancestry_distance\tprimary_analysis\treason\n${meta.cohort}\tTEST01\tTEST01\t0\tTRUE\tTRUE\tFALSE\tPASS\t0.5\tTRUE\tEligible\n${meta.cohort}\tTEST02\tTEST02\t0\tTRUE\tTRUE\tFALSE\tPASS\t0.6\tTRUE\tEligible\n' > ${meta.cohort}.participant_decisions.tsv
    printf 'TEST01\tTEST01\nTEST02\tTEST02\n' > ${meta.cohort}.score_eligible.keep
    printf '"${task.process}:${meta.cohort}":\n  R: stub\n  data.table: stub\n' > versions.yml
    """
}
