process REFERENCE_ANCESTRY {
    tag "${meta.cohort}"
    label 'process_high'

    container 'ghcr.io/paulyrp/dnaprs-imputation:1.1.0'

    input:
    tuple val(meta), path(target_dir), path(target_qc), val(reference), path(reference_dir)
    path ancestry_script
    path classifier_script
    val ancestry_pcs
    val ancestry_percentile

    output:
    tuple val(meta), path("${meta.cohort}.target_ancestry.tsv"), emit: target
    tuple val(meta), path("${meta.cohort}.reference_projection.tsv"), emit: reference
    tuple val(meta), path("${meta.cohort}.ancestry_summary.tsv"), emit: summary
    tuple val(meta), path("${meta.cohort}.reference_ancestry.log"), emit: logs
    tuple val("${task.process}"), val('plink2'), eval("command -v plink2 >/dev/null && plink2 --version 2>&1 | head -n 1 | cut -d ' ' -f 2 | sed 's/^v//' || printf stub"), emit: versions_plink2, topic: versions
    tuple val("${task.process}"), val('R'), eval("Rscript -e 'cat(as.character(getRversion()))' 2>/dev/null || printf stub"), emit: versions_r, topic: versions
    tuple val("${task.process}"), val('data.table'), eval("Rscript -e 'cat(as.character(packageVersion(\"data.table\")))' 2>/dev/null || printf stub"), emit: versions_data_table, topic: versions
    script:
    """
    bash ${ancestry_script} \
        '${meta.cohort}' '${target_dir}' '${reference_dir}' \
        '${ancestry_pcs}' '${ancestry_percentile}' '${task.cpus}' '${classifier_script}'
    """

    stub:
    """
    printf 'cohort\tFID\tIID\tPC1\tPC2\tancestry_distance\tancestry_threshold\tancestry_percentile\tancestry_flag\n${meta.cohort}\tTEST01\tTEST01\t0.01\t0.02\t0.5\t3\t${ancestry_percentile}\tPASS\n${meta.cohort}\tTEST02\tTEST02\t-0.01\t-0.02\t0.6\t3\t${ancestry_percentile}\tPASS\n' > ${meta.cohort}.target_ancestry.tsv
    printf 'IID\tFID\tPC1\tPC2\tpopulation\tsuper_population\tancestry_distance\tancestry_flag\nTEST01\tTEST01\t0.01\t0.02\tGBR\tEUR\t0.5\tEUR_REFERENCE\n' > ${meta.cohort}.reference_projection.tsv
    printf 'cohort\tmatched_variants\tpruned_variants\treference_participants\teuropean_reference_participants\ttarget_participants\tpcs\tpercentile\tdistance_threshold\tancestry_pass\tancestry_outlier\tstatus\n${meta.cohort}\t1\t1\t2\t2\t2\t2\t${ancestry_percentile}\t3\t2\t0\tPASS\n' > ${meta.cohort}.ancestry_summary.tsv
    printf 'Reference ancestry stub\n' > ${meta.cohort}.reference_ancestry.log
    """
}
