process REFERENCE_ANCESTRY {
    tag "${meta.cohort}"
    label 'process_high'

    container 'docker.io/zhiliz/sbayesrc:0.2.6@sha256:5a6139e6c3ab471799c059fe7032870f22bacc07dfadf14fc0517314db5ea4bb'

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
    path 'versions.yml', emit: versions

    script:
    """
    bash ${ancestry_script} \
        '${meta.cohort}' '${target_dir}' '${reference_dir}' \
        '${ancestry_pcs}' '${ancestry_percentile}' '${task.cpus}' '${classifier_script}'

    cat > versions.yml <<-VERSIONS
    "${task.process}:${meta.cohort}":
        plink2: \$(plink2 --version 2>&1 | head -n 1 | cut -d ' ' -f 2 | sed 's/^v//')
        R: \$(Rscript -e 'cat(as.character(getRversion()))')
        data.table: \$(Rscript -e "cat(as.character(packageVersion('data.table')))" 2>/dev/null)
    VERSIONS
    """

    stub:
    """
    printf 'cohort\tFID\tIID\tPC1\tPC2\tancestry_distance\tancestry_threshold\tancestry_percentile\tancestry_flag\n${meta.cohort}\tTEST01\tTEST01\t0.01\t0.02\t0.5\t3\t${ancestry_percentile}\tPASS\n${meta.cohort}\tTEST02\tTEST02\t-0.01\t-0.02\t0.6\t3\t${ancestry_percentile}\tPASS\n' > ${meta.cohort}.target_ancestry.tsv
    printf 'IID\tFID\tPC1\tPC2\tpopulation\tsuper_population\tancestry_distance\tancestry_flag\nTEST01\tTEST01\t0.01\t0.02\tGBR\tEUR\t0.5\tEUR_REFERENCE\n' > ${meta.cohort}.reference_projection.tsv
    printf 'cohort\tmatched_variants\tpruned_variants\treference_participants\teuropean_reference_participants\ttarget_participants\tpcs\tpercentile\tdistance_threshold\tancestry_pass\tancestry_outlier\tstatus\n${meta.cohort}\t1\t1\t2\t2\t2\t2\t${ancestry_percentile}\t3\t2\t0\tPASS\n' > ${meta.cohort}.ancestry_summary.tsv
    printf 'Reference ancestry stub\n' > ${meta.cohort}.reference_ancestry.log
    printf '"${task.process}:${meta.cohort}":\n  plink2: stub\n  R: stub\n  data.table: stub\n' > versions.yml
    """
}
