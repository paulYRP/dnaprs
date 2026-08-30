process COMBINE_SCORES {
    tag 'all scores'
    label 'process_medium'
    label 'process_r'

    input:
    path score_files
    path combine_script

    output:
    path 'prs_scores_long.tsv', emit: scores_long
    path 'prs_scores_wide.tsv', emit: scores_wide
    path 'score_qc.tsv', emit: score_qc
    path 'method_concordance.tsv', emit: concordance
    path 'versions.yml', emit: versions

    script:
    scoreARG = score_files.collect { score_file -> score_file.toString() }.join(',')
    """
    Rscript ${combine_script} --scores '${scoreARG}'

    cat > versions.yml <<-VERSIONS
    "${task.process}":
        R: \$(Rscript -e 'cat(as.character(getRversion()))')
        data.table: \$(Rscript -e "cat(as.character(packageVersion('data.table')))" 2>/dev/null)
    VERSIONS
    """

    stub:
    """
    cp \$(find . -maxdepth 1 -name '*.score.tsv' | head -n 1) prs_scores_long.tsv
    printf 'cohort\trole\tFID\tIID\tMDD_PRS_PLINK_CT\nTEST\ttarget\tTEST01\tTEST01\t-0.7071\nTEST\ttarget\tTEST02\tTEST02\t0.7071\n' > prs_scores_wide.tsv
    printf 'cohort\trole\ttrait_id\tprs_name\tmethod\tparticipants\tfinite_scores\traw_mean\traw_sd\tz_mean\tz_sd\tused_variants\tstatus\nTEST\ttarget\tMDD\tMDD_PRS\tplink_ct\t2\t2\t0.15\t0.07\t0\t1\t1\tPASS\n' > score_qc.tsv
    printf 'cohort\ttrait_id\tprs_name\tmethod_1\tmethod_2\tparticipants\tpearson_r\tspearman_r\n' > method_concordance.tsv
    printf '"${task.process}":\n  R: stub\n  data.table: stub\n' > versions.yml
    """
}
