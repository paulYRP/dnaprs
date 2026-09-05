process COMBINE_SCORES {
    tag 'all scores'
    label 'process_single'

    container 'ghcr.io/paulyrp/dnaprs-analysis:1.0.1'

    input:
    path score_files
    path participant_decisions
    path combine_script

    output:
    path 'prs_scores_long.tsv', emit: scores_long
    path 'prs_scores_wide.tsv', emit: scores_wide
    path 'score_qc.tsv', emit: score_qc
    path 'method_concordance.tsv', emit: concordance
    tuple val("${task.process}"), val('R'), eval("Rscript -e 'cat(as.character(getRversion()))' 2>/dev/null || printf stub"), emit: versions_r, topic: versions
    tuple val("${task.process}"), val('data.table'), eval("Rscript -e 'cat(as.character(packageVersion(\"data.table\")))' 2>/dev/null || printf stub"), emit: versions_data_table, topic: versions
    script:
    scoreARG = score_files.collect { score_file -> score_file.toString() }.join(',')
    decisionARG = participant_decisions.collect { decision_file -> decision_file.toString() }.join(',')
    """
    Rscript ${combine_script} --scores '${scoreARG}' --participant-decisions '${decisionARG}'
    """

    stub:
    """
    score_files=( ./*.score.tsv )
    awk 'BEGIN{FS=OFS="\t"} NR==1{print \$0,"sample_missingness_pass","heterozygosity_z","heterozygosity_pass","sex_check_pass","technical_pass","score_eligible","related_flag","ancestry_flag","ancestry_distance","primary_analysis"; next} {print \$0,"TRUE","0","TRUE","TRUE","TRUE","TRUE","FALSE","NOT_ASSESSED","NA","TRUE"}' "\${score_files[0]}" > prs_scores_long.tsv
    printf 'cohort\trole\tFID\tIID\tMDD_PRS\tsample_missingness_pass\theterozygosity_z\theterozygosity_pass\tsex_check_pass\ttechnical_pass\tscore_eligible\trelated_flag\tancestry_flag\tancestry_distance\tprimary_analysis\nTEST\ttarget\tTEST01\tTEST01\t-0.7071\tTRUE\t0\tTRUE\tTRUE\tTRUE\tTRUE\tFALSE\tNOT_ASSESSED\tNA\tTRUE\nTEST\ttarget\tTEST02\tTEST02\t0.7071\tTRUE\t0\tTRUE\tTRUE\tTRUE\tTRUE\tFALSE\tNOT_ASSESSED\tNA\tTRUE\n' > prs_scores_wide.tsv
    printf 'cohort\trole\ttrait_id\tprs_name\tmethod\tparticipants\tfinite_scores\traw_mean\traw_sd\tz_mean\tz_sd\tused_variants\tstatus\nTEST\ttarget\tMDD\tMDD_PRS\tplink_ct\t2\t2\t0.15\t0.07\t0\t1\t1\tPASS\n' > score_qc.tsv
    printf 'cohort\ttrait_id\tprs_name\tmethod_1\tmethod_2\tparticipants\tpearson_r\tspearman_r\n' > method_concordance.tsv
    """
}
