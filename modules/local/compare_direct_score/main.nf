process COMPARE_DIRECT_SCORE {
    tag "${target.cohort}:${gwas.trait_id}"
    label 'process_single'

    container 'ghcr.io/paulyrp/dnaprs-analysis:1.0.1'

    input:
    tuple val(target), val(gwas), path(primary_score), path(primary_used), path(direct_score), path(direct_used)
    path comparison_script

    output:
    tuple val(target), val(gwas), path("${target.cohort}.${gwas.trait_id}.plink_ct.sensitivity.tsv"), path("${target.cohort}.${gwas.trait_id}.plink_ct.sensitivity_qc.tsv"), emit: comparison
    tuple val("${task.process}"), val('R'), eval("Rscript -e 'cat(as.character(getRversion()))' 2>/dev/null || printf stub"), emit: versions_r, topic: versions
    tuple val("${task.process}"), val('data.table'), eval("Rscript -e 'cat(as.character(packageVersion(\"data.table\")))' 2>/dev/null || printf stub"), emit: versions_data_table, topic: versions
    script:
    """
    Rscript ${comparison_script} \
        --primary-score '${primary_score}' \
        --primary-used '${primary_used}' \
        --direct-score '${direct_score}' \
        --direct-used '${direct_used}'
    """

    stub:
    """
    printf 'cohort\ttrait_id\tprs_name\tFID\tIID\traw_prs_imputed\traw_prs_direct\timputed_prs_z\tdirect_prs_z\tz_difference\n${target.cohort}\t${gwas.trait_id}\t${gwas.prs_name}\tTEST01\tTEST01\t0.10\t0.08\t-0.7071\t-0.7071\t0\n${target.cohort}\t${gwas.trait_id}\t${gwas.prs_name}\tTEST02\tTEST02\t0.20\t0.18\t0.7071\t0.7071\t0\n' > ${target.cohort}.${gwas.trait_id}.plink_ct.sensitivity.tsv
    printf 'cohort\ttrait_id\tprs_name\tparticipants\timputed_scoring_variants\tdirect_scoring_variants\tshared_variants\timputed_only_variants\tdirect_only_variants\tpearson_r\tspearman_r\tmean_z_difference\tstatus\n${target.cohort}\t${gwas.trait_id}\t${gwas.prs_name}\t2\t1\t1\t1\t0\t0\t1\t1\t0\tPASS\n' > ${target.cohort}.${gwas.trait_id}.plink_ct.sensitivity_qc.tsv
    """
}
