process PARSE_PLINK_SCORE {
    tag "${target.cohort}:${gwas.trait_id}"
    label 'process_single'

    container 'ghcr.io/paulyrp/dnaprs-analysis:1.0.0'

    input:
    tuple val(target), val(gwas), path(score), path(used), path(score_log), path(weight), path(target_pvar)
    path parser_script
    path audit_script

    output:
    tuple val(target), val(gwas), path("${target.cohort}.${gwas.trait_id}.${target.score_method ?: 'plink_ct'}.score.tsv"), path("${target.cohort}.${gwas.trait_id}.${target.score_method ?: 'plink_ct'}.score_qc.tsv"), emit: scores
    tuple val(target), val(gwas), path("${target.cohort}.${gwas.trait_id}.${target.score_method ?: 'plink_ct'}.variant_compatibility.tsv"), path("${target.cohort}.${gwas.trait_id}.${target.score_method ?: 'plink_ct'}.variant_coverage.tsv"), emit: compatibility
    tuple val(target), val(gwas), path(score_log), emit: logs
    tuple val("${task.process}"), val('R'), eval("Rscript -e 'cat(as.character(getRversion()))' 2>/dev/null || printf stub"), emit: versions_r, topic: versions
    tuple val("${task.process}"), val('data.table'), eval("Rscript -e 'cat(as.character(packageVersion(\"data.table\")))' 2>/dev/null || printf stub"), emit: versions_data_table, topic: versions
    script:
    """
    Rscript ${parser_script} \
        --score '${score}' \
        --used '${used}' \
        --cohort '${target.cohort}' \
        --role '${target.role}' \
        --trait-id '${gwas.trait_id}' \
        --prs-name '${gwas.prs_name}' \
        --method '${target.score_method ?: 'plink_ct'}'
    Rscript ${audit_script} \
        --weights '${weight}' \
        --pvar '${target_pvar}' \
        --used '${used}' \
        --cohort '${target.cohort}' \
        --trait-id '${gwas.trait_id}' \
        --method '${target.score_method ?: 'plink_ct'}' \
        --scoring-stage '${target.scoring_stage ?: 'direct'}'
    """

    stub:
    """
    printf 'cohort\trole\ttrait_id\tprs_name\tmethod\tFID\tIID\traw_prs\tallele_count\tused_variants\n${target.cohort}\t${target.role}\t${gwas.trait_id}\t${gwas.prs_name}\t${target.score_method ?: 'plink_ct'}\tTEST01\tTEST01\t0.10\t2\t1\n${target.cohort}\t${target.role}\t${gwas.trait_id}\t${gwas.prs_name}\t${target.score_method ?: 'plink_ct'}\tTEST02\tTEST02\t0.20\t2\t1\n' > ${target.cohort}.${gwas.trait_id}.${target.score_method ?: 'plink_ct'}.score.tsv
    printf 'cohort\trole\ttrait_id\tprs_name\tmethod\tparticipants\tused_variants\n${target.cohort}\t${target.role}\t${gwas.trait_id}\t${gwas.prs_name}\t${target.score_method ?: 'plink_ct'}\t2\t1\n' > ${target.cohort}.${gwas.trait_id}.${target.score_method ?: 'plink_ct'}.score_qc.tsv
    printf 'cohort\ttrait_id\tmethod\tscoring_stage\tvariant_id\teffect_allele\ttarget_ref\ttarget_alt\tallele_state\tused\treason\n${target.cohort}\t${gwas.trait_id}\t${target.score_method ?: 'plink_ct'}\t${target.scoring_stage ?: 'direct'}\t1:100:A:G\tA\tA\tG\tDIRECT\tTRUE\tUsed by PLINK scoring\n' > ${target.cohort}.${gwas.trait_id}.${target.score_method ?: 'plink_ct'}.variant_compatibility.tsv
    printf 'cohort\ttrait_id\tmethod\tscoring_stage\trequested_variants\ttarget_compatible_variants\tused_variants\tused_fraction\treview_required\tstatus\n${target.cohort}\t${gwas.trait_id}\t${target.score_method ?: 'plink_ct'}\t${target.scoring_stage ?: 'direct'}\t1\t1\t1\t1\tFALSE\tPASS\n' > ${target.cohort}.${gwas.trait_id}.${target.score_method ?: 'plink_ct'}.variant_coverage.tsv
    """
}
