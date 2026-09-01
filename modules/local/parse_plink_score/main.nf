process PARSE_PLINK_SCORE {
    tag "${target.cohort}:${gwas.trait_id}"
    label 'process_low'
    label 'process_r'

    container 'ghcr.io/paulyrp/dnaprs-analysis:1.0.0'

    input:
    tuple val(target), val(gwas), path(score), path(used), path(score_log)
    path parser_script

    output:
    tuple val(target), val(gwas), path("${target.cohort}.${gwas.trait_id}.${target.score_method ?: 'plink_ct'}.score.tsv"), path("${target.cohort}.${gwas.trait_id}.${target.score_method ?: 'plink_ct'}.score_qc.tsv"), emit: scores
    tuple val(target), val(gwas), path(score_log), emit: logs
    path 'versions.yml', emit: versions

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

    cat > versions.yml <<-VERSIONS
    "${task.process}:${target.cohort}:${gwas.trait_id}":
        R: \$(Rscript -e 'cat(as.character(getRversion()))')
        data.table: \$(Rscript -e "cat(as.character(packageVersion('data.table')))" 2>/dev/null)
    VERSIONS
    """

    stub:
    """
    printf 'cohort\trole\ttrait_id\tprs_name\tmethod\tFID\tIID\traw_prs\tallele_count\tused_variants\n${target.cohort}\t${target.role}\t${gwas.trait_id}\t${gwas.prs_name}\t${target.score_method ?: 'plink_ct'}\tTEST01\tTEST01\t0.10\t2\t1\n${target.cohort}\t${target.role}\t${gwas.trait_id}\t${gwas.prs_name}\t${target.score_method ?: 'plink_ct'}\tTEST02\tTEST02\t0.20\t2\t1\n' > ${target.cohort}.${gwas.trait_id}.${target.score_method ?: 'plink_ct'}.score.tsv
    printf 'cohort\trole\ttrait_id\tprs_name\tmethod\tparticipants\tused_variants\n${target.cohort}\t${target.role}\t${gwas.trait_id}\t${gwas.prs_name}\t${target.score_method ?: 'plink_ct'}\t2\t1\n' > ${target.cohort}.${gwas.trait_id}.${target.score_method ?: 'plink_ct'}.score_qc.tsv
    printf '"${task.process}:${target.cohort}:${gwas.trait_id}":\n  R: stub\n  data.table: stub\n' > versions.yml
    """
}
