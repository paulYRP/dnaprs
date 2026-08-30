process PLINK_SCORE {
    tag "${target.cohort}:${gwas.trait_id}"
    label 'process_high'
    label 'process_plink'
    label 'process_r'

    input:
    tuple val(target), path(target_dir), path(target_qc), val(gwas), path(weight), path(weight_qc), path(harmonisation_qc), path(clump_log), val(reference), path(reference_freq), path(reference_log)
    path parser_script

    output:
    tuple val(target), val(gwas), path("${target.cohort}.${gwas.trait_id}.plink_ct.score.tsv"), path("${target.cohort}.${gwas.trait_id}.plink_ct.score_qc.tsv"), emit: scores
    tuple val(target), val(gwas), path("${target.cohort}.${gwas.trait_id}.plink_ct.log"), emit: logs
    path 'versions.yml', emit: versions

    script:
    """
    plink2 \
        --pfile '${target_dir}/${target.cohort}' \
        --read-freq '${reference_freq}' \
        --score '${weight}' 1 2 3 header-read list-variants cols=fid,nallele,denom,dosagesum,scoreavgs,scoresums \
        --threads '${task.cpus}' \
        --out '${target.cohort}.${gwas.trait_id}.plink_ct'

    Rscript ${parser_script} \
        --score '${target.cohort}.${gwas.trait_id}.plink_ct.sscore' \
        --used '${target.cohort}.${gwas.trait_id}.plink_ct.sscore.vars' \
        --cohort '${target.cohort}' \
        --role '${target.role}' \
        --trait-id '${gwas.trait_id}' \
        --prs-name '${gwas.prs_name}'

    cat > versions.yml <<-VERSIONS
    "${task.process}:${target.cohort}:${gwas.trait_id}":
        plink2: \$(plink2 --version 2>&1 | head -n 1 | cut -d ' ' -f 2 | sed 's/^v//')
        R: \$(Rscript -e 'cat(as.character(getRversion()))')
    VERSIONS
    """

    stub:
    """
    printf 'cohort\trole\ttrait_id\tprs_name\tmethod\tFID\tIID\traw_prs\tallele_count\tused_variants\n${target.cohort}\t${target.role}\t${gwas.trait_id}\t${gwas.prs_name}\tplink_ct\tTEST01\tTEST01\t0.10\t2\t1\n${target.cohort}\t${target.role}\t${gwas.trait_id}\t${gwas.prs_name}\tplink_ct\tTEST02\tTEST02\t0.20\t2\t1\n' > ${target.cohort}.${gwas.trait_id}.plink_ct.score.tsv
    printf 'cohort\trole\ttrait_id\tprs_name\tmethod\tparticipants\tused_variants\n${target.cohort}\t${target.role}\t${gwas.trait_id}\t${gwas.prs_name}\tplink_ct\t2\t1\n' > ${target.cohort}.${gwas.trait_id}.plink_ct.score_qc.tsv
    printf 'PLINK score stub\n' > ${target.cohort}.${gwas.trait_id}.plink_ct.log
    printf '"${task.process}:${target.cohort}:${gwas.trait_id}":\n  plink2: stub\n  R: stub\n' > versions.yml
    """
}
