process PLINK_SCORE {
    tag "${target.cohort}:${gwas.trait_id}"
    label 'process_high'

    container 'ghcr.io/paulyrp/dnaprs-plink2:2.0.0-a.6.12'

    input:
    tuple val(target), path(target_dir), path(target_qc), path(participant_decisions), path(participant_keep), val(gwas), path(weight), path(weight_qc), path(harmonisation_qc), path(clump_log), val(reference), path(reference_freq), path(reference_log)
    output:
    tuple val(target), val(gwas), path("${target.cohort}.${gwas.trait_id}.${target.score_method ?: 'plink_ct'}.sscore"), path("${target.cohort}.${gwas.trait_id}.${target.score_method ?: 'plink_ct'}.sscore.vars"), path("${target.cohort}.${gwas.trait_id}.${target.score_method ?: 'plink_ct'}.log"), path("${target.cohort}.${gwas.trait_id}.${target.score_method ?: 'plink_ct'}.weights.tsv"), path("${target.cohort}.${gwas.trait_id}.${target.score_method ?: 'plink_ct'}.target.pvar"), emit: raw
    tuple val("${task.process}"), val('plink2'), eval("command -v plink2 >/dev/null && plink2 --version 2>&1 | head -n 1 | cut -d ' ' -f 2 | sed 's/^v//' || printf stub"), emit: versions_plink2, topic: versions
    script:
    """
    plink2 \
        --pfile '${target_dir}/${target.cohort}' \
        --read-freq '${reference_freq}' \
        --keep '${participant_keep}' \
        --score '${weight}' 1 2 3 header-read list-variants cols=fid,nallele,denom,dosagesum,scoreavgs,scoresums \
        --threads '${task.cpus}' \
        --out '${target.cohort}.${gwas.trait_id}.${target.score_method ?: 'plink_ct'}'
    cp '${weight}' '${target.cohort}.${gwas.trait_id}.${target.score_method ?: 'plink_ct'}.weights.tsv'
    cp '${target_dir}/${target.cohort}.pvar' '${target.cohort}.${gwas.trait_id}.${target.score_method ?: 'plink_ct'}.target.pvar'
    """

    stub:
    """
    printf '#FID\tIID\tALLELE_CT\tNAMED_ALLELE_DOSAGE_SUM\tSCORE1_AVG\tSCORE1_SUM\nTEST01\tTEST01\t2\t1\t0.05\t0.10\nTEST02\tTEST02\t2\t2\t0.10\t0.20\n' > ${target.cohort}.${gwas.trait_id}.${target.score_method ?: 'plink_ct'}.sscore
    printf '1:100:A:G\n' > ${target.cohort}.${gwas.trait_id}.${target.score_method ?: 'plink_ct'}.sscore.vars
    printf 'PLINK score stub\n' > ${target.cohort}.${gwas.trait_id}.${target.score_method ?: 'plink_ct'}.log
    cp ${weight} ${target.cohort}.${gwas.trait_id}.${target.score_method ?: 'plink_ct'}.weights.tsv
    cp ${target_dir}/${target.cohort}.pvar ${target.cohort}.${gwas.trait_id}.${target.score_method ?: 'plink_ct'}.target.pvar
    """
}
