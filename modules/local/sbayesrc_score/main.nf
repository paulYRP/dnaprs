process SBAYESRC_SCORE {
    tag "${target.cohort}:${gwas.trait_id}"
    label 'process_high'

    container 'docker.io/zhiliz/sbayesrc:0.2.6@sha256:5a6139e6c3ab471799c059fe7032870f22bacc07dfadf14fc0517314db5ea4bb'

    input:
    tuple val(target), path(target_dir), path(target_qc), path(participant_decisions), path(participant_keep), val(gwas), path(weight), path(parameter), path(model_qc), path(impute_qc), path(tidy_qc), path(harmonisation_qc)
    path sbayesrc_script

    output:
    tuple val(target), val(gwas), path("${target.cohort}.${gwas.trait_id}.sbayesrc.score.tsv"), path("${target.cohort}.${gwas.trait_id}.sbayesrc.score_qc.tsv"), emit: scores
    tuple val(target), val(gwas), path("${target.cohort}.${gwas.trait_id}.sbayesrc.score.log"), emit: logs
    tuple val("${task.process}"), val('SBayesRC'), eval("Rscript -e 'cat(as.character(packageVersion(\"SBayesRC\")))' 2>/dev/null || printf stub"), emit: versions_sbayesrc, topic: versions
    tuple val("${task.process}"), val('plink2'), eval("command -v plink2 >/dev/null && plink2 --version 2>&1 | head -n 1 | cut -d ' ' -f 2 | sed 's/^v//' || printf stub"), emit: versions_plink2, topic: versions
    tuple val("${task.process}"), val('R'), eval("Rscript -e 'cat(as.character(getRversion()))' 2>/dev/null || printf stub"), emit: versions_r, topic: versions
    script:
    """
    export OMP_NUM_THREADS='${task.cpus}'
    Rscript ${sbayesrc_script} \
        --action score \
        --input '${weight}' \
        --target-dir '${target_dir}' \
        --cohort '${target.cohort}' \
        --role '${target.role}' \
        --trait-id '${gwas.trait_id}' \
        --prs-name '${gwas.prs_name}' \
        --keep '${participant_keep}' \
        --plink \$(command -v plink2)
    """

    stub:
    """
    printf 'cohort\trole\ttrait_id\tprs_name\tmethod\tFID\tIID\traw_prs\tallele_count\tused_variants\n${target.cohort}\t${target.role}\t${gwas.trait_id}\t${gwas.prs_name}\tsbayesrc\tTEST01\tTEST01\t0.15\t2\tNA\n${target.cohort}\t${target.role}\t${gwas.trait_id}\t${gwas.prs_name}\tsbayesrc\tTEST02\tTEST02\t0.25\t2\tNA\n' > ${target.cohort}.${gwas.trait_id}.sbayesrc.score.tsv
    printf 'cohort\trole\ttrait_id\tprs_name\tmethod\tparticipants\trequested_variants\tused_variants\tused_fraction\treview_required\tfinite_scores\tstatus\n${target.cohort}\t${target.role}\t${gwas.trait_id}\t${gwas.prs_name}\tsbayesrc\t2\t1\t1\t1\tFALSE\t2\tPASS\n' > ${target.cohort}.${gwas.trait_id}.sbayesrc.score_qc.tsv
    printf 'SBayesRC score stub\n' > ${target.cohort}.${gwas.trait_id}.sbayesrc.score.log
    """
}
