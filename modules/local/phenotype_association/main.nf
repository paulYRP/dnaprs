process PHENOTYPE_ASSOCIATION {
    tag "${score_job.cohort}:${score_job.trait_id}:${score_job.method}:${phenotype_model.baseName}"
    label 'process_medium'
    label 'process_r'

    container 'ghcr.io/paulyrp/dnaprs-analysis:1.0.0'

    input:
    tuple val(score_job), path(scores_long), path(phenotype_model)
    path phenotype_file
    path association_script
    val seed

    output:
    path '*.phenotype_associations.tsv', emit: associations
    path '*.phenotype_models_fitted.tsv', emit: fitted_models
    path '*.phenotype_plot_data.tsv', emit: plot_data
    path '*.phenotype_permutations.tsv', emit: permutations
    path '*.phenotype_influence.tsv', emit: influence
    path '*.phenoPRS.csv', emit: phenotype_prs
    path '*.versions.yml', emit: versions

    script:
    model_id = phenotype_model.baseName
    job_id = "${model_id}__${score_job.cohort}__${score_job.trait_id}__${score_job.method}"
    """
    Rscript ${association_script} \
        --scores '${scores_long}' \
        --phenotype '${phenotype_file}' \
        --models '${phenotype_model}' \
        --cohort '${score_job.cohort}' \
        --trait-id '${score_job.trait_id}' \
        --method '${score_job.method}' \
        --seed '${seed}'

    mv phenotype_associations.tsv ${job_id}.phenotype_associations.tsv
    mv phenotype_models_fitted.tsv ${job_id}.phenotype_models_fitted.tsv
    mv phenotype_plot_data.tsv ${job_id}.phenotype_plot_data.tsv
    mv phenotype_permutations.tsv ${job_id}.phenotype_permutations.tsv
    mv phenotype_influence.tsv ${job_id}.phenotype_influence.tsv
    mv phenoPRS.csv ${job_id}.phenoPRS.csv

    cat > ${job_id}.versions.yml <<-VERSIONS
    "${task.process}:${job_id}":
        R: \$(Rscript -e 'cat(as.character(getRversion()))')
        glm2: \$(Rscript -e "if (requireNamespace('glm2', quietly=TRUE)) cat(as.character(packageVersion('glm2'))) else cat('not-used')" 2>/dev/null)
        lme4: \$(Rscript -e "if (requireNamespace('lme4', quietly=TRUE)) cat(as.character(packageVersion('lme4'))) else cat('not-used')" 2>/dev/null)
    VERSIONS
    """

    stub:
    model_id = phenotype_model.baseName
    job_id = "${model_id}__${score_job.cohort}__${score_job.trait_id}__${score_job.method}"
    """
    printf 'model_id\tcohort\trole\ttrait_id\tprs_name\tmethod\tfamily\tn\tbeta\tstd_error\tci_low\tci_high\tp_value\tnull_fit\tfull_fit\tincremental_fit\tfit_metric\tpermutation_scheme\tpermutations\tempirical_p\texpected_direction\tdirection_match\tstatus\n' > ${job_id}.phenotype_associations.tsv
    printf 'model_id\tcohort\tmethod\tformula\tnull_formula\testimator\texpected_direction\tprimary\n' > ${job_id}.phenotype_models_fitted.tsv
    printf 'model_id\toutcome\tcohort\trole\ttrait_id\tprs_name\tmethod\tfamily\testimator\tIID\tobserved\tfitted_null\tfitted_full\tresidual_full\tadjusted_outcome\tadjusted_prs\n' > ${job_id}.phenotype_plot_data.tsv
    printf 'model_id\tcohort\trole\ttrait_id\tprs_name\tmethod\tfamily\testimator\tpermutation_id\tpermuted_beta\tobserved_beta\tpermutation_scheme\tstatus\treason\n' > ${job_id}.phenotype_permutations.tsv
    printf 'model_id\tcohort\trole\ttrait_id\tprs_name\tmethod\tfamily\testimator\tIID\tfull_beta\tbeta_without\tbeta_change\tstatus\treason\n' > ${job_id}.phenotype_influence.tsv
    printf 'IID\tMDD_PRS\tMDD_PRS_SBAYESRC\nTEST01\t-0.7071\t-0.7071\nTEST02\t0.7071\t0.7071\n' > ${job_id}.phenoPRS.csv
    printf '"${task.process}:${job_id}":\n  R: stub\n  glm2: stub\n  lme4: stub\n' > ${job_id}.versions.yml
    """
}
