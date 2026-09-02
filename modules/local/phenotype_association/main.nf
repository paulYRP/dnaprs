process PHENOTYPE_ASSOCIATION {
    tag "${score_job.cohort}:${score_job.trait_id}:${score_job.method}:${model_id}"
    label 'process_single'

    container 'ghcr.io/paulyrp/dnaprs-analysis:1.0.0'

    input:
    tuple val(score_job), path(scores_long), path(phenotype_models), val(model_id)
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
    tuple val("${task.process}"), val('R'), eval("Rscript -e 'cat(as.character(getRversion()))' 2>/dev/null || printf stub"), emit: versions_r, topic: versions
    tuple val("${task.process}"), val('glm2'), eval("Rscript -e 'if (requireNamespace(\"glm2\", quietly=TRUE)) cat(as.character(packageVersion(\"glm2\"))) else cat(\"not-used\")' 2>/dev/null || printf stub"), emit: versions_glm2, topic: versions
    tuple val("${task.process}"), val('lme4'), eval("Rscript -e 'if (requireNamespace(\"lme4\", quietly=TRUE)) cat(as.character(packageVersion(\"lme4\"))) else cat(\"not-used\")' 2>/dev/null || printf stub"), emit: versions_lme4, topic: versions
    script:
    job_id = "${model_id}__${score_job.cohort}__${score_job.trait_id}__${score_job.method}"
    """
    Rscript ${association_script} \
        --scores '${scores_long}' \
        --phenotype '${phenotype_file}' \
        --models '${phenotype_models}' \
        --model-id '${model_id}' \
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
    """

    stub:
    job_id = "${model_id}__${score_job.cohort}__${score_job.trait_id}__${score_job.method}"
    """
    printf 'model_id\tcohort\trole\ttrait_id\tprs_name\tmethod\tfamily\tn\tbeta\tstd_error\tci_low\tci_high\tp_value\tnull_fit\tfull_fit\tincremental_fit\tfit_metric\tpermutation_scheme\tpermutations\tempirical_p\texpected_direction\tdirection_match\tinput_rows\tcomplete_cases\texcluded_missing\tconvergence\tsingular\tseparation\tdiagnostics\tstatus\n' > ${job_id}.phenotype_associations.tsv
    printf 'model_id\tcohort\tmethod\tformula\tnull_formula\testimator\texpected_direction\tprimary\tconvergence\tsingular\tseparation\tdiagnostics\tstatus\n' > ${job_id}.phenotype_models_fitted.tsv
    printf 'model_id\toutcome\tcohort\trole\ttrait_id\tprs_name\tmethod\tfamily\testimator\tIID\tobserved\tfitted_null\tfitted_full\tresidual_full\tadjusted_outcome\tadjusted_prs\n' > ${job_id}.phenotype_plot_data.tsv
    printf 'model_id\tcohort\trole\ttrait_id\tprs_name\tmethod\tfamily\testimator\tpermutation_id\tpermuted_beta\tobserved_beta\tpermutation_scheme\tstatus\treason\n' > ${job_id}.phenotype_permutations.tsv
    printf 'model_id\tcohort\trole\ttrait_id\tprs_name\tmethod\tfamily\testimator\tIID\tfull_beta\tbeta_without\tbeta_change\tstatus\treason\n' > ${job_id}.phenotype_influence.tsv
    printf 'IID\tMDD_PRS\tMDD_PRS_SBAYESRC\nTEST01\t-0.7071\t-0.7071\nTEST02\t0.7071\t0.7071\n' > ${job_id}.phenoPRS.csv
    """
}
