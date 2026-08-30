process PHENOTYPE_ASSOCIATION {
    tag 'prespecified models'
    label 'process_medium'
    label 'process_r'

    input:
    path scores_long
    path phenotype_file
    path phenotype_models
    path association_script

    output:
    path 'phenotype_associations.tsv', emit: associations
    path 'phenotype_models_fitted.tsv', emit: fitted_models
    path 'phenotype_plot_data.tsv', emit: plot_data
    path 'versions.yml', emit: versions

    script:
    """
    Rscript ${association_script} \
        --scores '${scores_long}' \
        --phenotype '${phenotype_file}' \
        --models '${phenotype_models}'

    cat > versions.yml <<-VERSIONS
    "${task.process}":
        R: \$(Rscript -e 'cat(as.character(getRversion()))')
        glm2: \$(Rscript -e "if (requireNamespace('glm2', quietly=TRUE)) cat(as.character(packageVersion('glm2'))) else cat('not-used')" 2>/dev/null)
        lme4: \$(Rscript -e "if (requireNamespace('lme4', quietly=TRUE)) cat(as.character(packageVersion('lme4'))) else cat('not-used')" 2>/dev/null)
    VERSIONS
    """

    stub:
    """
    printf 'model_id\tcohort\trole\ttrait_id\tprs_name\tmethod\tfamily\tn\tbeta\tstd_error\tci_low\tci_high\tp_value\tnull_fit\tfull_fit\tincremental_fit\tfit_metric\tstatus\n' > phenotype_associations.tsv
    printf 'model_id\tcohort\tmethod\tformula\tnull_formula\testimator\texpected_direction\tprimary\n' > phenotype_models_fitted.tsv
    printf 'model_id\toutcome\tcohort\trole\ttrait_id\tprs_name\tmethod\tfamily\testimator\tIID\tobserved\tfitted_null\tfitted_full\tresidual_full\tadjusted_outcome\tadjusted_prs\n' > phenotype_plot_data.tsv
    printf '"${task.process}":\n  R: stub\n  glm2: stub\n  lme4: stub\n' > versions.yml
    """
}
