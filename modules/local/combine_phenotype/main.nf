process COMBINE_PHENOTYPE {
    tag 'phenotype results'
    label 'process_low'
    label 'process_r'

    container 'ghcr.io/paulyrp/dnaprs-analysis:1.0.0'

    input:
    path associations
    path fitted_models
    path plot_data
    path permutations
    path influence
    path phenotype_prs
    path combine_script

    output:
    path 'phenotype_associations.tsv', emit: associations
    path 'phenotype_models_fitted.tsv', emit: fitted_models
    path 'phenotype_plot_data.tsv', emit: plot_data
    path 'phenotype_permutations.tsv', emit: permutations
    path 'phenotype_influence.tsv', emit: influence
    path 'phenoPRS.csv', emit: phenotype_prs
    path 'versions.yml', emit: versions

    script:
    """
    Rscript ${combine_script} --input-dir .

    cat > versions.yml <<-VERSIONS
    "${task.process}":
        R: \$(Rscript -e 'cat(as.character(getRversion()))')
        data.table: \$(Rscript -e "cat(as.character(packageVersion('data.table')))")
    VERSIONS
    """

    stub:
    """
    cp ${associations[0]} phenotype_associations.tsv
    cp ${fitted_models[0]} phenotype_models_fitted.tsv
    cp ${plot_data[0]} phenotype_plot_data.tsv
    cp ${permutations[0]} phenotype_permutations.tsv
    cp ${influence[0]} phenotype_influence.tsv
    cp ${phenotype_prs[0]} phenoPRS.csv
    printf '"${task.process}":\n  R: stub\n  data.table: stub\n' > versions.yml
    """
}
