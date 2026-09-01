process REFERENCE_PLAN {
    tag reference_mode
    label 'process_single'
    label 'process_r'

    container 'ghcr.io/paulyrp/dnaprs-analysis:1.0.0'

    input:
    path provided_references
    path catalogue
    val methods
    val target_imputation
    val reference_mode
    path plan_script

    output:
    path 'reference_assets.tsv', emit: assets
    path 'reference_plan.tsv', emit: plan
    path 'versions.yml', emit: versions

    script:
    method_arg = methods.join(',')
    """
    Rscript ${plan_script} \
        --provided '${provided_references}' \
        --catalogue '${catalogue}' \
        --methods '${method_arg}' \
        --target-imputation '${target_imputation}' \
        --mode '${reference_mode}'
    printf '"${task.process}":\n  R: %s\n' "\$(Rscript -e 'cat(as.character(getRversion()))')" > versions.yml
    """

    stub:
    """
    head -n 1 ${catalogue} > reference_assets.tsv
    grep -E '^beagle_jar\t|^unbref3_jar\t' ${catalogue} >> reference_assets.tsv || true
    printf 'reference_type\tsource\tstatus\nreferences\tstub\tPLANNED\n' > reference_plan.tsv
    printf '"${task.process}":\n  R: stub\n' > versions.yml
    """
}
