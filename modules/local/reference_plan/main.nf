process REFERENCE_PLAN {
    tag "${reference_mode}"
    label 'process_single'

    container 'ghcr.io/paulyrp/dnaprs-analysis:1.0.0'

    input:
    path provided_references
    path catalogue
    path run_plan
    val reference_mode
    path plan_script

    output:
    path 'reference_assets.tsv', emit: assets
    path 'reference_plan.tsv', emit: plan
    tuple val("${task.process}"), val('R'), eval("Rscript -e 'cat(as.character(getRversion()))' 2>/dev/null || printf stub"), emit: versions_r, topic: versions
    script:
    """
    Rscript ${plan_script} \
        --provided '${provided_references}' \
        --catalogue '${catalogue}' \
        --run-plan '${run_plan}' \
        --mode '${reference_mode}'
    """

    stub:
    """
    head -n 1 ${catalogue} > reference_assets.tsv
    grep -E '^beagle_jar\t|^unbref3_jar\t' ${catalogue} >> reference_assets.tsv || true
    printf 'reference_type\tsource\tstatus\nreferences\tstub\tPLANNED\n' > reference_plan.tsv
    """
}
