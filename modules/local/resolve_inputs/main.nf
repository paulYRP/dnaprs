process RESOLVE_INPUTS {
    tag 'raw inputs'
    label 'process_single'

    container 'ghcr.io/paulyrp/dnaprs-analysis:1.0.0'

    input:
    val input_spec
    val gwas_spec
    val reference_spec
    val models_spec
    val genome_build
    val methods
    val target_imputation
    val stop_after
    val reference_only
    val reference_mode
    val reference_bundle
    val phenotype
    val outcome
    val covariates
    val model_type
    val participant_id
    val group_column
    val control_value
    val case_value
    val beagle_jar
    val unbref3_jar
    path target_source, stageAs: 'target_source/*'
    path gwas_source, stageAs: 'gwas_source/*'
    path reference_source, stageAs: 'reference_source/*'
    path resolver_script

    output:
    path 'targets.tsv', emit: targets
    path 'gwas.tsv', emit: gwas
    path 'references.tsv', emit: references
    path 'models.tsv', emit: models
    path 'run_plan.tsv', emit: run_plan
    path 'reference_settings.yml', emit: settings
    path 'input_resolution.tsv', emit: checks
    tuple val("${task.process}"), val('R'), eval("Rscript -e 'cat(as.character(getRversion()))' 2>/dev/null || printf stub"), emit: versions_r, topic: versions
    tuple val("${task.process}"), val('jsonlite'), eval("Rscript -e 'cat(as.character(packageVersion(\"jsonlite\")))' 2>/dev/null || printf stub"), emit: versions_jsonlite, topic: versions
    tuple val("${task.process}"), val('openssl'), eval("Rscript -e 'cat(as.character(packageVersion(\"openssl\")))' 2>/dev/null || printf stub"), emit: versions_openssl, topic: versions
    tuple val("${task.process}"), val('yaml'), eval("Rscript -e 'cat(as.character(packageVersion(\"yaml\")))' 2>/dev/null || printf stub"), emit: versions_yaml, topic: versions
    script:
    method_arg = methods.join(',')
    """
    Rscript ${resolver_script} \
        --input-spec '${input_spec}' \
        --gwas-spec '${gwas_spec}' \
        --reference-spec '${reference_spec}' \
        --models-spec '${models_spec}' \
        --genome-build '${genome_build}' \
        --methods '${method_arg}' \
        --target-imputation '${target_imputation}' \
        --stop-after '${stop_after}' \
        --reference-only '${reference_only}' \
        --reference-mode '${reference_mode}' \
        --reference-bundle '${reference_bundle}' \
        --phenotype '${phenotype}' \
        --outcome '${outcome}' \
        --covariates '${covariates}' \
        --model-type '${model_type}' \
        --participant-id '${participant_id}' \
        --group-column '${group_column}' \
        --control-value '${control_value}' \
        --case-value '${case_value}' \
        --beagle-jar '${beagle_jar}' \
        --unbref3-jar '${unbref3_jar}' \
        --target-staged ${target_source} \
        --gwas-staged ${gwas_source} \
        --reference-staged ${reference_source}
    """

    stub:
    """
    awk -F '\t' -v OFS='\t' -v base='${projectDir}' '
        NR > 1 {
            if (substr(\$4, 1, 1) != "/") \$4 = base "/" \$4
            if (\$5 != "" && substr(\$5, 1, 1) != "/") \$5 = base "/" \$5
            if (\$6 != "" && substr(\$6, 1, 1) != "/") \$6 = base "/" \$6
        }
        { print }
    ' ${projectDir}/tests/data/target_manifest.tsv > targets.tsv
    awk -F '\t' -v OFS='\t' -v base='${projectDir}' '
        NR > 1 && substr(\$3, 1, 1) != "/" { \$3 = base "/" \$3 }
        { print }
    ' ${projectDir}/tests/data/gwas_manifest.tsv > gwas.tsv
    awk -F '\t' -v OFS='\t' -v base='${projectDir}/tests/data' '
        NR > 1 {
            if (substr(\$5, 1, 1) != "/") \$5 = base "/" \$5
            if (\$6 != "" && substr(\$6, 1, 1) != "/") \$6 = base "/" \$6
        }
        { print }
    ' ${projectDir}/tests/data/reference_manifest.tsv > references.tsv
    if [[ -n '${phenotype}' ]]; then
        cp ${projectDir}/tests/data/phenotype_models.tsv models.tsv
    else
        head -n 1 ${projectDir}/tests/data/phenotype_models.tsv > models.tsv
    fi
    printf 'setting\tvalue\nstop_after\t${stop_after}\nrun_imputation\t${target_imputation}\nrun_prs\ttrue\nrun_phenotype\ttrue\nreference_only\t${reference_only}\nrequired_reference_roles\tdbsnp,reference_fasta,imputation_panel,population_panel,related_samples,unbref3_jar,genetic_map,beagle_jar,sbayesrc_ld_source,annotation_source\n' > run_plan.tsv
    printf 'genome: ${genome_build}\nmethods:\n  - plink_ct\nreference_only: ${reference_only}\n' > reference_settings.yml
    printf 'kind\tid\tsource\tresolution\ntarget\tTEST\tstub\tstub\n' > input_resolution.tsv
    """
}
