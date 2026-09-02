process VALIDATE_MANIFESTS {
    tag "$run_name"
    label 'process_single'

    container 'ghcr.io/paulyrp/dnaprs-analysis:1.0.0'

    input:
    val run_name
    path target_manifest, stageAs: 'input_targets.tsv'
    path gwas_manifest, stageAs: 'input_gwas.tsv'
    path reference_manifest, stageAs: 'input_references.tsv'
    path phenotype_file
    path phenotype_models, stageAs: 'input_models.tsv'
    val methods
    val genome_build
    val seed
    val report_enabled
    val target_imputation
    val launch_dir
    val reference_base
    path run_plan, stageAs: 'input_run_plan.tsv'
    path validator_script

    output:
    path 'targets.tsv', emit: target_manifest
    path 'gwas.tsv', emit: gwas_manifest
    path 'references.tsv', emit: reference_manifest
    path 'models.tsv', emit: phenotype_models
    path 'run_plan.tsv', emit: run_plan
    path 'run_settings.yml', emit: run_settings
    path 'input_checksums.tsv', emit: input_checksums
    path 'input_checks.tsv', emit: input_checks
    path 'reference_integrity.tsv', emit: reference_integrity
    path 'versions.yml', emit: versions

    script:
    method_arg = methods.join(',')
    """
    Rscript ${validator_script} \
        --run-name '${run_name}' \
        --target-manifest '${target_manifest}' \
        --gwas-manifest '${gwas_manifest}' \
        --reference-manifest '${reference_manifest}' \
        --phenotype-file '${phenotype_file}' \
        --phenotype-models '${phenotype_models}' \
        --methods '${method_arg}' \
        --genome-build '${genome_build}' \
        --seed '${seed}' \
        --report-enabled '${report_enabled}' \
        --target-imputation '${target_imputation}' \
        --launch-dir '${launch_dir}' \
        --reference-base '${reference_base}' \
        --run-plan '${run_plan}'

    cat > versions.yml <<-VERSIONS
    "${task.process}:${run_name}":
        R: \$(Rscript -e 'cat(as.character(getRversion()))')
    VERSIONS
    """

    stub:
    """
    awk -F '\t' -v OFS='\t' -v base='${launch_dir}' '
        function rooted(value) {
            return value == "" || substr(value, 1, 1) == "/" || value ~ /^[A-Za-z]:/ ? value : base "/" value
        }
        NR > 1 {
            \$4 = rooted(\$4)
            \$5 = rooted(\$5)
            \$6 = rooted(\$6)
            \$11 = rooted(\$11)
            \$12 = rooted(\$12)
        }
        { print }
    ' ${target_manifest} > targets.tsv
    awk -F '\t' -v OFS='\t' -v base='${launch_dir}' '
        function rooted(value) {
            return value == "" || substr(value, 1, 1) == "/" || value ~ /^[A-Za-z]:/ ? value : base "/" value
        }
        NR > 1 { \$3 = rooted(\$3) }
        { print }
    ' ${gwas_manifest} > gwas.tsv
    awk -F '\t' -v OFS='\t' -v base='${reference_base}' '
        function rooted(value) {
            return value == "" || substr(value, 1, 1) == "/" || value ~ /^[A-Za-z]:/ ? value : base "/" value
        }
        NR > 1 {
            \$5 = rooted(\$5)
            \$6 = rooted(\$6)
        }
        { print }
    ' ${reference_manifest} > references.tsv
    cp ${phenotype_models} models.tsv
    cp ${run_plan} run_plan.tsv
    printf "run_name: '${run_name}'\ngenome_build: '${genome_build}'\nmethods:\n" > run_settings.yml
    printf 'path\talgorithm\tchecksum\nmanifest\tSHA-256\tstub\n' > input_checksums.tsv
    printf 'check\tvalue\tstatus\ntarget_cohorts\t1\tPASS\ngwas_traits\t1\tPASS\nreference_resources\t1\tPASS\n' > input_checks.tsv
    printf 'reference_id\treference_type\tdeclared_sha256\tobserved_sha256\tstatus\nTEST\tdbsnp\tstub\tstub\tAUTHENTICATED\n' > reference_integrity.tsv
    printf '"${task.process}:${run_name}":\n  R: stub\n' > versions.yml
    """
}
