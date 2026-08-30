process VALIDATE_MANIFESTS {
    tag "$run_name"
    label 'process_low'
    label 'process_r'

    input:
    val run_name
    path target_manifest
    path gwas_manifest
    path reference_manifest
    path phenotype_file
    path phenotype_models
    val methods
    val genome_build
    val seed
    val report_enabled
    val launch_dir
    path validator_script

    output:
    path 'target_manifest.validated.tsv', emit: target_manifest
    path 'gwas_manifest.validated.tsv', emit: gwas_manifest
    path 'reference_manifest.validated.tsv', emit: reference_manifest
    path 'phenotype_models.validated.tsv', emit: phenotype_models
    path 'resolved_params.yaml', emit: resolved_params
    path 'input_checksums.tsv', emit: input_checksums
    path 'preflight_qc.tsv', emit: preflight_qc
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
        --launch-dir '${launch_dir}'

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
        }
        { print }
    ' ${target_manifest} > target_manifest.validated.tsv
    awk -F '\t' -v OFS='\t' -v base='${launch_dir}' '
        function rooted(value) {
            return value == "" || substr(value, 1, 1) == "/" || value ~ /^[A-Za-z]:/ ? value : base "/" value
        }
        NR > 1 { \$3 = rooted(\$3) }
        { print }
    ' ${gwas_manifest} > gwas_manifest.validated.tsv
    awk -F '\t' -v OFS='\t' -v base='${launch_dir}' '
        function rooted(value) {
            return value == "" || substr(value, 1, 1) == "/" || value ~ /^[A-Za-z]:/ ? value : base "/" value
        }
        NR > 1 { \$3 = rooted(\$3) }
        { print }
    ' ${reference_manifest} > reference_manifest.validated.tsv
    cp ${phenotype_models} phenotype_models.validated.tsv
    printf "run_name: '${run_name}'\ngenome_build: '${genome_build}'\nmethods:\n" > resolved_params.yaml
    printf 'path\talgorithm\tchecksum\nmanifest\tSHA-256\tstub\n' > input_checksums.tsv
    printf 'check\tvalue\tstatus\ntarget_cohorts\t1\tPASS\ngwas_traits\t1\tPASS\nreference_resources\t1\tPASS\n' > preflight_qc.tsv
    printf '"${task.process}:${run_name}":\n  R: stub\n' > versions.yml
    """
}
