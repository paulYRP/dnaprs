#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    paulYRP/dnaprs
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/paulYRP/dnaprs
----------------------------------------------------------------------------------------
*/

include { DNAPRS                  } from './workflows/dnaprs'
include { RESOLVE_INPUTS          } from './modules/local/resolve_inputs/main'
include { REFERENCE_PLAN          } from './modules/local/reference_plan/main'
include { REFERENCE_ASSET         } from './modules/local/reference_asset/main'
include { ASSEMBLE_REFERENCES     } from './modules/local/assemble_references/main'
include { COLLECT_VERSIONS as COLLECT_REFERENCE_VERSIONS } from './modules/local/collect_versions/main'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_dnaprs_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_dnaprs_pipeline'

workflow NFCORE_DNAPRS {
    take:
    run_name
    target_manifest
    gwas_manifest
    reference_manifest
    phenotype_file
    phenotype_models
    phenotype_enabled
    methods
    genome_build
    seed
    report_enabled
    imputation_variant_missingness
    direct_variant_missingness
    sample_missingness
    target_maf
    target_hwe
    ancestry_pcs
    ancestry_percentile
    target_imputation
    imputation_dr2
    stop_after
    launch_dir
    reference_base
    input_checks
    input_versions
    run_plan

    main:
    resolved_phenotype_file = phenotype_file ?: file("${projectDir}/assets/empty_phenotype.tsv")
    resolved_phenotype_models = phenotype_models ?: file("${projectDir}/assets/empty_models.tsv")
    if (stop_after == 'imputation' && !target_imputation) {
        error "stop_after=imputation requires target_imputation=true"
    }
    script_files = [
        validate: file("${projectDir}/bin/validate_manifests.R"),
        genotype_eda: file("${projectDir}/bin/genotype_eda.sh"),
        harmonise: file("${projectDir}/bin/harmonise_gwas.R"),
        align_plink_gwas: file("${projectDir}/bin/align_plink_gwas.R"),
        prepare_target: file("${projectDir}/bin/prepare_target.sh"),
        target_adapter: file("${projectDir}/bin/target_adapter.pl"),
        marker_resolver: file("${projectDir}/bin/resolve_target_markers.R"),
        target_qc: file("${projectDir}/bin/target_qc.sh"),
        participant_decisions: file("${projectDir}/bin/participant_decisions.R"),
        reference_ancestry: file("${projectDir}/bin/reference_ancestry.sh"),
        classify_ancestry: file("${projectDir}/bin/classify_ancestry.R"),
        target_impute_chromosome: file("${projectDir}/bin/target_impute_chromosome.sh"),
        assemble_target_imputation: file("${projectDir}/bin/assemble_target_imputation.sh"),
        prepare_plink_reference: file("${projectDir}/bin/prepare_plink_reference.sh"),
        prepare_sbayesrc_reference: file("${projectDir}/bin/prepare_sbayesrc_reference.sh"),
        ct_weights: file("${projectDir}/bin/build_ct_weights.R"),
        parse_plink: file("${projectDir}/bin/parse_plink_score.R"),
        audit_scoring: file("${projectDir}/bin/audit_scoring_variants.R"),
        compare_direct_score: file("${projectDir}/bin/compare_direct_score.R"),
        sbayesrc: file("${projectDir}/bin/run_sbayesrc.R"),
        combine: file("${projectDir}/bin/combine_scores.R"),
        combine_phenotype: file("${projectDir}/bin/combine_phenotype.R"),
        collect_versions: file("${projectDir}/bin/collect_versions.R"),
        association: file("${projectDir}/bin/phenotype_association.R"),
        generation_qc: file("${projectDir}/bin/summarise_generation_qc.R"),
    ]
    report_source = file("${projectDir}/assets/report")

    DNAPRS(
        run_name,
        target_manifest,
        gwas_manifest,
        reference_manifest,
        resolved_phenotype_file,
        resolved_phenotype_models,
        phenotype_enabled,
        methods,
        genome_build,
        seed,
        report_enabled,
        imputation_variant_missingness,
        direct_variant_missingness,
        sample_missingness,
        target_maf,
        target_hwe,
        ancestry_pcs,
        ancestry_percentile,
        target_imputation,
        imputation_dr2,
        stop_after,
        launch_dir,
        reference_base,
        script_files,
        report_source,
        input_checks,
        input_versions,
        run_plan,
    )

    emit:
    DNAPRS.out
}

workflow {
    main:
    run_outdir = file("${params.outdir}/${params.run_name}")
    run_outdir_file = run_outdir.toFile()
    if (run_outdir_file.exists() && !run_outdir_file.isDirectory()) {
        error "Run output path exists but is not a directory: ${run_outdir}"
    }
    // Nextflow creates these configured provenance files before workflow code starts.
    // They belong to the current launch and must not make a new output directory look
    // stale. Every other existing file remains protected by the collision guard.
    bootstrap_provenance = [
        'reports/provenance/execution_report.html',
        'reports/provenance/execution_timeline.html',
        'reports/provenance/execution_trace.txt',
        'reports/provenance/pipeline_dag.html',
    ] as Set
    launch_started_ms = workflow.start.toInstant().toEpochMilli()
    existing_run_files = []
    if (run_outdir_file.isDirectory()) {
        run_outdir_file.eachFileRecurse(groovy.io.FileType.FILES) { existing_file ->
            def relative_path = run_outdir_file.toPath()
                .relativize(existing_file.toPath())
                .toString()
                .replace('\\', '/')
            def current_bootstrap = bootstrap_provenance.contains(relative_path) &&
                existing_file.lastModified() >= launch_started_ms
            if (!current_bootstrap) existing_run_files << relative_path
        }
    }
    if (existing_run_files && !workflow.resume && !params.overwrite) {
        error "Run output directory contains existing results: ${run_outdir}. Use -resume for the same run or set --overwrite true after reviewing its contents."
    }
    selected_methods = params.methods.tokenize(',').collect { it.trim() }.findAll { it }.unique()
    if (!selected_methods || selected_methods.any { !['plink_ct', 'sbayesrc'].contains(it) }) {
        error "--methods must contain one or more of: plink_ct,sbayesrc"
    }
    empty_input = file("${projectDir}/assets/empty_input")
    target_source = params.reference_only ? empty_input : (params.input instanceof CharSequence ? file(params.input, checkIfExists: true) : empty_input)
    run_prs_requested = !params.reference_only && ['prs', 'phenotype', 'report'].contains(params.stop_after)
    gwas_source = run_prs_requested ? (params.gwas instanceof CharSequence ? file(params.gwas, checkIfExists: true) : empty_input) : empty_input
    reference_source = params.references instanceof CharSequence && params.references ? file(params.references, checkIfExists: true) : empty_input
    phenotype_file_input = params.phenotype ? file(params.phenotype, checkIfExists: true) : null

    input_spec = groovy.json.JsonOutput.toJson(params.input).getBytes('UTF-8').encodeBase64().toString()
    gwas_spec = groovy.json.JsonOutput.toJson(params.gwas).getBytes('UTF-8').encodeBase64().toString()
    reference_spec = groovy.json.JsonOutput.toJson(params.references ?: '').getBytes('UTF-8').encodeBase64().toString()
    selected_models = params.models ?: []
    if (params.model_id) {
        selected_models = selected_models.findAll { model -> model.id == params.model_id }
        if (selected_models.size() != 1) {
            error "--model_id '${params.model_id}' must match exactly one item in the YAML models list."
        }
    }
    models_spec = groovy.json.JsonOutput.toJson(selected_models).getBytes('UTF-8').encodeBase64().toString()

    PIPELINE_INITIALISATION(
        params.version,
        params.validate_params,
        params.monochrome_logs,
        [],
        run_outdir,
        params.help,
        params.help_full,
        params.show_hidden,
    )

    RESOLVE_INPUTS(
        input_spec,
        gwas_spec,
        reference_spec,
        models_spec,
        params.genome,
        selected_methods,
        params.target_imputation,
        params.stop_after,
        params.reference_only,
        params.reference_mode,
        params.reference_bundle,
        params.phenotype ?: '',
        params.outcome ?: '',
        params.covariates ?: '',
        params.model_type ?: '',
        params.participant_id ?: '',
        params.group_column ?: '',
        params.control_value ?: '',
        params.case_value ?: '',
        params.beagle_jar ?: '',
        params.unbref3_jar ?: '',
        target_source,
        gwas_source,
        reference_source,
        file("${projectDir}/bin/resolve_inputs.R"),
    )

    resolved_references = RESOLVE_INPUTS.out.references
    resolved_input_checks = RESOLVE_INPUTS.out.checks
    resolved_input_versions = RESOLVE_INPUTS.out.versions
    if (params.reference_mode != 'local') {
        REFERENCE_PLAN(
            RESOLVE_INPUTS.out.references,
            file("${projectDir}/assets/reference_catalogue.tsv"),
            RESOLVE_INPUTS.out.run_plan,
            params.reference_mode,
            file("${projectDir}/bin/plan_references.R"),
        )
        reference_assets = REFERENCE_PLAN.out.assets
            .splitCsv(header: true, sep: '\t')
            .map { row ->
                def cached = row.asset_id == 'cache_complete' ? empty_input :
                    file("${params.reference_dir}/${params.reference_bundle}/${row.relative_path}", checkIfExists: false)
                def cached_input = row.asset_id != 'cache_complete' && java.nio.file.Files.exists(cached) ? cached : empty_input
                tuple(row, cached_input, cached.toString())
            }
        REFERENCE_ASSET(
            reference_assets,
            file("${projectDir}/assets/reference/human_g1k_v37.fasta.fai"),
            file("${projectDir}/bin/reference_asset.sh"),
        )
        ASSEMBLE_REFERENCES(
            REFERENCE_ASSET.out.asset.map { row, _asset -> row }.collect(),
            REFERENCE_ASSET.out.asset.map { _row, asset -> asset }.collect(),
            RESOLVE_INPUTS.out.references,
            params.reference_bundle,
            params.genome,
            file("${params.reference_dir}/${params.reference_bundle}", checkIfExists: false).toString(),
            file("${projectDir}/bin/assemble_references.R"),
        )
        resolved_references = ASSEMBLE_REFERENCES.out.references
        resolved_input_checks = resolved_input_checks
            .mix(REFERENCE_PLAN.out.plan)
            .mix(ASSEMBLE_REFERENCES.out.receipt)
        resolved_input_versions = resolved_input_versions
            .mix(REFERENCE_PLAN.out.versions)
            .mix(REFERENCE_ASSET.out.versions)
            .mix(ASSEMBLE_REFERENCES.out.versions)
    }

    if (params.reference_only) {
        COLLECT_REFERENCE_VERSIONS(
            resolved_input_versions.collect(),
            file("${projectDir}/bin/collect_versions.R"),
        )
        pipeline_results = resolved_references.map { result_file -> tuple('data/inputs', result_file) }
            .mix(RESOLVE_INPUTS.out.settings.map { result_file -> tuple('data/inputs', result_file) })
            .mix(RESOLVE_INPUTS.out.run_plan.map { result_file -> tuple('data/inputs', result_file) })
            .mix(resolved_input_checks.map { result_file -> tuple('logs/references', result_file) })
            .mix(COLLECT_REFERENCE_VERSIONS.out.versions.map { result_file -> tuple('data/pipeline_info', result_file) })
    } else {
        NFCORE_DNAPRS(
            params.run_name,
            RESOLVE_INPUTS.out.targets,
            RESOLVE_INPUTS.out.gwas,
            resolved_references,
            phenotype_file_input,
            RESOLVE_INPUTS.out.models,
            params.phenotype ? true : false,
            selected_methods,
            params.genome,
            params.seed,
            params.report_enabled,
            params.imputation_variant_missingness,
            params.direct_variant_missingness,
            params.sample_missingness,
            params.maf_filter,
            params.hwe_filter,
            params.ancestry_pcs,
            params.ancestry_percentile,
            params.target_imputation,
            params.imputation_dr2,
            params.stop_after,
            launchDir,
            launchDir,
            resolved_input_checks,
            resolved_input_versions,
            RESOLVE_INPUTS.out.run_plan,
        )
        pipeline_results = NFCORE_DNAPRS.out
    }

    PIPELINE_COMPLETION(
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        run_outdir,
        params.monochrome_logs,
    )

    publish:
    results = pipeline_results
}

output {
    results {
        path { publish_path, _result_file -> publish_path ? "${params.run_name}/${publish_path}" : "${params.run_name}" }
        index {
            path "${params.run_name}/data/inputs/published_files.json"
        }
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
