#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-core/dnaprs
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/nf-core/dnaprs
    Website: https://nf-co.re/dnaprs
----------------------------------------------------------------------------------------
*/

include { DNAPRS                  } from './workflows/dnaprs'
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
    methods
    genome_build
    seed
    report_enabled
    launch_dir

    main:
    resolved_phenotype_file = phenotype_file ?: file("${projectDir}/assets/empty_phenotype.tsv")
    resolved_phenotype_models = phenotype_models ?: file("${projectDir}/assets/empty_models.tsv")
    script_files = [
        validate: file("${projectDir}/bin/validate_manifests.R"),
        harmonise: file("${projectDir}/bin/harmonise_gwas.R"),
        prepare_target: file("${projectDir}/bin/prepare_target.sh"),
        ct_weights: file("${projectDir}/bin/build_ct_weights.R"),
        parse_plink: file("${projectDir}/bin/parse_plink_score.R"),
        sbayesrc: file("${projectDir}/bin/run_sbayesrc.R"),
        combine: file("${projectDir}/bin/combine_scores.R"),
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
        methods,
        genome_build,
        seed,
        report_enabled,
        launch_dir,
        script_files,
        report_source,
    )

    emit:
    results = DNAPRS.out
}

workflow {
    main:
    run_outdir = file("${params.outdir}/${params.run_name}")

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

    NFCORE_DNAPRS(
        params.run_name,
        params.input,
        params.gwas_manifest,
        params.reference_manifest,
        params.phenotype_file,
        params.phenotype_models,
        params.methods.tokenize(','),
        params.genome_build,
        params.seed,
        params.report_enabled,
        params.manifest_base ?: launchDir,
    )

    PIPELINE_COMPLETION(
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        run_outdir,
        params.monochrome_logs,
    )

    publish:
    results = NFCORE_DNAPRS.out.results
}

output {
    results {
        path { publish_path, _result_file -> publish_path ? "${params.run_name}/${publish_path}" : "${params.run_name}" }
        index {
            path "${params.run_name}/run/published_files.json"
        }
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
