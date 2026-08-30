include { VALIDATE_MANIFESTS } from '../modules/local/validate_manifests/main'
include { HARMONISE_GWAS } from '../modules/local/harmonise_gwas/main'
include { PREPARE_TARGET } from '../modules/local/prepare_target/main'
include { PLINK_REFERENCE_FREQ } from '../modules/local/plink_reference_freq/main'
include { PLINK_CLUMP } from '../modules/local/plink_clump/main'
include { BUILD_CT_WEIGHTS } from '../modules/local/build_ct_weights/main'
include { PLINK_SCORE } from '../modules/local/plink_score/main'
include { SBAYESRC_TIDY } from '../modules/local/sbayesrc_tidy/main'
include { SBAYESRC_IMPUTE } from '../modules/local/sbayesrc_impute/main'
include { SBAYESRC_MODEL } from '../modules/local/sbayesrc_model/main'
include { SBAYESRC_SCORE } from '../modules/local/sbayesrc_score/main'
include { COMBINE_SCORES } from '../modules/local/combine_scores/main'
include { SUMMARISE_GENERATION_QC } from '../modules/local/summarise_generation_qc/main'
include { PHENOTYPE_ASSOCIATION } from '../modules/local/phenotype_association/main'
include { RENDER_REPORT } from '../modules/local/render_report/main'
include { COLLECT_VERSIONS } from '../modules/local/collect_versions/main'

workflow DNAPRS {
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
    script_files
    report_source

    main:
    VALIDATE_MANIFESTS(
        run_name,
        target_manifest,
        gwas_manifest,
        reference_manifest,
        phenotype_file,
        phenotype_models,
        methods,
        genome_build,
        seed,
        report_enabled,
        launch_dir,
        script_files.validate,
    )

    gwas_rows = VALIDATE_MANIFESTS.out.gwas_manifest
        .splitCsv(header: true, sep: '\t')
        .map { row -> tuple(row, file(row.path)) }
    target_rows = VALIDATE_MANIFESTS.out.target_manifest
        .splitCsv(header: true, sep: '\t')
    reference_rows = VALIDATE_MANIFESTS.out.reference_manifest
        .splitCsv(header: true, sep: '\t')

    HARMONISE_GWAS(gwas_rows, script_files.harmonise)
    PREPARE_TARGET(target_rows, script_files.prepare_target)

    score_files = channel.empty()
    generation_qc_files = HARMONISE_GWAS.out.harmonised.map { _meta, _cojo, _clump_input, harmonisation_qc -> harmonisation_qc }
    result_files = VALIDATE_MANIFESTS.out.target_manifest.map { result_file -> tuple('run/preflight', result_file) }
        .mix(VALIDATE_MANIFESTS.out.gwas_manifest.map { result_file -> tuple('run/preflight', result_file) })
        .mix(VALIDATE_MANIFESTS.out.reference_manifest.map { result_file -> tuple('run/preflight', result_file) })
        .mix(VALIDATE_MANIFESTS.out.phenotype_models.map { result_file -> tuple('run/preflight', result_file) })
        .mix(VALIDATE_MANIFESTS.out.resolved_params.map { result_file -> tuple('run/preflight', result_file) })
        .mix(VALIDATE_MANIFESTS.out.input_checksums.map { result_file -> tuple('run/preflight', result_file) })
        .mix(VALIDATE_MANIFESTS.out.preflight_qc.map { result_file -> tuple('qc/preflight', result_file) })
        .mix(HARMONISE_GWAS.out.harmonised.map { meta, cojo, _clump_input, _harmonisation_qc -> tuple("gwas/${meta.trait_id}", cojo) })
        .mix(HARMONISE_GWAS.out.harmonised.map { meta, _cojo, _clump_input, harmonisation_qc -> tuple("qc/gwas/${meta.trait_id}", harmonisation_qc) })
        .mix(PREPARE_TARGET.out.prepared.map { meta, _target_dir, target_qc -> tuple("qc/target/${meta.cohort}", target_qc) })
    version_files = VALIDATE_MANIFESTS.out.versions
        .mix(HARMONISE_GWAS.out.versions)
        .mix(PREPARE_TARGET.out.versions)

    if (methods.contains('plink_ct')) {
        plink_reference = reference_rows
            .filter { row -> row.reference_type == 'plink_ld' }
            .first()
        PLINK_REFERENCE_FREQ(plink_reference)

        clump_input = HARMONISE_GWAS.out.harmonised.combine(plink_reference)
        PLINK_CLUMP(clump_input)
        BUILD_CT_WEIGHTS(PLINK_CLUMP.out.clumped, script_files.ct_weights)

        plink_score_input = PREPARE_TARGET.out.prepared
            .combine(BUILD_CT_WEIGHTS.out.weights)
            .combine(PLINK_REFERENCE_FREQ.out.frequency)
        PLINK_SCORE(plink_score_input, script_files.parse_plink)
        score_files = score_files.mix(PLINK_SCORE.out.scores.map { _target, _gwas, score, _score_qc -> score })
        generation_qc_files = generation_qc_files.mix(
            BUILD_CT_WEIGHTS.out.weights.map { _meta, _weight, weight_qc, _harmonisation_qc, _clump_log -> weight_qc }
        )

        result_files = result_files
            .mix(PLINK_REFERENCE_FREQ.out.frequency.map { _reference, frequency, _reference_log -> tuple('reference/plink_ct', frequency) })
            .mix(PLINK_REFERENCE_FREQ.out.frequency.map { _reference, _frequency, reference_log -> tuple('logs/plink_ct/reference', reference_log) })
            .mix(BUILD_CT_WEIGHTS.out.weights.map { meta, weight, _weight_qc, _harmonisation_qc, _clump_log -> tuple("plink_ct/${meta.trait_id}", weight) })
            .mix(BUILD_CT_WEIGHTS.out.weights.map { meta, _weight, weight_qc, _harmonisation_qc, _clump_log -> tuple("qc/plink_ct/${meta.trait_id}", weight_qc) })
            .mix(BUILD_CT_WEIGHTS.out.weights.map { meta, _weight, _weight_qc, _harmonisation_qc, clump_log -> tuple("qc/plink_ct/${meta.trait_id}", clump_log) })
            .mix(PLINK_SCORE.out.scores.map { target, gwas, score, _score_qc -> tuple("scores/${target.cohort}/${gwas.trait_id}", score) })
            .mix(PLINK_SCORE.out.scores.map { target, gwas, _score, score_qc -> tuple("qc/scores/${target.cohort}/${gwas.trait_id}", score_qc) })
            .mix(PLINK_SCORE.out.logs.map { target, gwas, score_log -> tuple("logs/plink_ct/${target.cohort}/${gwas.trait_id}", score_log) })
        version_files = version_files
            .mix(PLINK_REFERENCE_FREQ.out.versions)
            .mix(PLINK_CLUMP.out.versions)
            .mix(BUILD_CT_WEIGHTS.out.versions)
            .mix(PLINK_SCORE.out.versions)
    }

    if (methods.contains('sbayesrc')) {
        sbayesrc_ld = reference_rows
            .filter { row -> row.reference_type == 'sbayesrc_ld' }
            .first()
        annotation = reference_rows
            .filter { row -> row.reference_type == 'annotation' }
            .first()

        tidy_input = HARMONISE_GWAS.out.harmonised.combine(sbayesrc_ld)
        SBAYESRC_TIDY(tidy_input, script_files.sbayesrc)
        impute_input = SBAYESRC_TIDY.out.tidy.combine(sbayesrc_ld)
        SBAYESRC_IMPUTE(impute_input, script_files.sbayesrc)
        model_input = SBAYESRC_IMPUTE.out.imputed
            .combine(sbayesrc_ld)
            .combine(annotation)
        SBAYESRC_MODEL(model_input, script_files.sbayesrc, seed)
        sbayesrc_score_input = PREPARE_TARGET.out.prepared.combine(SBAYESRC_MODEL.out.model)
        SBAYESRC_SCORE(sbayesrc_score_input, script_files.sbayesrc)
        score_files = score_files.mix(SBAYESRC_SCORE.out.scores.map { _target, _gwas, score, _score_qc -> score })
        generation_qc_files = generation_qc_files
            .mix(SBAYESRC_TIDY.out.tidy.map { _meta, _tidy, tidy_qc, _harmonisation_qc -> tidy_qc })
            .mix(SBAYESRC_IMPUTE.out.imputed.map { _meta, _imputed, impute_qc, _tidy_qc, _harmonisation_qc -> impute_qc })
            .mix(SBAYESRC_MODEL.out.model.map { _meta, _weight, _parameter, model_qc, _impute_qc, _tidy_qc, _harmonisation_qc -> model_qc })

        result_files = result_files
            .mix(SBAYESRC_TIDY.out.tidy.map { meta, tidy, _tidy_qc, _harmonisation_qc -> tuple("sbayesrc/${meta.trait_id}/summary", tidy) })
            .mix(SBAYESRC_TIDY.out.tidy.map { meta, _tidy, tidy_qc, _harmonisation_qc -> tuple("qc/sbayesrc/${meta.trait_id}", tidy_qc) })
            .mix(SBAYESRC_TIDY.out.logs.map { meta, stage_log -> tuple("logs/sbayesrc/${meta.trait_id}", stage_log) })
            .mix(SBAYESRC_IMPUTE.out.imputed.map { meta, imputed, _impute_qc, _tidy_qc, _harmonisation_qc -> tuple("sbayesrc/${meta.trait_id}/summary", imputed) })
            .mix(SBAYESRC_IMPUTE.out.imputed.map { meta, _imputed, impute_qc, _tidy_qc, _harmonisation_qc -> tuple("qc/sbayesrc/${meta.trait_id}", impute_qc) })
            .mix(SBAYESRC_IMPUTE.out.logs.map { meta, stage_log -> tuple("logs/sbayesrc/${meta.trait_id}", stage_log) })
            .mix(SBAYESRC_MODEL.out.model.map { meta, weight, _parameter, _model_qc, _impute_qc, _tidy_qc, _harmonisation_qc -> tuple("sbayesrc/${meta.trait_id}/model", weight) })
            .mix(SBAYESRC_MODEL.out.model.map { meta, _weight, parameter, _model_qc, _impute_qc, _tidy_qc, _harmonisation_qc -> tuple("sbayesrc/${meta.trait_id}/model", parameter) })
            .mix(SBAYESRC_MODEL.out.model.map { meta, _weight, _parameter, model_qc, _impute_qc, _tidy_qc, _harmonisation_qc -> tuple("qc/sbayesrc/${meta.trait_id}", model_qc) })
            .mix(SBAYESRC_MODEL.out.logs.map { meta, stage_log -> tuple("logs/sbayesrc/${meta.trait_id}", stage_log) })
            .mix(SBAYESRC_SCORE.out.scores.map { target, gwas, score, _score_qc -> tuple("scores/${target.cohort}/${gwas.trait_id}", score) })
            .mix(SBAYESRC_SCORE.out.scores.map { target, gwas, _score, score_qc -> tuple("qc/scores/${target.cohort}/${gwas.trait_id}", score_qc) })
            .mix(SBAYESRC_SCORE.out.logs.map { target, gwas, stage_log -> tuple("logs/sbayesrc/${gwas.trait_id}/${target.cohort}", stage_log) })
        version_files = version_files
            .mix(SBAYESRC_TIDY.out.versions)
            .mix(SBAYESRC_IMPUTE.out.versions)
            .mix(SBAYESRC_MODEL.out.versions)
            .mix(SBAYESRC_SCORE.out.versions)
    }

    COMBINE_SCORES(score_files.collect(), script_files.combine)
    SUMMARISE_GENERATION_QC(
        generation_qc_files.collect(),
        COMBINE_SCORES.out.score_qc,
        script_files.generation_qc,
    )
    PHENOTYPE_ASSOCIATION(
        COMBINE_SCORES.out.scores_long,
        phenotype_file,
        phenotype_models,
        script_files.association,
    )

    result_files = result_files
        .mix(COMBINE_SCORES.out.scores_long.map { result_file -> tuple('scores/combined', result_file) })
        .mix(COMBINE_SCORES.out.scores_wide.map { result_file -> tuple('scores/combined', result_file) })
        .mix(COMBINE_SCORES.out.score_qc.map { result_file -> tuple('qc/scores', result_file) })
        .mix(COMBINE_SCORES.out.concordance.map { result_file -> tuple('qc/scores', result_file) })
        .mix(SUMMARISE_GENERATION_QC.out.variant_flow.map { result_file -> tuple('qc/generation', result_file) })
        .mix(PHENOTYPE_ASSOCIATION.out.associations.map { result_file -> tuple('phenotype', result_file) })
        .mix(PHENOTYPE_ASSOCIATION.out.fitted_models.map { result_file -> tuple('phenotype', result_file) })
        .mix(PHENOTYPE_ASSOCIATION.out.plot_data.map { result_file -> tuple('phenotype', result_file) })
    version_files = version_files
        .mix(COMBINE_SCORES.out.versions)
        .mix(SUMMARISE_GENERATION_QC.out.versions)
        .mix(PHENOTYPE_ASSOCIATION.out.versions)

    COLLECT_VERSIONS(version_files.collect())
    result_files = result_files.mix(COLLECT_VERSIONS.out.versions.map { result_file -> tuple('pipeline_info', result_file) })

    if (report_enabled) {
        output_manifest = result_files
            .map { publish_path, result_file -> "${publish_path}\t${result_file.name}" }
            .collectFile(
                name: 'output_files.tsv',
                seed: 'publish_path\tfile_name\n',
                sort: true,
                newLine: true,
            )
        report_file_list = result_files
            .map { _publish_path, result_file -> result_file }
            .collect()
        RENDER_REPORT(
            report_file_list,
            output_manifest,
            report_source,
        )
        published_files = result_files
            .mix(RENDER_REPORT.out.pages.flatten().map { report_file -> tuple('', report_file) })
            .mix(RENDER_REPORT.out.libraries.map { report_file -> tuple('', report_file) })
            .mix(RENDER_REPORT.out.assets.map { report_file -> tuple('', report_file) })
            .mix(RENDER_REPORT.out.downloads.map { report_file -> tuple('', report_file) })
            .mix(RENDER_REPORT.out.figures.map { report_file -> tuple('', report_file) })
            .mix(RENDER_REPORT.out.provenance.flatten().map { report_file -> tuple('', report_file) })
    } else {
        published_files = result_files
    }

    emit:
    published_files.map { publish_path, result_file -> tuple(publish_path, result_file) }
}
