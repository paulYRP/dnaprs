include { VALIDATE_MANIFESTS } from '../modules/local/validate_manifests/main'
include { HARMONISE_GWAS } from '../modules/local/harmonise_gwas/main'
include { PREPARE_TARGET } from '../modules/local/prepare_target/main'
include { TARGET_QC } from '../modules/local/target_qc/main'
include { PARTICIPANT_DECISIONS } from '../modules/local/participant_decisions/main'
include { REFERENCE_ANCESTRY } from '../modules/local/reference_ancestry/main'
include { TARGET_IMPUTE_CHROMOSOME } from '../modules/local/target_impute/main'
include { ASSEMBLE_TARGET_IMPUTATION } from '../modules/local/assemble_target_imputation/main'
include { PREPARE_PLINK_REFERENCE } from '../modules/local/prepare_plink_reference/main'
include { PREPARE_SBAYESRC_REFERENCE } from '../modules/local/prepare_sbayesrc_reference/main'
include { GENOTYPE_EDA } from '../modules/local/genotype_eda/main'
include { GENOTYPE_EDA as TARGET_QC_REVIEW } from '../modules/local/genotype_eda/main'
include { PLINK_REFERENCE_FREQ } from '../modules/local/plink_reference_freq/main'
include { ALIGN_PLINK_GWAS } from '../modules/local/align_plink_gwas/main'
include { PLINK_CLUMP } from '../modules/local/plink_clump/main'
include { BUILD_CT_WEIGHTS } from '../modules/local/build_ct_weights/main'
include { PLINK_SCORE } from '../modules/local/plink_score/main'
include { PARSE_PLINK_SCORE } from '../modules/local/parse_plink_score/main'
include { PLINK_SCORE as PLINK_DIRECT_SCORE } from '../modules/local/plink_score/main'
include { PARSE_PLINK_SCORE as PARSE_PLINK_DIRECT_SCORE } from '../modules/local/parse_plink_score/main'
include { COMPARE_DIRECT_SCORE } from '../modules/local/compare_direct_score/main'
include { SBAYESRC_TIDY } from '../modules/local/sbayesrc_tidy/main'
include { SBAYESRC_IMPUTE } from '../modules/local/sbayesrc_impute/main'
include { SBAYESRC_MODEL } from '../modules/local/sbayesrc_model/main'
include { SBAYESRC_SCORE } from '../modules/local/sbayesrc_score/main'
include { COMBINE_SCORES } from '../modules/local/combine_scores/main'
include { SUMMARISE_GENERATION_QC } from '../modules/local/summarise_generation_qc/main'
include { PHENOTYPE_ASSOCIATION } from '../modules/local/phenotype_association/main'
include { COMBINE_PHENOTYPE } from '../modules/local/combine_phenotype/main'
include { RENDER_REPORT } from '../modules/local/render_report/main'
include { REPORT_SOFTWARE } from '../modules/local/report_software/main'
include { PUBLIC_FIGURES } from '../modules/local/public_figures/main'
include { COLLECT_VERSIONS } from '../modules/local/collect_versions/main'

// Resolve every target companion as a Nextflow path so containers receive only declared
// inputs and task hashes include the actual genotype files. The stable source path remains
// in meta.source_genotype for provenance; meta.genotype names the staged task-local input.
def targetBaseName(value) {
    value.toString().replace('\\', '/').tokenize('/').last()
}

// User paths are relative to the directory from which Nextflow was launched. Resolve
// them explicitly because this workflow can be loaded from a remote project cache.
def userInputPath(value, launchDir) {
    def declared = value.toString()
    if (declared ==~ /^[A-Za-z][A-Za-z0-9+.-]*:\/\/.*/) return file(declared)

    def candidate = java.nio.file.Paths.get(declared)
    candidate.isAbsolute() ? candidate.normalize() : launchDir.resolve(candidate).normalize()
}

def requiredUserInput(value, launchDir, description) {
    def resolved = userInputPath(value, launchDir)
    try {
        file(resolved, checkIfExists: true)
    } catch (Exception _cause) {
        def message = "${description} cannot be accessed: '${value}' (resolved to '${resolved}')"
        log.error(message)
        error message
    }
}

def addTargetInput(sourceFiles, value, launchDir, description) {
    sourceFiles << requiredUserInput(value, launchDir, description)
}

def addTargetGenotype(sourceFiles, value, format, launchDir, cohort) {
    if (format == 'pgen') {
        def prefix = value.replaceFirst(/(?i)\.pgen$/, '')
        addTargetInput(sourceFiles, "${prefix}.pgen", launchDir, "Target '${cohort}' PGEN file")
        addTargetInput(sourceFiles, "${prefix}.pvar", launchDir, "Target '${cohort}' PVAR companion")
        addTargetInput(sourceFiles, "${prefix}.psam", launchDir, "Target '${cohort}' PSAM companion")
    } else if (format == 'bed') {
        def prefix = value.replaceFirst(/(?i)\.bed$/, '')
        addTargetInput(sourceFiles, "${prefix}.bed", launchDir, "Target '${cohort}' BED file")
        addTargetInput(sourceFiles, "${prefix}.bim", launchDir, "Target '${cohort}' BIM companion")
        addTargetInput(sourceFiles, "${prefix}.fam", launchDir, "Target '${cohort}' FAM companion")
    } else if (format == 'ped') {
        def prefix = value.replaceFirst(/(?i)\.ped$/, '')
        addTargetInput(sourceFiles, "${prefix}.ped", launchDir, "Target '${cohort}' PED file")
        addTargetInput(sourceFiles, "${prefix}.map", launchDir, "Target '${cohort}' MAP companion")
    } else {
        addTargetInput(sourceFiles, value, launchDir, "Target '${cohort}' genotype")
    }
}

def targetInputFiles(row, launchDir) {
    def sourceGenotype = row.genotype.toString()
    def sourceSample = row.sample?.toString() ?: ''
    def sourceKeep = row.keep?.toString() ?: ''
    def sourceAssayManifest = row.assay_manifest?.toString() ?: ''
    def sourceMarkerMap = row.marker_map?.toString() ?: ''
    def format = (row.source_format ?: row.format).toString()
    def sourceFiles = []
    def hasChromosomePattern = sourceGenotype.contains('{chr}') ||
        sourceGenotype.contains('{CHR}') ||
        sourceGenotype.contains('{chromosome}') ||
        sourceGenotype.contains('#')

    if (sourceGenotype.contains('*') || sourceGenotype.contains('?')) {
        error "Target '${row.cohort}' uses a raw glob. Use {chr}, {CHR}, {chromosome}, or # as the chromosome placeholder."
    }

    if (hasChromosomePattern) {
        (1..22).each { chromosome ->
            def value = sourceGenotype
                .replace('{chr}', chromosome.toString())
                .replace('{CHR}', chromosome.toString())
                .replace('{chromosome}', chromosome.toString())
                .replace('#', chromosome.toString())
            def primary = format == 'pgen' ? "${value.replaceFirst(/(?i)\.pgen$/, '')}.pgen" :
                format == 'bed' ? "${value.replaceFirst(/(?i)\.bed$/, '')}.bed" : value
            if (java.nio.file.Files.exists(userInputPath(primary, launchDir))) {
                addTargetGenotype(sourceFiles, value, format, launchDir, row.cohort)
            }
        }
        if (!sourceFiles) error "Target '${row.cohort}' chromosome placeholder matched no genotype files."
    } else {
        addTargetGenotype(sourceFiles, sourceGenotype, format, launchDir, row.cohort)
    }
    if (sourceSample) addTargetInput(sourceFiles, sourceSample, launchDir, "Target '${row.cohort}' sample file")
    if (sourceKeep) addTargetInput(sourceFiles, sourceKeep, launchDir, "Target '${row.cohort}' keep file")
    if (sourceAssayManifest) addTargetInput(sourceFiles, sourceAssayManifest, launchDir, "Target '${row.cohort}' assay manifest")
    if (sourceMarkerMap) addTargetInput(sourceFiles, sourceMarkerMap, launchDir, "Target '${row.cohort}' marker map")

    def stagedMeta = row + [
        source_genotype: sourceGenotype,
        source_sample: sourceSample,
        source_keep: sourceKeep,
        source_assay_manifest: sourceAssayManifest,
        source_marker_map: sourceMarkerMap,
        genotype: targetBaseName(sourceGenotype),
        sample: sourceSample ? targetBaseName(sourceSample) : '',
        keep: sourceKeep ? targetBaseName(sourceKeep) : '',
        assay_manifest: sourceAssayManifest ? targetBaseName(sourceAssayManifest) : '',
        marker_map: sourceMarkerMap ? targetBaseName(sourceMarkerMap) : '',
        source_format: format,
        format: format,
    ]
    tuple(stagedMeta, sourceFiles.unique())
}

def resolveReferenceInputPath(row, value, launchDir, field) {
    def matchedFiles = []
    def hasChromosomePattern = value.contains('{chr}') || value.contains('{CHR}') ||
        value.contains('{chromosome}') || value.contains('#')
    if (value.contains('*') || value.contains('?')) {
        error "Reference '${row.reference_id}' uses a raw glob. Use {chr}, {CHR}, {chromosome}, or #."
    }
    if (hasChromosomePattern) {
        (1..22).each { chromosome ->
            def resolved = value
                .replace('{chr}', chromosome.toString())
                .replace('{CHR}', chromosome.toString())
                .replace('{chromosome}', chromosome.toString())
                .replace('#', chromosome.toString())
            def resolvedPath = userInputPath(resolved, launchDir)
            if (java.nio.file.Files.exists(resolvedPath)) {
                matchedFiles << file(resolvedPath, checkIfExists: true)
            }
        }
        if (!matchedFiles) error "Reference '${row.reference_id}' chromosome placeholder matched no files."
    } else {
        matchedFiles << requiredUserInput(value, launchDir, "Reference '${row.reference_id}' ${field}")
    }
    matchedFiles
}

// Stage reference files or directories as declared Nextflow inputs. This keeps
// external bundles portable across Docker and Apptainer and adds their content to
// task hashes without copying or modifying the source bundle.
def referenceInputFiles(row, launchDir) {
    def sourcePath = row.path.toString()
    def sourceCompanion = row.companion?.toString() ?: ''
    def sourceFiles = resolveReferenceInputPath(row, sourcePath, launchDir, 'path')
    if (sourceCompanion) sourceFiles.addAll(resolveReferenceInputPath(row, sourceCompanion, launchDir, 'companion'))
    def stagedMeta = row + [
        source_path: sourcePath,
        source_companion: sourceCompanion,
        path: targetBaseName(sourcePath),
        companion: sourceCompanion ? targetBaseName(sourceCompanion) : '',
    ]
    tuple(stagedMeta, sourceFiles.unique())
}

workflow DNAPRS {
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
    script_files
    report_source
    input_checks
    run_plan

    main:
    run_prs = ['prs', 'phenotype', 'report'].contains(stop_after)
    run_phenotype = phenotype_enabled && ['phenotype', 'report'].contains(stop_after)

    // Declare every path inspected by the validator as a real Nextflow input.
    // This lets Docker and Apptainer mount data outside projectDir while the
    // generated records retain stable source paths for provenance.
    validation_target_assets = target_manifest
        .splitCsv(header: true, sep: '\t')
        .flatMap { row -> targetInputFiles(row, launch_dir)[1] }
        .ifEmpty { file("${projectDir}/assets/empty_input") }
        .collect()
    validation_gwas_assets = gwas_manifest
        .splitCsv(header: true, sep: '\t')
        .map { row -> requiredUserInput(row.path, launch_dir, "GWAS '${row.trait_id}' path") }
        .ifEmpty { file("${projectDir}/assets/empty_input") }
        .collect()
    validation_reference_assets = reference_manifest
        .splitCsv(header: true, sep: '\t')
        .flatMap { row -> referenceInputFiles(row, launch_dir)[1] }
        .ifEmpty { file("${projectDir}/assets/empty_input") }
        .collect()

    VALIDATE_MANIFESTS(
        run_name,
        target_manifest,
        gwas_manifest,
        reference_manifest,
        phenotype_file,
        phenotype_models,
        validation_target_assets,
        validation_gwas_assets,
        validation_reference_assets,
        methods,
        genome_build,
        seed,
        report_enabled,
        target_imputation,
        launch_dir,
        reference_base,
        run_plan,
        script_files.validate,
    )

    gwas_rows = VALIDATE_MANIFESTS.out.gwas_manifest
        .splitCsv(header: true, sep: '\t')
        .map { row -> tuple(row, requiredUserInput(row.path, launch_dir, "GWAS '${row.trait_id}' path")) }
    target_rows = VALIDATE_MANIFESTS.out.target_manifest
        .splitCsv(header: true, sep: '\t')
    target_inputs = target_rows.map { row -> targetInputFiles(row, launch_dir) }
    reference_rows = VALIDATE_MANIFESTS.out.reference_manifest
        .splitCsv(header: true, sep: '\t')
    reference_inputs = reference_rows.map { row -> referenceInputFiles(row, launch_dir) }

    dbsnp_source = reference_inputs
        .filter { row, _files -> row.reference_type == 'dbsnp' }
        .ifEmpty { tuple([reference_id: 'NOT_REQUIRED', path: ''], file("${projectDir}/assets/empty_input")) }
        .first()
    reference_fasta_source = reference_inputs
        .filter { row, _files -> row.reference_type == 'reference_fasta' }
        .ifEmpty { tuple([reference_id: 'NOT_REQUIRED', path: ''], file("${projectDir}/assets/empty_input")) }
        .first()

    if (run_prs) {
        HARMONISE_GWAS(gwas_rows, script_files.harmonise)
    }
    // Describe the untouched source genotypes before marker renaming, allele correction,
    // duplicate handling, or technical filtering. The EDA process performs only a
    // task-local format import and never edits the source files.
    GENOTYPE_EDA(target_inputs, script_files.genotype_eda, script_files.target_adapter)
    PREPARE_TARGET(
        target_inputs,
        dbsnp_source,
        reference_fasta_source,
        script_files.prepare_target,
        script_files.target_adapter,
        script_files.marker_resolver,
    )
    TARGET_QC(
        PREPARE_TARGET.out.prepared,
        script_files.target_qc,
        imputation_variant_missingness,
        direct_variant_missingness,
        sample_missingness,
        target_maf,
        target_hwe,
    )
    target_qc_review_inputs = TARGET_QC.out.imputation_ready.map { meta, target_dir, _qc ->
        def review_meta = meta + [
            source_cohort: meta.cohort,
            cohort: "${meta.cohort}_qc_review",
            format: 'pgen',
            genotype: "${target_dir.name}/${meta.cohort}",
            sample: '',
            dosage: 'DS',
            input_stage: 'qc_completed',
        ]
        tuple(review_meta, target_dir)
    }
    TARGET_QC_REVIEW(target_qc_review_inputs, script_files.genotype_eda, script_files.target_adapter)
    plink_panel_source = reference_inputs
        .filter { row, _files -> row.reference_type == 'imputation_panel' }
        .first()
    population_panel_source = reference_inputs
        .filter { row, _files -> row.reference_type == 'population_panel' }
        .first()
    related_samples_source = reference_inputs
        .filter { row, _files -> row.reference_type == 'related_samples' }
        .first()
    unbref3_jar_source = reference_inputs
        .filter { row, _files -> row.reference_type == 'unbref3_jar' }
        .first()
    PREPARE_PLINK_REFERENCE(
        plink_panel_source,
        population_panel_source,
        related_samples_source,
        unbref3_jar_source,
        script_files.prepare_plink_reference,
        genome_build,
    )
    plink_reference = PREPARE_PLINK_REFERENCE.out.reference.map { source, reference_dir ->
        def prepared = source + [
            reference_id: "${source.reference_id}_PLINK_LD",
            reference_type: 'plink_ld',
            reference_stage: 'prepared',
            source_format: 'pgen',
            path: "${reference_dir.name}/eur_reference",
        ]
        tuple(prepared, reference_dir)
    }
    ancestry_input = TARGET_QC.out.imputation_ready.combine(plink_reference)
    REFERENCE_ANCESTRY(
        ancestry_input,
        script_files.reference_ancestry,
        script_files.classify_ancestry,
        ancestry_pcs,
        ancestry_percentile,
    )
    participant_decision_inputs = TARGET_QC.out.sample_decisions
        .map { meta, decisions -> tuple(meta.cohort, meta, decisions) }
        .join(TARGET_QC_REVIEW.out.tables.map { meta, tables -> tuple(meta.source_cohort, tables) }, failOnDuplicate: true, failOnMismatch: true)
        .join(GENOTYPE_EDA.out.tables.map { meta, tables -> tuple(meta.cohort, tables) }, failOnDuplicate: true, failOnMismatch: true)
        .join(REFERENCE_ANCESTRY.out.target.map { meta, ancestry -> tuple(meta.cohort, ancestry) }, failOnDuplicate: true, failOnMismatch: true)
        .map { _cohort, meta, decisions, review_tables, raw_tables, ancestry ->
            tuple(meta, decisions, review_tables + raw_tables, ancestry)
        }
    PARTICIPANT_DECISIONS(participant_decision_inputs, script_files.participant_decisions)

    checkpoint_files = PREPARE_TARGET.out.checkpoint
        .filter { meta, _target_dir, _manifest -> ['raw', 'corrected'].contains(meta.input_stage) }
        .map { meta, target_dir, _manifest ->
            def stage = meta.input_stage == 'raw' ? 'corrected' : meta.input_stage
            tuple("checkpoints/${stage}", target_dir)
        }
        .mix(TARGET_QC.out.imputation_checkpoint.map { meta, target_dir, _manifest ->
            tuple("target/prepared/${meta.cohort}/imputation_ready", target_dir)
        })
        .mix(TARGET_QC.out.direct_checkpoint.map { meta, target_dir, _manifest ->
            tuple("target/prepared/${meta.cohort}/direct_ready", target_dir)
        })
        .mix(PREPARE_PLINK_REFERENCE.out.reference.map { _meta, reference_dir ->
            tuple('reference/plink_ct/prepared', reference_dir)
        })

    score_files = channel.empty()
    score_job_records = channel.empty()
    generation_qc_files = channel.empty()
    result_files = VALIDATE_MANIFESTS.out.target_manifest.map { result_file -> tuple('inputs', result_file) }
        .mix(VALIDATE_MANIFESTS.out.gwas_manifest.map { result_file -> tuple('inputs', result_file) })
        .mix(VALIDATE_MANIFESTS.out.reference_manifest.map { result_file -> tuple('inputs', result_file) })
        .mix(VALIDATE_MANIFESTS.out.phenotype_models.map { result_file -> tuple('inputs', result_file) })
        .mix(VALIDATE_MANIFESTS.out.run_plan.map { result_file -> tuple('inputs', result_file) })
        .mix(VALIDATE_MANIFESTS.out.run_settings.map { result_file -> tuple('inputs', result_file) })
        .mix(VALIDATE_MANIFESTS.out.input_checksums.map { result_file -> tuple('inputs', result_file) })
        .mix(VALIDATE_MANIFESTS.out.input_checks.map { result_file -> tuple('inputs/checks', result_file) })
        .mix(VALIDATE_MANIFESTS.out.reference_integrity.map { result_file -> tuple('inputs/checks', result_file) })
        .mix(input_checks.map { result_file -> tuple('inputs', result_file) })
        .mix(GENOTYPE_EDA.out.tables.flatMap { meta, tables -> tables.collect { result_file -> tuple("genotype_eda/${meta.cohort}", result_file) } })
        .mix(GENOTYPE_EDA.out.logs.map { meta, stage_log -> tuple("logs/genotype_eda/${meta.cohort}", stage_log) })
        .mix(TARGET_QC_REVIEW.out.tables.flatMap { meta, tables -> tables.collect { result_file -> tuple("target_qc/${meta.source_cohort}/sample_review", result_file) } })
        .mix(TARGET_QC_REVIEW.out.logs.map { meta, stage_log -> tuple("logs/target_qc/${meta.source_cohort}/sample_review", stage_log) })
        .mix(PREPARE_TARGET.out.prep.map { meta, prep_qc -> tuple("target_prep/${meta.cohort}", prep_qc) })
        .mix(PREPARE_TARGET.out.marker_decisions.map { meta, marker_decisions -> tuple("target_prep/${meta.cohort}", marker_decisions) })
        .mix(PREPARE_TARGET.out.checkpoint
            .filter { meta, _target_dir, _manifest -> ['raw', 'corrected'].contains(meta.input_stage) }
            .map { meta, _target_dir, manifest ->
                def stage = meta.input_stage == 'raw' ? 'corrected' : meta.input_stage
                tuple("checkpoints/${stage}", manifest)
            })
        .mix(TARGET_QC.out.imputation_ready.map { meta, _target_dir, target_qc -> tuple("target_qc/${meta.cohort}", target_qc) })
        .mix(TARGET_QC.out.sample_decisions.map { meta, decisions -> tuple("target_qc/${meta.cohort}", decisions) })
        .mix(TARGET_QC.out.variant_decisions.map { meta, decisions -> tuple("target_qc/${meta.cohort}", decisions) })
        .mix(PARTICIPANT_DECISIONS.out.decisions.map { meta, decisions, _keep -> tuple("target_qc/${meta.cohort}", decisions) })
        .mix(PREPARE_PLINK_REFERENCE.out.summary.map { _meta, summary -> tuple('reference/plink_ct', summary) })
        .mix(PREPARE_PLINK_REFERENCE.out.source_qc.map { _meta, source_qc -> tuple('reference/plink_ct', source_qc) })
        .mix(PREPARE_PLINK_REFERENCE.out.logs.map { _meta, log -> tuple('logs/reference/plink_ct', log) })
        .mix(REFERENCE_ANCESTRY.out.target.map { meta, ancestry -> tuple("target_qc/${meta.cohort}/ancestry", ancestry) })
        .mix(REFERENCE_ANCESTRY.out.reference.map { meta, projection -> tuple("target_qc/${meta.cohort}/ancestry", projection) })
        .mix(REFERENCE_ANCESTRY.out.summary.map { meta, summary -> tuple("target_qc/${meta.cohort}/ancestry", summary) })
        .mix(REFERENCE_ANCESTRY.out.logs.map { meta, log -> tuple("logs/target_qc/${meta.cohort}/ancestry", log) })
        .mix(TARGET_QC.out.imputation_checkpoint.map { meta, _target_dir, manifest ->
            tuple("target/prepared/${meta.cohort}/imputation_ready", manifest)
        })
        .mix(TARGET_QC.out.direct_checkpoint.map { meta, _target_dir, manifest ->
            tuple("target/prepared/${meta.cohort}/direct_ready", manifest)
        })
    score_target_files = TARGET_QC.out.direct_ready.map { meta, target_dir, target_qc ->
        tuple(meta + [scoring_stage: 'direct'], target_dir, target_qc)
    }
    if (target_imputation && stop_after != 'target_qc') {
        imputation_panel = reference_inputs
            .filter { row, _files -> row.reference_type == 'imputation_panel' }
            .first()
        genetic_map = reference_inputs
            .filter { row, _files -> row.reference_type == 'genetic_map' }
            .first()
        beagle_jar = reference_inputs
            .filter { row, _files -> row.reference_type == 'beagle_jar' }
            .first()
        target_impute_chromosomes = TARGET_QC.out.imputation_ready.flatMap { meta, target_dir, target_qc ->
            def available = (1..22).findAll { chromosome ->
                java.nio.file.Files.exists(target_dir.resolve("${meta.cohort}_chr${chromosome}.pgen")) &&
                    java.nio.file.Files.exists(target_dir.resolve("${meta.cohort}_chr${chromosome}.pvar")) &&
                    java.nio.file.Files.exists(target_dir.resolve("${meta.cohort}_chr${chromosome}.psam"))
            }
            if (!available) error "Target '${meta.cohort}' has no chromosome PGEN files for imputation."
            def group_key = groupKey(meta, available.size())
            available.collect { chromosome -> tuple(
                group_key,
                chromosome,
                target_dir.resolve("${meta.cohort}_chr${chromosome}.pgen"),
                target_dir.resolve("${meta.cohort}_chr${chromosome}.pvar"),
                target_dir.resolve("${meta.cohort}_chr${chromosome}.psam"),
                target_qc,
            ) }
        }
        target_impute_input = target_impute_chromosomes
            .combine(imputation_panel)
            .combine(genetic_map)
            .combine(beagle_jar)
            .combine(reference_fasta_source)
        TARGET_IMPUTE_CHROMOSOME(target_impute_input, script_files.target_impute_chromosome, imputation_dr2)
        target_impute_gather = TARGET_IMPUTE_CHROMOSOME.out.chromosomes
            .groupTuple(sort: 'deep')
            .map { group_key, chromosomes, chromosome_dirs, chromosome_manifests, chromosome_qc, chromosome_dr2, chromosome_logs ->
                def order = (0..<chromosomes.size()).toList().sort { index -> chromosomes[index] as int }
                tuple(
                    group_key.getGroupTarget(),
                    order.collect { chromosomes[it] },
                    order.collect { chromosome_dirs[it] },
                    order.collect { chromosome_manifests[it] },
                    order.collect { chromosome_qc[it] },
                    order.collect { chromosome_dr2[it] },
                    order.collect { chromosome_logs[it] },
                )
            }
        ASSEMBLE_TARGET_IMPUTATION(target_impute_gather, script_files.assemble_target_imputation)
        score_target_files = ASSEMBLE_TARGET_IMPUTATION.out.prepared.map { meta, target_dir, target_qc ->
            tuple(meta + [scoring_stage: 'imputed'], target_dir, target_qc)
        }
        result_files = result_files
            .mix(ASSEMBLE_TARGET_IMPUTATION.out.manifest.map { meta, manifest -> tuple("target_imputation/${meta.cohort}", manifest) })
            .mix(ASSEMBLE_TARGET_IMPUTATION.out.qc.map { meta, qc -> tuple("target_imputation/${meta.cohort}", qc) })
            .mix(ASSEMBLE_TARGET_IMPUTATION.out.dr2.map { meta, dr2 -> tuple("target_imputation/${meta.cohort}", dr2) })
            .mix(ASSEMBLE_TARGET_IMPUTATION.out.checkpoint.map { _meta, _target_dir, manifest -> tuple('checkpoints/imputed', manifest) })
            .mix(ASSEMBLE_TARGET_IMPUTATION.out.logs.map { meta, log -> tuple("logs/target_imputation/${meta.cohort}", log) })
        checkpoint_files = checkpoint_files.mix(
            ASSEMBLE_TARGET_IMPUTATION.out.checkpoint.map { _meta, target_dir, _manifest -> tuple('checkpoints/imputed', target_dir) }
        )
    }

    participant_decision_rows = PARTICIPANT_DECISIONS.out.decisions
        .map { meta, decisions, keep -> tuple(meta.cohort, decisions, keep) }
    score_targets = score_target_files
        .map { meta, target_dir, target_qc -> tuple(meta.cohort, meta, target_dir, target_qc) }
        .join(participant_decision_rows, failOnDuplicate: true, failOnMismatch: true)
        .map { _cohort, meta, target_dir, target_qc, decisions, keep ->
            tuple(meta, target_dir, target_qc, decisions, keep)
        }
    direct_score_targets = TARGET_QC.out.direct_ready
        .map { meta, target_dir, target_qc -> tuple(meta.cohort, meta + [score_method: 'plink_ct_direct', scoring_stage: 'direct'], target_dir, target_qc) }
        .join(participant_decision_rows, failOnDuplicate: true, failOnMismatch: true)
        .map { _cohort, meta, target_dir, target_qc, decisions, keep ->
            tuple(meta, target_dir, target_qc, decisions, keep)
        }

    if (run_prs) {
        generation_qc_files = generation_qc_files.mix(
            HARMONISE_GWAS.out.harmonised.map { _meta, _cojo, _clump_input, harmonisation_qc -> harmonisation_qc }
        )
        result_files = result_files
            .mix(HARMONISE_GWAS.out.harmonised.map { meta, cojo, _clump_input, _harmonisation_qc -> tuple("gwas/${meta.trait_id}", cojo) })
            .mix(HARMONISE_GWAS.out.harmonised.map { meta, _cojo, _clump_input, harmonisation_qc -> tuple("qc/gwas/${meta.trait_id}", harmonisation_qc) })
    }

    if (run_prs && methods.contains('plink_ct')) {
        PLINK_REFERENCE_FREQ(plink_reference)

        plink_alignment_input = HARMONISE_GWAS.out.harmonised.combine(plink_reference)
        ALIGN_PLINK_GWAS(plink_alignment_input, script_files.align_plink_gwas)
        PLINK_CLUMP(ALIGN_PLINK_GWAS.out.aligned.combine(plink_reference))
        BUILD_CT_WEIGHTS(PLINK_CLUMP.out.clumped, script_files.ct_weights)

        plink_score_input = score_targets
            .combine(BUILD_CT_WEIGHTS.out.weights)
            .combine(PLINK_REFERENCE_FREQ.out.frequency)
        PLINK_SCORE(plink_score_input)
        PARSE_PLINK_SCORE(PLINK_SCORE.out.raw, script_files.parse_plink, script_files.audit_scoring)
        score_files = score_files.mix(PARSE_PLINK_SCORE.out.scores.map { _target, _gwas, score, _score_qc -> score })
        score_job_records = score_job_records.mix(PARSE_PLINK_SCORE.out.scores.map { target, gwas, _score, _score_qc ->
            [cohort: target.cohort, trait_id: gwas.trait_id, prs_name: gwas.prs_name, method: 'plink_ct']
        })
        generation_qc_files = generation_qc_files.mix(
            ALIGN_PLINK_GWAS.out.qc.map { _meta, alignment_qc -> alignment_qc },
            BUILD_CT_WEIGHTS.out.weights.map { _meta, _weight, weight_qc, _harmonisation_qc, _clump_log -> weight_qc }
        )

        result_files = result_files
            .mix(PLINK_REFERENCE_FREQ.out.frequency.map { _reference, frequency, _reference_log -> tuple('reference/plink_ct', frequency) })
            .mix(PLINK_REFERENCE_FREQ.out.frequency.map { _reference, _frequency, reference_log -> tuple('logs/plink_ct/reference', reference_log) })
            .mix(ALIGN_PLINK_GWAS.out.qc.map { meta, alignment_qc -> tuple("qc/plink_ct/${meta.trait_id}", alignment_qc) })
            .mix(BUILD_CT_WEIGHTS.out.weights.map { meta, weight, _weight_qc, _harmonisation_qc, _clump_log -> tuple("plink_ct/${meta.trait_id}", weight) })
            .mix(BUILD_CT_WEIGHTS.out.weights.map { meta, _weight, weight_qc, _harmonisation_qc, _clump_log -> tuple("qc/plink_ct/${meta.trait_id}", weight_qc) })
            .mix(BUILD_CT_WEIGHTS.out.weights.map { meta, _weight, _weight_qc, _harmonisation_qc, clump_log -> tuple("qc/plink_ct/${meta.trait_id}", clump_log) })
            .mix(PARSE_PLINK_SCORE.out.scores.map { target, gwas, score, _score_qc -> tuple("scores/${target.cohort}/${gwas.trait_id}", score) })
            .mix(PARSE_PLINK_SCORE.out.scores.map { target, gwas, _score, score_qc -> tuple("qc/scores/${target.cohort}/${gwas.trait_id}", score_qc) })
            .mix(PARSE_PLINK_SCORE.out.compatibility.map { target, gwas, audit, _coverage -> tuple("qc/scores/${target.cohort}/${gwas.trait_id}", audit) })
            .mix(PARSE_PLINK_SCORE.out.compatibility.map { target, gwas, _audit, coverage -> tuple("qc/scores/${target.cohort}/${gwas.trait_id}", coverage) })
            .mix(PARSE_PLINK_SCORE.out.logs.map { target, gwas, score_log -> tuple("logs/plink_ct/${target.cohort}/${gwas.trait_id}", score_log) })
        if (target_imputation) {
            plink_direct_input = direct_score_targets
                .combine(BUILD_CT_WEIGHTS.out.weights)
                .combine(PLINK_REFERENCE_FREQ.out.frequency)
            PLINK_DIRECT_SCORE(plink_direct_input)
            PARSE_PLINK_DIRECT_SCORE(PLINK_DIRECT_SCORE.out.raw, script_files.parse_plink, script_files.audit_scoring)

            primary_sensitivity = PARSE_PLINK_SCORE.out.scores
                .map { target, gwas, score, _qc -> tuple("${target.cohort}\t${gwas.trait_id}", target, gwas, score) }
            direct_sensitivity = PARSE_PLINK_DIRECT_SCORE.out.scores
                .map { _target, gwas, score, _qc -> tuple("${_target.cohort}\t${gwas.trait_id}", score) }
            primary_used = PLINK_SCORE.out.raw
                .map { target, gwas, _score, used, _log, _weight, _pvar -> tuple("${target.cohort}\t${gwas.trait_id}", used) }
            direct_used = PLINK_DIRECT_SCORE.out.raw
                .map { target, gwas, _score, used, _log, _weight, _pvar -> tuple("${target.cohort}\t${gwas.trait_id}", used) }
            sensitivity_input = primary_sensitivity
                .join(primary_used, failOnDuplicate: true, failOnMismatch: true)
                .join(direct_sensitivity, failOnDuplicate: true, failOnMismatch: true)
                .join(direct_used, failOnDuplicate: true, failOnMismatch: true)
                .map { _key, target, gwas, primary_score, primary_variants, direct_score, direct_variants ->
                    tuple(target, gwas, primary_score, primary_variants, direct_score, direct_variants)
                }
            COMPARE_DIRECT_SCORE(sensitivity_input, script_files.compare_direct_score)

            result_files = result_files
                .mix(PARSE_PLINK_DIRECT_SCORE.out.scores.map { target, gwas, score, _qc -> tuple("scores/${target.cohort}/${gwas.trait_id}/sensitivity", score) })
                .mix(PARSE_PLINK_DIRECT_SCORE.out.scores.map { target, gwas, _score, qc -> tuple("qc/scores/${target.cohort}/${gwas.trait_id}/sensitivity", qc) })
                .mix(PARSE_PLINK_DIRECT_SCORE.out.compatibility.map { target, gwas, audit, _coverage -> tuple("qc/scores/${target.cohort}/${gwas.trait_id}/sensitivity", audit) })
                .mix(PARSE_PLINK_DIRECT_SCORE.out.compatibility.map { target, gwas, _audit, coverage -> tuple("qc/scores/${target.cohort}/${gwas.trait_id}/sensitivity", coverage) })
                .mix(PLINK_DIRECT_SCORE.out.raw.map { target, gwas, _score, _used, log, _weight, _pvar -> tuple("logs/plink_ct/${target.cohort}/${gwas.trait_id}/sensitivity", log) })
                .mix(COMPARE_DIRECT_SCORE.out.comparison.map { target, gwas, comparison, _qc -> tuple("scores/${target.cohort}/${gwas.trait_id}/sensitivity", comparison) })
                .mix(COMPARE_DIRECT_SCORE.out.comparison.map { target, gwas, _comparison, qc -> tuple("qc/scores/${target.cohort}/${gwas.trait_id}/sensitivity", qc) })
        }
    }

    if (run_prs && methods.contains('sbayesrc')) {
        sbayesrc_ld_source = reference_inputs
            .filter { row, _files -> row.reference_type == 'sbayesrc_ld_source' }
            .first()
        annotation_source = reference_inputs
            .filter { row, _files -> row.reference_type == 'annotation_source' }
            .first()
        PREPARE_SBAYESRC_REFERENCE(
            sbayesrc_ld_source,
            annotation_source,
            script_files.prepare_sbayesrc_reference,
        )
        sbayesrc_ld = PREPARE_SBAYESRC_REFERENCE.out.ld.map { source, ld_dir ->
            tuple(source + [
                reference_type: 'sbayesrc_ld',
                reference_stage: 'prepared',
                source_format: 'directory',
                path: ld_dir.name,
            ], ld_dir)
        }
        annotation = PREPARE_SBAYESRC_REFERENCE.out.annotation.map { source, annotation_file ->
            tuple(source + [
                reference_type: 'annotation',
                reference_stage: 'prepared',
                source_format: 'tsv',
                path: annotation_file.name,
            ], annotation_file)
        }

        tidy_input = HARMONISE_GWAS.out.harmonised.combine(sbayesrc_ld)
        SBAYESRC_TIDY(tidy_input, script_files.sbayesrc)
        impute_input = SBAYESRC_TIDY.out.tidy.combine(sbayesrc_ld)
        SBAYESRC_IMPUTE(impute_input, script_files.sbayesrc)
        model_input = SBAYESRC_IMPUTE.out.imputed
            .combine(sbayesrc_ld)
            .combine(annotation)
        SBAYESRC_MODEL(model_input, script_files.sbayesrc, seed)
        sbayesrc_score_input = score_targets.combine(SBAYESRC_MODEL.out.model)
        SBAYESRC_SCORE(sbayesrc_score_input, script_files.sbayesrc)
        score_files = score_files.mix(SBAYESRC_SCORE.out.scores.map { _target, _gwas, score, _score_qc -> score })
        score_job_records = score_job_records.mix(SBAYESRC_SCORE.out.scores.map { target, gwas, _score, _score_qc ->
            [cohort: target.cohort, trait_id: gwas.trait_id, prs_name: gwas.prs_name, method: 'sbayesrc']
        })
        generation_qc_files = generation_qc_files
            .mix(SBAYESRC_TIDY.out.tidy.map { _meta, _tidy, tidy_qc, _harmonisation_qc -> tidy_qc })
            .mix(SBAYESRC_IMPUTE.out.imputed.map { _meta, _imputed, impute_qc, _tidy_qc, _harmonisation_qc -> impute_qc })
            .mix(SBAYESRC_MODEL.out.model.map { _meta, _weight, _parameter, model_qc, _impute_qc, _tidy_qc, _harmonisation_qc -> model_qc })

        result_files = result_files
            .mix(PREPARE_SBAYESRC_REFERENCE.out.summary.map { summary -> tuple('reference/sbayesrc', summary) })
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
        checkpoint_files = checkpoint_files
            .mix(PREPARE_SBAYESRC_REFERENCE.out.ld.map { _meta, ld_dir -> tuple('reference/sbayesrc/prepared', ld_dir) })
            .mix(PREPARE_SBAYESRC_REFERENCE.out.annotation.map { _meta, annotation_file -> tuple('reference/sbayesrc/prepared', annotation_file) })
    }

    if (run_prs) {
        COMBINE_SCORES(
            score_files.collect(),
            PARTICIPANT_DECISIONS.out.decisions.map { _meta, decisions, _keep -> decisions }.collect(),
            script_files.combine,
        )
        SUMMARISE_GENERATION_QC(
            generation_qc_files.collect(),
            COMBINE_SCORES.out.score_qc,
            script_files.generation_qc,
        )
        result_files = result_files
            .mix(COMBINE_SCORES.out.scores_long.map { result_file -> tuple('scores/combined', result_file) })
            .mix(COMBINE_SCORES.out.scores_wide.map { result_file -> tuple('scores/combined', result_file) })
            .mix(COMBINE_SCORES.out.score_qc.map { result_file -> tuple('qc/scores', result_file) })
            .mix(COMBINE_SCORES.out.concordance.map { result_file -> tuple('qc/scores', result_file) })
            .mix(SUMMARISE_GENERATION_QC.out.variant_flow.map { result_file -> tuple('qc/generation', result_file) })
        if (run_phenotype) {
            phenotype_model_records = phenotype_models
                .splitCsv(header: true, sep: '\t')
                .map { model -> tuple('resolved_models', model.prs_name, model.model_id) }
            phenotype_model_file = phenotype_models
                .map { models -> tuple('resolved_models', models) }
            phenotype_model_jobs = phenotype_model_records
                .combine(phenotype_model_file, by: 0)
                .map { _key, prs_name, model_id, models -> tuple(prs_name, model_id, models) }
            score_model_jobs = score_job_records
                .combine(COMBINE_SCORES.out.scores_long)
                .map { score_job, scores -> tuple(score_job.prs_name, score_job, scores) }
            phenotype_jobs = score_model_jobs
                .combine(phenotype_model_jobs, by: 0)
                .map { _prs_name, score_job, scores, model_id, models -> tuple(score_job, scores, models, model_id) }
            PHENOTYPE_ASSOCIATION(
                phenotype_jobs,
                phenotype_file,
                script_files.association,
                seed,
            )
            COMBINE_PHENOTYPE(
                PHENOTYPE_ASSOCIATION.out.associations.collect(),
                PHENOTYPE_ASSOCIATION.out.fitted_models.collect(),
                PHENOTYPE_ASSOCIATION.out.plot_data.collect(),
                PHENOTYPE_ASSOCIATION.out.permutations.collect(),
                PHENOTYPE_ASSOCIATION.out.influence.collect(),
                PHENOTYPE_ASSOCIATION.out.phenotype_prs.collect(),
                PHENOTYPE_ASSOCIATION.out.phenotype_with_prs.collect(),
                PHENOTYPE_ASSOCIATION.out.participant_level.collect(),
                PHENOTYPE_ASSOCIATION.out.timepoint_completeness.collect(),
                script_files.combine_phenotype,
            )
            result_files = result_files
                .mix(COMBINE_PHENOTYPE.out.associations.map { result_file -> tuple('phenotype', result_file) })
                .mix(COMBINE_PHENOTYPE.out.fitted_models.map { result_file -> tuple('phenotype', result_file) })
                .mix(COMBINE_PHENOTYPE.out.plot_data.map { result_file -> tuple('phenotype', result_file) })
                .mix(COMBINE_PHENOTYPE.out.permutations.map { result_file -> tuple('phenotype', result_file) })
                .mix(COMBINE_PHENOTYPE.out.influence.map { result_file -> tuple('phenotype', result_file) })
                .mix(COMBINE_PHENOTYPE.out.phenotype_prs.map { result_file -> tuple('phenotype', result_file) })
                .mix(COMBINE_PHENOTYPE.out.phenotype_with_prs.map { result_file -> tuple('phenotype', result_file) })
                .mix(COMBINE_PHENOTYPE.out.participant_level.map { result_file -> tuple('phenotype', result_file) })
                .mix(COMBINE_PHENOTYPE.out.timepoint_completeness.map { result_file -> tuple('phenotype', result_file) })
        }
    }

    if (report_enabled && stop_after == 'report') {
        // The report consumes software_versions.yml, so its environment is recorded
        // in a small upstream task to avoid a provenance dependency cycle.
        REPORT_SOFTWARE()
    }

    // Local modules publish the current nf-core (process, tool, version) tuples.
    // Collapse repeated scatter-task records into one deterministic YAML document
    // before it becomes a report input.
    topic_versions = channel.topic('versions')
        .distinct()
        .map { process, tool, version ->
            def process_name = process.substring(process.lastIndexOf(':') + 1)
            [process_name, "  ${tool}: ${version.toString().trim()}"]
        }
        .groupTuple(by: 0)
        .map { process, tool_versions ->
            def unique_tool_versions = tool_versions.unique().sort()
            "${process}:\n${unique_tool_versions.join('\n')}"
        }

    topic_versions_file = topic_versions.collectFile(
        name: 'topic_versions.yml',
        sort: true,
        newLine: true
    )

    COLLECT_VERSIONS(topic_versions_file, script_files.collect_versions)
    result_files = result_files.mix(COLLECT_VERSIONS.out.versions.map { result_file -> tuple('pipeline_info', result_file) })

    if (report_enabled && stop_after == 'report') {
        output_manifest = result_files
            .map { publish_path, result_file -> "${publish_path}\t${result_file.name}" }
            .collectFile(
                name: 'output_files.tsv',
                seed: 'publish_path\tfile_name',
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
        PUBLIC_FIGURES(RENDER_REPORT.out.figures)
        published_files = result_files
            .mix(RENDER_REPORT.out.pages.flatten().map { report_file -> tuple('', report_file) })
            .mix(RENDER_REPORT.out.libraries.map { report_file -> tuple('', report_file) })
            .mix(RENDER_REPORT.out.assets.map { report_file -> tuple('', report_file) })
            .mix(RENDER_REPORT.out.downloads.map { report_file -> tuple('', report_file) })
            .mix(RENDER_REPORT.out.figures.map { report_file -> tuple('', report_file) })
            .mix(RENDER_REPORT.out.provenance.flatten().map { report_file -> tuple('', report_file) })
            .mix(PUBLIC_FIGURES.out.figures.map { report_file -> tuple('', report_file) })
            .mix(checkpoint_files)
    } else {
        published_files = result_files.mix(checkpoint_files)
    }

    emit:
    published_files.map { publish_path, result_file ->
        if (!publish_path) {
            tuple('', result_file)
        } else if (publish_path == 'logs' || publish_path.startsWith('logs/')) {
            tuple(publish_path, result_file)
        } else {
            tuple("data/${publish_path}", result_file)
        }
    }
}
