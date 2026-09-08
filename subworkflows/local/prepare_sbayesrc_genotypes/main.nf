include { PREPARE_SBAYESRC_GENOTYPES } from '../../../modules/local/prepare_sbayesrc_genotypes/main'
include { ASSEMBLE_SBAYESRC_GENOTYPES } from '../../../modules/local/assemble_sbayesrc_genotypes/main'

workflow PREPARE_SBAYESRC_TARGETS {
    take:
    targets
    preparation_script

    main:
    routes = targets.branch { record ->
        imputed: record[0].scoring_stage == 'imputed'
        other: true
    }
    chromosome_inputs = routes.imputed.flatMap { meta, target_dir, _qc, _decisions, keep ->
        def available = (1..22).findAll { chromosome ->
            java.nio.file.Files.exists(target_dir.resolve("${meta.cohort}_chr${chromosome}.vcf.gz"))
        }
        if (!available) error "Cohort '${meta.cohort}' has no completed chromosome VCFs for SBayesRC genotype preparation."
        def group_key = groupKey(meta, available.size())
        available.collect { chromosome -> tuple(group_key, chromosome, target_dir.resolve("${meta.cohort}_chr${chromosome}.vcf.gz"), keep) }
    }
    PREPARE_SBAYESRC_GENOTYPES(chromosome_inputs, preparation_script)
    gathered = PREPARE_SBAYESRC_GENOTYPES.out.chromosomes
        .groupTuple()
        .map { group_key, chromosomes, directories, qc ->
            def order = (0..<chromosomes.size()).toList().sort { index -> chromosomes[index] as int }
            tuple(group_key.getGroupTarget().cohort,
                order.collect { index -> chromosomes[index] }, order.collect { index -> directories[index] }, order.collect { index -> qc[index] })
        }
        .join(routes.imputed.map { meta, _dir, qc, decisions, keep -> tuple(meta.cohort, meta, qc, decisions, keep) }, failOnDuplicate: true, failOnMismatch: true)
        .map { _cohort, chromosomes, directories, chromosome_qc, meta, qc, decisions, keep ->
            tuple(meta, chromosomes, directories, chromosome_qc, qc, decisions, keep)
        }
    ASSEMBLE_SBAYESRC_GENOTYPES(gathered, preparation_script)

    emit:
    prepared = ASSEMBLE_SBAYESRC_GENOTYPES.out.prepared.mix(routes.other)
    qc = ASSEMBLE_SBAYESRC_GENOTYPES.out.qc
    logs = PREPARE_SBAYESRC_GENOTYPES.out.logs
        .map { group_key, log -> tuple(group_key.getGroupTarget(), log) }
}
