include { PREPARE_PLINK_REFERENCE as PREPARE_REFERENCE_GROUP } from '../../../modules/local/prepare_plink_reference/main'
include { ASSEMBLE_PLINK_REFERENCE } from '../../../modules/local/assemble_plink_reference/main'

workflow PREPARE_PLINK_REFERENCE {
    take:
    panel_source
    population_source
    related_source
    unbref3_source
    prepare_script
    genome_build
    assemble_script

    main:
    grouped_sources = panel_source.flatMap { panel, paths ->
        def groups = ReferenceGroups.balance(paths instanceof List ? paths : [paths])
        def meta = panel + [
            group_count: groups.size(),
            chromosomes: groups.collectMany { group -> group.records*.chromosome }.sort(),
        ]
        groups.collect { group -> tuple(meta, group.id, group.records*.chromosome, group.records*.path) }
    }
    PREPARE_REFERENCE_GROUP(grouped_sources, population_source, related_source, unbref3_source, prepare_script, genome_build)
    gathered = PREPARE_REFERENCE_GROUP.out.prepared
        .map { panel, group_id, directory, qc, log -> tuple(groupKey(panel, panel.group_count), group_id, directory, qc, log) }
        .groupTuple(remainder: true)
        .map { key, ids, directories, qcs, logs ->
            def panel = key.getGroupTarget()
            if (ids.sort(false) != (1..panel.group_count).toList()) error "Reference '${panel.reference_id}' has missing or repeated chromosome groups."
            def order = (0..<ids.size()).toList().sort { index -> ids[index] }
            tuple(panel, order.collect { index -> directories[index] }, order.collect { index -> qcs[index] }, order.collect { index -> logs[index] })
        }
    ASSEMBLE_PLINK_REFERENCE(gathered, population_source, related_source, assemble_script, genome_build)

    emit:
    reference = ASSEMBLE_PLINK_REFERENCE.out.reference
    summary = ASSEMBLE_PLINK_REFERENCE.out.summary
    source_qc = ASSEMBLE_PLINK_REFERENCE.out.source_qc
    logs = ASSEMBLE_PLINK_REFERENCE.out.logs
    identity = ASSEMBLE_PLINK_REFERENCE.out.identity
}
