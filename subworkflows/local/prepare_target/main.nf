include { IMPORT_TARGET } from '../../../modules/local/import_target/main'
include { RESOLVE_TARGET_MARKERS } from '../../../modules/local/resolve_target_markers/main'
include { FINALISE_TARGET } from '../../../modules/local/finalise_target/main'

workflow PREPARE_TARGET {
    take:
    target_files
    dbsnp
    reference_fasta
    import_script
    adapter_script
    marker_resolver
    resolve_script
    finalise_script

    main:
    IMPORT_TARGET(target_files, import_script, adapter_script)
    RESOLVE_TARGET_MARKERS(IMPORT_TARGET.out.imported, dbsnp, resolve_script, marker_resolver)
    FINALISE_TARGET(RESOLVE_TARGET_MARKERS.out.resolved, reference_fasta, finalise_script)

    emit:
    prepared = FINALISE_TARGET.out.prepared
    prep = FINALISE_TARGET.out.prep
    marker_decisions = FINALISE_TARGET.out.marker_decisions
    checkpoint = FINALISE_TARGET.out.checkpoint
    duplicate_markers = RESOLVE_TARGET_MARKERS.out.duplicate_markers
}
