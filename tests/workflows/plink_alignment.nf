include { BUILD_PLINK_COMPATIBILITY } from '../../modules/local/build_plink_compatibility/main'
include { ALIGN_PLINK_GWAS } from '../../modules/local/align_plink_gwas/main'

workflow TEST_PLINK_ALIGNMENT {
    take:
    records
    align_script
    index_script

    main:
    records = records instanceof List ? channel.of(records) : records
    index_input = records.map { target, target_dir, _qc, _meta, _cojo, _clump, _harm_qc, reference, reference_dir ->
        tuple(target, target_dir.resolve('ALIGN.pvar'), reference, reference_dir.resolve('eur_reference.pvar'))
    }.distinct()
    gwas = records.map { _target, _dir, _qc, meta, cojo, clump, qc, _reference, _ref_dir -> tuple(meta, cojo, clump, qc) }
    BUILD_PLINK_COMPATIBILITY(index_input, index_script)
    ALIGN_PLINK_GWAS(BUILD_PLINK_COMPATIBILITY.out.index.combine(gwas), align_script)

    emit:
    index = BUILD_PLINK_COMPATIBILITY.out.index
    aligned = ALIGN_PLINK_GWAS.out.aligned
    qc = ALIGN_PLINK_GWAS.out.qc
}
