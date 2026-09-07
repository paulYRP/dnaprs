process IMPORT_TARGET {
    tag "${meta.cohort}"
    label 'process_high'
    label 'process_long'
    container 'ghcr.io/paulyrp/dnaprs-imputation:1.1.0'

    input:
    tuple val(meta), path(target_files)
    path import_script
    path adapter_script

    output:
    tuple val(meta), path("${meta.cohort}.imported"), path("${meta.cohort}.initial_marker_decisions.tsv"), emit: imported
    tuple val("${task.process}"), val('plink2'), eval("command -v plink2 >/dev/null && plink2 --version 2>&1 | head -n 1 | cut -d ' ' -f 2 | sed 's/^v//' || printf stub"), emit: versions_plink2, topic: versions
    script:
    """
    PLINK_MEMORY_MB=${Math.max(640, (task.memory.toMega() * 0.8).intValue())} bash ${import_script} '${meta.cohort}' '${meta.format}' '${meta.genotype}' \
        '${meta.sample ?: ''}' '${meta.keep ?: ''}' '${meta.dosage ?: 'DS'}' \
        '${task.cpus}' '${meta.input_stage ?: 'raw'}' '${meta.assay_manifest ?: ''}' \
        '${meta.marker_map ?: ''}' '${adapter_script}'
    """

    stub:
    """
    mkdir -p ${meta.cohort}.imported
    printf 'stub\n' > ${meta.cohort}.imported/${meta.cohort}.pgen
    printf '#CHROM\tPOS\tID\tREF\tALT\n1\t100\t1:100:A:G\tA\tG\n' > ${meta.cohort}.imported/${meta.cohort}.pvar
    printf '#FID\tIID\nTEST01\tTEST01\nTEST02\tTEST02\n' > ${meta.cohort}.imported/${meta.cohort}.psam
    printf 'source_id\tfinal_id\n1:100:A:G\t1:100:A:G\n' > ${meta.cohort}.initial_marker_decisions.tsv
    """
}
