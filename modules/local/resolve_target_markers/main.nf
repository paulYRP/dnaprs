process RESOLVE_TARGET_MARKERS {
    tag "${meta.cohort}"
    label 'process_high'
    label 'process_long'
    container 'ghcr.io/paulyrp/dnaprs-imputation:1.1.0'

    input:
    tuple val(meta), path(imported), path(initial)
    tuple val(dbsnp), path(dbsnp_files)
    path resolve_script
    path marker_resolver

    output:
    tuple val(meta), path(imported), path("${meta.cohort}.marker_decisions.tsv"), path("${meta.cohort}.retained_markers.txt"), path("${meta.cohort}.rename_markers.tsv"), emit: resolved
    tuple val(meta), path("${meta.cohort}.duplicate_markers.txt"), emit: duplicate_markers
    tuple val("${task.process}"), val('bcftools'), eval("command -v bcftools >/dev/null && bcftools --version 2>/dev/null | head -n 1 | cut -d ' ' -f 2 || printf stub"), emit: versions_bcftools, topic: versions
    tuple val("${task.process}"), val('plink2'), eval("command -v plink2 >/dev/null && plink2 --version 2>&1 | head -n 1 | cut -d ' ' -f 2 | sed 's/^v//' || printf stub"), emit: versions_plink2, topic: versions
    tuple val("${task.process}"), val('R'), eval("Rscript -e 'cat(as.character(getRversion()))' 2>/dev/null || printf stub"), emit: versions_r, topic: versions
    tuple val("${task.process}"), val('data.table'), eval("Rscript -e 'cat(as.character(packageVersion(\"data.table\")))' 2>/dev/null || printf stub"), emit: versions_data_table, topic: versions
    script:
    """
    PLINK_MEMORY_MB=${Math.max(640, (task.memory.toMega() * 0.8).intValue())} bash ${resolve_script} '${meta.cohort}' '${imported}' '${initial}' \
        '${meta.input_stage ?: 'raw'}' '${dbsnp.path}' '${marker_resolver}' '${task.cpus}'
    """

    stub:
    """
    printf 'source_id\tfinal_id\tsource_chr\tsource_pos\tfinal_chr\tfinal_pos\tsource_ref\tsource_alt\tfinal_ref\tfinal_alt\tdecision\treason\tcall_count\tcall_rate\tassay_ref_a\tassay_ref_b\tassay_pair\tassay_status\n1:100:A:G\t1:100:A:G\t1\t100\t1\t100\tA\tG\tA\tG\tINHERITED\tStub target\t2\t1\tA\tG\tA/G\tNOT_SUPPLIED\n' > ${meta.cohort}.marker_decisions.tsv
    printf '1:100:A:G\n' > ${meta.cohort}.retained_markers.txt
    printf '1:100:A:G\t1:100:A:G\n' > ${meta.cohort}.rename_markers.tsv
    touch ${meta.cohort}.duplicate_markers.txt
    """
}
