process BUILD_PLINK_COMPATIBILITY {
    tag "${target.cohort}:${reference.reference_id}"
    label 'process_low'
    label 'process_high_memory'
    container 'ghcr.io/paulyrp/dnaprs-analysis:1.0.1'

    input:
    tuple val(target), path(target_pvar, stageAs: 'target/target.pvar'), val(reference), path(reference_pvar, stageAs: 'reference/reference.pvar')
    path index_script

    output:
    tuple val(target), val(reference), path("${target.cohort}.${reference.reference_id}.compatibility.tsv"), emit: index
    tuple val("${task.process}"), val('R'), eval("Rscript -e 'cat(as.character(getRversion()))' 2>/dev/null || printf stub"), emit: versions_r, topic: versions
    tuple val("${task.process}"), val('data.table'), eval("Rscript -e 'cat(as.character(packageVersion(\"data.table\")))' 2>/dev/null || printf stub"), emit: versions_data_table, topic: versions
    script:
    """
    Rscript ${index_script} --target-pvar '${target_pvar}' --reference-pvar '${reference_pvar}' \
        --threads '${task.cpus}' --output-index '${target.cohort}.${reference.reference_id}.compatibility.tsv'
    """

    stub:
    """
    printf 'canonical_key\\tCHR\\tBP\\tcanonical_pair\\ttarget_id\\ttarget_ref\\ttarget_alt\\ttarget_key_count\\treference_id\\treference_ref\\treference_alt\\treference_key_count\\tunique_compatible\\n1:100:A/G\\t1\\t100\\tA/G\\t1:100:A:G\\tA\\tG\\t1\\t1:100:A:G\\tA\\tG\\t1\\tTRUE\\n' > ${target.cohort}.${reference.reference_id}.compatibility.tsv
    """
}
