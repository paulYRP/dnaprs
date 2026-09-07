process PREPARE_PLINK_REFERENCE {
    tag "${panel.reference_id}:${group_id}"
    label 'process_low'
    container 'ghcr.io/paulyrp/dnaprs-imputation:1.1.0'

    input:
    tuple val(panel), val(group_id), val(chromosomes), path(panel_files, stageAs: 'source/*')
    tuple val(population), path(population_files)
    tuple val(related), path(related_files)
    tuple val(unbref3), path(unbref3_files)
    path prepare_script
    val genome_build

    output:
    tuple val(panel), val(group_id), path("${panel.reference_id}.group${group_id}"), path("${panel.reference_id}.group${group_id}.source_qc.tsv"), path("${panel.reference_id}.group${group_id}.prepare.log"), emit: prepared
    tuple val("${task.process}"), val('bcftools'), eval("command -v bcftools >/dev/null && bcftools --version 2>/dev/null | head -n 1 | cut -d ' ' -f 2 || printf stub"), emit: versions_bcftools, topic: versions
    tuple val("${task.process}"), val('plink2'), eval("command -v plink2 >/dev/null && plink2 --version 2>&1 | head -n 1 | cut -d ' ' -f 2 | sed 's/^v//' || printf stub"), emit: versions_plink2, topic: versions
    tuple val("${task.process}"), val('unbref3'), eval("printf 27Feb25.75f"), emit: versions_unbref3, topic: versions
    script:
    def paths = panel_files instanceof List ? panel_files : [panel_files]
    """
    printf '%s\\n' ${paths.collect { source -> "'$source'" }.join(' ')} > reference_sources.txt
    PLINK_MEMORY_MB=${Math.max(640, (task.memory.toMega() * 0.8).intValue())} UNBREF3_JAR='${unbref3.path}' bash ${prepare_script} \
        reference_sources.txt '${population.path}' '${related.path}' \
        '${panel.reference_id}.group${group_id}' '${task.cpus}' '${genome_build}' '${chromosomes.join(',')}'
    """

    stub:
    def prefix = "${panel.reference_id}.group${group_id}"
    """
    mkdir -p ${prefix}/chromosomes
    printf 'chromosome\\tsource\\tsamples\\tvariants\\tstatus\\tsource_sha256\\n' > ${prefix}.source_qc.tsv
    for chromosome in ${chromosomes.join(' ')}; do
        for population in eur all; do
            printf 'stub\\n' > ${prefix}/chromosomes/\${population}_chr\${chromosome}.pgen
            printf '#CHROM\\tPOS\\tID\\tREF\\tALT\\n%s\\t100\\t%s:100:A:G\\tA\\tG\\n' "\$chromosome" "\$chromosome" > ${prefix}/chromosomes/\${population}_chr\${chromosome}.pvar
            printf '#FID\\tIID\\nTEST01\\tTEST01\\nTEST02\\tTEST02\\n' > ${prefix}/chromosomes/\${population}_chr\${chromosome}.psam
        done
        printf '%s\\tstub\\t2\\t1\\tPASS\\tstub\\n' "\$chromosome" >> ${prefix}.source_qc.tsv
    done
    printf 'Reference group stub\\n' > ${prefix}.prepare.log
    """
}
