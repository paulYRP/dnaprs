process ASSEMBLE_SBAYESRC_GENOTYPES {
    tag "${meta.cohort}"
    label 'process_single'

    container 'docker.io/zhiliz/sbayesrc:0.2.6'

    input:
    tuple val(meta), val(chromosomes), path(chromosome_dirs), path(chromosome_qc), path(target_qc), path(participant_decisions), path(participant_keep)
    path preparation_script

    output:
    tuple val(meta), path("${meta.cohort}.sbayesrc_genotypes"), path(target_qc), path(participant_decisions), path(participant_keep), emit: prepared
    tuple val(meta), path("${meta.cohort}.sbayesrc.genotype_qc.tsv"), emit: qc
    tuple val("${task.process}"), val('R'), eval("Rscript -e 'cat(as.character(getRversion()))' 2>/dev/null || printf stub"), emit: versions_r, topic: versions

    script:
    records = chromosomes.indices.collect { i -> "${chromosomes[i]}\t${chromosome_dirs[i]}\t${chromosome_qc[i]}" }
    """
    printf '%s\n' 'chromosome\tdirectory\tqc' ${records.collect { record -> "'${record}'" }.join(' ')} > chromosomes.tsv
    Rscript ${preparation_script} --action assemble --cohort '${meta.cohort}' \
        --manifest chromosomes.tsv --keep '${participant_keep}' --threads '${task.cpus}'
    """

    stub:
    """
    mkdir -p ${meta.cohort}.sbayesrc_genotypes
    ${chromosome_dirs.collect { directory -> "cp '${directory}'/*.pgen '${directory}'/*.pvar '${directory}'/*.psam ${meta.cohort}.sbayesrc_genotypes/" }.join('\n')}
    head -n 1 '${chromosome_qc[0]}' > ${meta.cohort}.sbayesrc.genotype_qc.tsv
    ${chromosome_qc.collect { qc -> "tail -n +2 '${qc}' >> ${meta.cohort}.sbayesrc.genotype_qc.tsv" }.join('\n')}
    """
}
