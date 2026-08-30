process PREPARE_TARGET {
    tag "${meta.cohort}"
    label 'process_high'
    label 'process_plink'

    input:
    val meta
    path prepare_script

    output:
    tuple val(meta), path("${meta.cohort}"), path("${meta.cohort}.target_qc.tsv"), emit: prepared
    path 'versions.yml', emit: versions

    script:
    """
    bash ${prepare_script} \
        '${meta.cohort}' \
        '${meta.format}' \
        '${meta.genotype}' \
        '${meta.sample ?: ''}' \
        '${meta.keep ?: ''}' \
        '${meta.dosage ?: 'DS'}' \
        '${task.cpus}'

    cat > versions.yml <<-VERSIONS
    "${task.process}:${meta.cohort}":
        plink2: \$(plink2 --version 2>&1 | head -n 1 | cut -d ' ' -f 2 | sed 's/^v//')
    VERSIONS
    """

    stub:
    """
    mkdir -p ${meta.cohort}
    printf 'stub\n' > ${meta.cohort}/${meta.cohort}.pgen
    printf '#CHROM\tPOS\tID\tREF\tALT\n1\t100\t1:100:A:G\tA\tG\n' > ${meta.cohort}/${meta.cohort}.pvar
    printf '#FID\tIID\nTEST01\tTEST01\nTEST02\tTEST02\n' > ${meta.cohort}/${meta.cohort}.psam
    cp ${meta.cohort}/${meta.cohort}.pgen ${meta.cohort}/${meta.cohort}_chr1.pgen
    cp ${meta.cohort}/${meta.cohort}.pvar ${meta.cohort}/${meta.cohort}_chr1.pvar
    cp ${meta.cohort}/${meta.cohort}.psam ${meta.cohort}/${meta.cohort}_chr1.psam
    printf 'cohort\tparticipants\tvariants\tchromosomes\tstatus\n${meta.cohort}\t2\t1\t1\tPASS\n' > ${meta.cohort}.target_qc.tsv
    printf '"${task.process}:${meta.cohort}":\n  plink2: stub\n' > versions.yml
    """
}
