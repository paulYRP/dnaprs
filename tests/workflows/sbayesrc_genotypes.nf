include { PREPARE_SBAYESRC_TARGETS } from '../../subworkflows/local/prepare_sbayesrc_genotypes/main'

process SBAYESRC_TEST_INPUTS {
    input:
    val meta
    path vcf
    val chromosome_count

    output:
    tuple val(meta), path("${meta.cohort}.imputed"), path('target_qc.tsv'), path('decisions.tsv'), path('keep.tsv')

    script:
    """
    mkdir -p ${meta.cohort}.imputed
    for chromosome in \$(seq 1 '${chromosome_count}'); do
        sed "s/^1\t/\${chromosome}\t/; s/rs1_/rs\${chromosome}_/g; s/ID=1>/ID=\${chromosome}>/" '${vcf}' | gzip -c > ${meta.cohort}.imputed/${meta.cohort}_chr\${chromosome}.vcf.gz
    done
    printf 'cohort\tstatus\n${meta.cohort}\tPASS\n' > target_qc.tsv
    printf 'IID\tscore_eligible\nTEST01\tTRUE\nTEST02\tTRUE\n' > decisions.tsv
    printf 'TEST01\tTEST01\nTEST02\tTEST02\n' > keep.tsv
    """
}

workflow TEST_SBAYESRC_REUSE {
    take:
    vcf
    preparation_script
    chromosome_count

    main:
    cohorts = channel.of([cohort: 'FIRST', scoring_stage: 'imputed'], [cohort: 'SECOND', scoring_stage: 'imputed'])
    SBAYESRC_TEST_INPUTS(cohorts, vcf, chromosome_count)
    PREPARE_SBAYESRC_TARGETS(SBAYESRC_TEST_INPUTS.out, preparation_script)
    scoring_inputs = PREPARE_SBAYESRC_TARGETS.out.prepared.combine(channel.of('TRAIT1', 'TRAIT2'))

    emit:
    prepared = PREPARE_SBAYESRC_TARGETS.out.prepared
    qc = PREPARE_SBAYESRC_TARGETS.out.qc
    logs = PREPARE_SBAYESRC_TARGETS.out.logs
    traits = scoring_inputs
}
