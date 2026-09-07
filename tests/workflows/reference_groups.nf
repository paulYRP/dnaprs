include { PREPARE_PLINK_REFERENCE } from '../../subworkflows/local/prepare_plink_reference/main'

process MAKE_REFERENCE_FIXTURE {
    input:
    val chromosome_count
    val duplicate_chromosome
    path source_reference

    output:
    path 'panel', emit: panel

    script:
    """
    mkdir -p panel chromosome_sources
    ln -s ../chromosome_sources panel/chromosomes
    for chromosome in \$(seq 1 ${chromosome_count}); do
        if [[ -d '${source_reference}' ]]; then
            cp ${source_reference}/chr\$chromosome.bref3 chromosome_sources/
        else
            sed -E "s/ID=1,/ID=\$chromosome,/;s/^1[[:space:]]/\$chromosome\\t/" ${source_reference} > chromosome_sources/chr\$chromosome.vcf
        fi
    done
    if [[ '${duplicate_chromosome}' == 'true' ]]; then cp chromosome_sources/chr1.vcf chromosome_sources/chr1.duplicate.vcf; fi
    """
}

process STAGE_TEST_UNBREF3 {
    container 'ghcr.io/paulyrp/dnaprs-imputation:1.1.0'

    output:
    path 'unbref3.jar'

    script:
    """
    cp /opt/beagle/unbref3.jar unbref3.jar
    """
}

workflow TEST_REFERENCE_GROUPS {
    take:
    chromosome_count
    duplicate_chromosome
    source_reference
    population
    related
    unbref3
    prepare_script
    assemble_script

    main:
    MAKE_REFERENCE_FIXTURE(chromosome_count, duplicate_chromosome, source_reference)
    source = MAKE_REFERENCE_FIXTURE.out.panel.map { directory -> tuple([reference_id: 'GROUPS', path: 'panel'], [directory]) }
    PREPARE_PLINK_REFERENCE(source, population, related, unbref3, prepare_script, 'GRCh37', assemble_script)

    emit:
    reference = PREPARE_PLINK_REFERENCE.out.reference
    summary = PREPARE_PLINK_REFERENCE.out.summary
    source_qc = PREPARE_PLINK_REFERENCE.out.source_qc
}

workflow TEST_BREF3_REFERENCE_GROUPS {
    take:
    chromosome_count
    source_reference
    population
    related
    prepare_script
    assemble_script

    main:
    STAGE_TEST_UNBREF3()
    unbref3 = STAGE_TEST_UNBREF3.out.map { jar -> tuple([path: 'unbref3.jar'], jar) }
    TEST_REFERENCE_GROUPS(chromosome_count, false, source_reference, population, related, unbref3, prepare_script, assemble_script)

    emit:
    reference = TEST_REFERENCE_GROUPS.out.reference
    summary = TEST_REFERENCE_GROUPS.out.summary
    source_qc = TEST_REFERENCE_GROUPS.out.source_qc
}
