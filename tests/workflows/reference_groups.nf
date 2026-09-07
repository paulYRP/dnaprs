include { PREPARE_PLINK_REFERENCE } from '../../subworkflows/local/prepare_plink_reference/main'

process MAKE_REFERENCE_FIXTURE {
    input:
    val chromosome_count
    val duplicate_chromosome
    path source_vcf

    output:
    path 'panel', emit: panel

    script:
    """
    mkdir -p panel chromosome_sources
    ln -s ../chromosome_sources panel/chromosomes
    for chromosome in \$(seq 1 ${chromosome_count}); do
        sed -E "s/ID=1,/ID=\$chromosome,/;s/^1[[:space:]]/\$chromosome\\t/" ${source_vcf} > chromosome_sources/chr\$chromosome.vcf
    done
    if [[ '${duplicate_chromosome}' == 'true' ]]; then cp chromosome_sources/chr1.vcf chromosome_sources/chr1.duplicate.vcf; fi
    """
}

workflow TEST_REFERENCE_GROUPS {
    take:
    chromosome_count
    duplicate_chromosome
    source_vcf
    population
    related
    unbref3
    prepare_script
    assemble_script

    main:
    MAKE_REFERENCE_FIXTURE(chromosome_count, duplicate_chromosome, source_vcf)
    source = MAKE_REFERENCE_FIXTURE.out.panel.map { directory -> tuple([reference_id: 'GROUPS', path: 'panel'], [directory]) }
    PREPARE_PLINK_REFERENCE(source, population, related, unbref3, prepare_script, 'GRCh37', assemble_script)

    emit:
    reference = PREPARE_PLINK_REFERENCE.out.reference
    summary = PREPARE_PLINK_REFERENCE.out.summary
    source_qc = PREPARE_PLINK_REFERENCE.out.source_qc
}
