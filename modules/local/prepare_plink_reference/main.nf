process PREPARE_PLINK_REFERENCE {
    tag "${panel.reference_id}"
    label 'process_high'

    container 'ghcr.io/paulyrp/dnaprs-imputation:1.1.0'

    input:
    tuple val(panel), path(panel_files)
    tuple val(population), path(population_files)
    tuple val(related), path(related_files)
    tuple val(unbref3), path(unbref3_files)
    path prepare_script
    val genome_build

    output:
    tuple val(panel), path("${panel.reference_id}.plink_ld"), emit: reference
    tuple val(panel), path("${panel.reference_id}.plink_ld.summary.tsv"), emit: summary
    tuple val(panel), path("${panel.reference_id}.plink_ld.source_qc.tsv"), emit: source_qc
    tuple val(panel), path("${panel.reference_id}.plink_ld.prepare.log"), emit: logs
    path 'versions.yml', emit: versions

    script:
    """
    UNBREF3_JAR='${unbref3.path}' bash ${prepare_script} \
        '${panel.path}' '${population.path}' '${related.path}' \
        '${panel.reference_id}.plink_ld' '${task.cpus}' '${genome_build}'

    cat > versions.yml <<-VERSIONS
    "${task.process}:${panel.reference_id}":
        bcftools: \$(bcftools --version | head -n 1 | awk '{print \$2}')
        plink2: \$(plink2 --version 2>&1 | head -n 1 | cut -d ' ' -f 2 | sed 's/^v//')
        unbref3: 27Feb25.75f
    VERSIONS
    """

    stub:
    """
    mkdir -p ${panel.reference_id}.plink_ld
    printf 'stub\n' > ${panel.reference_id}.plink_ld/eur_reference.pgen
    printf '#CHROM\tPOS\tID\tREF\tALT\n1\t100\t1:100:A:G\tA\tG\n' > ${panel.reference_id}.plink_ld/eur_reference.pvar
    printf '#FID\tIID\nTEST01\tTEST01\nTEST02\tTEST02\n' > ${panel.reference_id}.plink_ld/eur_reference.psam
    cp ${panel.reference_id}.plink_ld/eur_reference.pgen ${panel.reference_id}.plink_ld/all_reference.pgen
    cp ${panel.reference_id}.plink_ld/eur_reference.pvar ${panel.reference_id}.plink_ld/all_reference.pvar
    cp ${panel.reference_id}.plink_ld/eur_reference.psam ${panel.reference_id}.plink_ld/all_reference.psam
    printf 'sample\tpop\tsuper_pop\nTEST01\tGBR\tEUR\nTEST02\tGBR\tEUR\n' > ${panel.reference_id}.plink_ld/population.tsv
    printf 'reference_type\tbuild\tancestry\tchromosomes\tsamples\tvariants\tstatus\nplink_ld\t${genome_build}\tEuropean\t1\t2\t1\tPASS\nancestry_reference\t${genome_build}\tMultiple\t1\t2\t1\tPASS\n' > ${panel.reference_id}.plink_ld.summary.tsv
    printf 'chromosome\tsource\tsamples\tvariants\tstatus\n1\tstub.bref3\t2\t1\tPASS\n' > ${panel.reference_id}.plink_ld.source_qc.tsv
    printf 'PLINK source-reference preparation stub\n' > ${panel.reference_id}.plink_ld.prepare.log
    printf '"${task.process}:${panel.reference_id}":\n  bcftools: stub\n  plink2: stub\n  unbref3: stub\n' > versions.yml
    """
}
