process TARGET_IMPUTE_CHROMOSOME {
    tag "${group_key.getGroupTarget().cohort}:chr${chromosome}"
    label 'process_high'
    label 'process_plink'

    container 'ghcr.io/paulyrp/dnaprs-imputation:1.1.0'

    input:
    tuple val(group_key), val(chromosome), path(target_pgen), path(target_pvar), path(target_psam), path(target_qc), val(panel), path(panel_files), val(genetic_map), path(genetic_map_files), val(beagle), path(beagle_files)
    path impute_script
    val imputation_dr2

    output:
    tuple val(group_key), val(chromosome), path("*.chr*.imputed"), path("*.chr*.imputation_manifest.tsv"), path("*.chr*.imputation_qc.tsv"), path("*.chr*.imputation_dr2.tsv"), path("*.chr*.target_imputation.log"), emit: chromosomes
    path "*.chr*.versions.yml", emit: versions

    script:
    meta = group_key.getGroupTarget()
    memory_mb = Math.max(1024, (task.memory.toMega() * 0.85) as int)
    """
    BEAGLE_JAR='${beagle.path}' bash ${impute_script} \
        '${meta.cohort}' '${chromosome}' '${target_pgen}' '${target_pvar}' '${target_psam}' \
        '${panel.path}' '${genetic_map.path}' '${imputation_dr2}' '${task.cpus}' '${memory_mb}' '${panel.build}'

    cat > ${meta.cohort}.chr${chromosome}.versions.yml <<-VERSIONS
    "${task.process}:${meta.cohort}:chr${chromosome}":
        beagle: 27Feb25.75f
        bcftools: \$(bcftools --version | head -n 1 | awk '{print \$2}')
        plink2: \$(plink2 --version 2>&1 | head -n 1 | cut -d ' ' -f 2 | sed 's/^v//')
    VERSIONS
    """

    stub:
    meta = group_key.getGroupTarget()
    """
    mkdir -p ${meta.cohort}.chr${chromosome}.imputed
    printf 'stub\n' > ${meta.cohort}.chr${chromosome}.imputed/${meta.cohort}_chr${chromosome}.pgen
    printf '#CHROM\tPOS\tID\tREF\tALT\n${chromosome}\t100\t${chromosome}:100:A:G\tA\tG\n' > ${meta.cohort}.chr${chromosome}.imputed/${meta.cohort}_chr${chromosome}.pvar
    printf '#FID\tIID\nTEST01\tTEST01\nTEST02\tTEST02\n' > ${meta.cohort}.chr${chromosome}.imputed/${meta.cohort}_chr${chromosome}.psam
    printf 'stub\n' > ${meta.cohort}.chr${chromosome}.imputed/${meta.cohort}_chr${chromosome}.vcf.gz
    printf 'stub\n' > ${meta.cohort}.chr${chromosome}.imputed/${meta.cohort}_chr${chromosome}.vcf.gz.tbi
    printf 'cohort\tchromosome\tvcf\tindex\tvcf_sha256\tindex_sha256\tbuild\tdr2_threshold\tstatus\n${meta.cohort}\t${chromosome}\t${meta.cohort}_chr${chromosome}.vcf.gz\t${meta.cohort}_chr${chromosome}.vcf.gz.tbi\tstub\tstub\t${panel.build}\t${imputation_dr2}\tPASS\n' > ${meta.cohort}.chr${chromosome}.imputation_manifest.tsv
    printf 'cohort\tchromosome\tinput_variants\tretained_variants\timputed_variants\tdr2_threshold\tsample_order\tunique_variant_keys\tdosage_range\tstatus\n${meta.cohort}\t${chromosome}\t1\t1\t0\t${imputation_dr2}\tPASS\tPASS\tPASS\tPASS\n' > ${meta.cohort}.chr${chromosome}.imputation_qc.tsv
    printf 'cohort\tchromosome\tdr2_bin\tvariants\n${meta.cohort}\t${chromosome}\t[0.9,1]\t1\n' > ${meta.cohort}.chr${chromosome}.imputation_dr2.tsv
    printf 'Target chromosome imputation stub\n' > ${meta.cohort}.chr${chromosome}.target_imputation.log
    printf '"${task.process}:${meta.cohort}:chr${chromosome}":\n  beagle: stub\n  bcftools: stub\n  plink2: stub\n' > ${meta.cohort}.chr${chromosome}.versions.yml
    """
}
