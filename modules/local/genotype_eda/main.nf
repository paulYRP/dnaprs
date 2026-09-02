process GENOTYPE_EDA {
    tag "${meta.cohort}"
    label 'process_high'

    container 'ghcr.io/paulyrp/dnaprs-plink2:2.0.0-a.6.12'

    input:
    tuple val(meta), path(target_files)
    path genotype_eda_script
    path adapter_script

    output:
    tuple val(meta), path("${meta.cohort}.*.tsv"), emit: tables
    tuple val(meta), path("${meta.cohort}.genotype_eda.log"), emit: logs
    tuple val("${task.process}"), val('plink2'), eval("command -v plink2 >/dev/null && plink2 --version 2>&1 | head -n 1 | cut -d ' ' -f 2 | sed 's/^v//' || printf stub"), emit: versions_plink2, topic: versions
    script:
    """
    bash ${genotype_eda_script} \
        '${meta.cohort}' \
        '${meta.format}' \
        '${meta.genotype}' \
        '${meta.sample ?: ''}' \
        '${meta.dosage ?: 'DS'}' \
        '${task.cpus}' \
        '${meta.input_stage ?: 'raw'}' \
        '${meta.role}' \
        '${meta.assay_manifest ?: ''}' \
        '${adapter_script}'
    """

    stub:
    """
    printf 'cohort\trole\tinput_stage\tformat\tparticipants\tvariants\tchromosomes\tautosomal_variants\tx_variants\ty_variants\tmitochondrial_variants\tunplaced_or_nonstandard_variants\trecorded_sex_participants\tphenotype_participants\tduplicated_participant_identifiers\tduplicated_variant_identifier_groups\tduplicated_variant_identifier_records\treview_items\tstatus\n${meta.cohort}\t${meta.role}\t${meta.input_stage ?: 'qc_completed'}\t${meta.format}\t2\t2\t1\t2\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\tPASS\n' > ${meta.cohort}.genotype_eda_summary.tsv
    printf 'cohort\tchromosome\tcategory\tvariants\n${meta.cohort}\t1\tautosome\t2\n' > ${meta.cohort}.chromosome_counts.tsv
    printf 'cohort\tchromosome\tbin_start\tbin_end\tvariants\n${meta.cohort}\t1\t1\t1000000\t2\n' > ${meta.cohort}.marker_density.tsv
    printf 'cohort\tidentifier_class\tvariants\tpercent\n${meta.cohort}\tcoordinate\t2\t100\n' > ${meta.cohort}.identifier_classes.tsv
    printf 'cohort\tallele_state\tvariants\tpercent\n${meta.cohort}\tacgt_snp\t2\t100\n' > ${meta.cohort}.allele_states.tsv
    printf 'cohort\tidentifier\toccurrences\n' > ${meta.cohort}.duplicate_identifiers.tsv
    printf '#FID\tIID\tMISSING_CT\tOBS_CT\tF_MISS\nTEST01\tTEST01\t0\t2\t0\nTEST02\tTEST02\t0\t2\t0\n' > ${meta.cohort}.sample_missingness.tsv
    printf '#CHROM\tID\tMISSING_CT\tOBS_CT\tF_MISS\n1\t1:100:A:G\t0\t2\t0\n1\t1:200:C:T\t0\t2\t0\n' > ${meta.cohort}.variant_missingness.tsv
    printf 'cohort\tmissingness_bin\tvariants\tpercent\n${meta.cohort}\t0\t2\t100\n' > ${meta.cohort}.variant_missingness_bins.tsv
    printf '#CHROM\tID\tREF\tALT\tALT_FREQS\tOBS_CT\n1\t1:100:A:G\tA\tG\t0.5\t4\n1\t1:200:C:T\tC\tT\t0.5\t4\n' > ${meta.cohort}.allele_frequency.tsv
    printf 'cohort\tmaf_bin\tvariants\tpercent\n${meta.cohort}\t(0.20,0.50]\t2\t100\n' > ${meta.cohort}.allele_frequency_bins.tsv
    printf 'cohort\tFID\tIID\tobserved_homozygotes\texpected_homozygotes\tobservations\tinbreeding_coefficient\theterozygosity_rate\tmissingness\theterozygosity_z\tstatus\n${meta.cohort}\tTEST01\tTEST01\t1\t1\t2\t0\t0.5\t0\t0\tPASS\n${meta.cohort}\tTEST02\tTEST02\t1\t1\t2\t0\t0.5\t0\t0\tPASS\n' > ${meta.cohort}.heterozygosity.tsv
    printf 'cohort\tFID\tIID\trecorded_sex\tgenetic_sex\tx_inbreeding_coefficient\tstatus\treason\n' > ${meta.cohort}.sex_check.tsv
    printf 'cohort\tFID1\tIID1\tFID2\tIID2\tvariants\tkinship\trelationship_category\treview_status\n' > ${meta.cohort}.relatedness.tsv
    printf 'cohort\tkinship_bin\tpairs\tpercent\n' > ${meta.cohort}.relatedness_bins.tsv
    printf 'cohort\tFID\tIID\tPC1\tPC2\tPC3\tPC4\tPC5\tPC6\tPC7\tPC8\tPC9\tPC10\n' > ${meta.cohort}.internal_pca.tsv
    printf 'cohort\tcomponent\teigenvalue\tpercent_of_reported_eigenvalues\n' > ${meta.cohort}.pca_eigenvalues.tsv
    printf 'cohort\tcheck\tstatus\tvalue\treason\n${meta.cohort}\tformat_import\tPASS\t2 participants; 2 variants\tThe supplied target was imported without changing the source files.\n${meta.cohort}\treported_sex\tNOT_RUN\t0 recorded; 0 X variants\tRecorded sex or X-chromosome variants were unavailable.\n' > ${meta.cohort}.genotype_eda_checks.tsv
    printf 'Genotype EDA stub completed.\n' > ${meta.cohort}.genotype_eda.log
    """
}
