process HARMONISE_GWAS {
    tag "${meta.trait_id}"
    label 'process_medium'
    label 'process_r'

    container 'ghcr.io/paulyrp/dnaprs-analysis:1.0.0'

    input:
    tuple val(meta), path(gwas_file)
    path harmonise_script

    output:
    tuple val(meta), path("${meta.trait_id}.cojo.ma"), path("${meta.trait_id}.clump.tsv"), path("${meta.trait_id}.harmonisation_qc.tsv"), emit: harmonised
    path 'versions.yml', emit: versions

    script:
    """
    Rscript ${harmonise_script} \
        --input '${gwas_file}' \
        --trait-id '${meta.trait_id}' \
        --prs-name '${meta.prs_name}' \
        --effect-type '${meta.effect_type}' \
        --sample-size '${meta.sample_size}' \
        --snp-col '${meta.snp_col}' \
        --chr-col '${meta.chr_col}' \
        --bp-col '${meta.bp_col}' \
        --effect-allele-col '${meta.effect_allele_col}' \
        --other-allele-col '${meta.other_allele_col}' \
        --beta-col '${meta.beta_col}' \
        --se-col '${meta.se_col}' \
        --p-col '${meta.p_col}' \
        --freq-col '${meta.freq_col ?: ''}' \
        --case-freq-col '${meta.case_freq_col ?: ''}' \
        --control-freq-col '${meta.control_freq_col ?: ''}' \
        --case-n-col '${meta.case_n_col ?: ''}' \
        --control-n-col '${meta.control_n_col ?: ''}' \
        --n-col '${meta.n_col ?: ''}' \
        --info-col '${meta.info_col ?: ''}' \
        --info-min '${meta.info_min ?: ''}' \
        --maf-min '${meta.maf_min ?: '0.01'}' \
        --source-format '${meta.source_format ?: 'auto'}'

    cat > versions.yml <<-VERSIONS
    "${task.process}:${meta.trait_id}":
        R: \$(Rscript -e 'cat(as.character(getRversion()))')
        data.table: \$(Rscript -e "cat(as.character(packageVersion('data.table')))" 2>/dev/null)
    VERSIONS
    """

    stub:
    """
    printf 'SNP\tCHR\tBP\tA1\tA2\tfreq\tb\tse\tp\tN\n1:100:A:G\t1\t100\tG\tA\t0.25\t0.10\t0.03\t0.001\t100000\n' > ${meta.trait_id}.cojo.ma
    printf 'ID\tCHR\tPOS\tA1\tP\n1:100:A:G\t1\t100\tG\t0.001\n' > ${meta.trait_id}.clump.tsv
    printf 'trait_id\tprs_name\tsource_format\tsource_variants\tharmonised_variants\tfiltered_structural\tfiltered_frequency\tfiltered_maf\tfiltered_info\tfiltered_ambiguous\tfiltered_duplicate\tmaf_min\tinfo_min\tduplicated_snp\tstructural_status\n${meta.trait_id}\t${meta.prs_name}\t${meta.source_format ?: 'auto'}\t1\t1\t0\t0\t0\t0\t0\t0\t0.01\t\t0\tPASS\n' > ${meta.trait_id}.harmonisation_qc.tsv
    printf '"${task.process}:${meta.trait_id}":\n  R: stub\n  data.table: stub\n' > versions.yml
    """
}
