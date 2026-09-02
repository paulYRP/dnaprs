process ALIGN_PLINK_GWAS {
    tag "${meta.trait_id}"
    label 'process_medium'

    container 'ghcr.io/paulyrp/dnaprs-analysis:1.0.0'

    input:
    tuple val(meta), path(cojo), path(clump_input), path(harmonisation_qc), val(reference), path(reference_files)
    path align_script

    output:
    tuple val(meta), path("${meta.trait_id}.plink.cojo.ma"), path("${meta.trait_id}.plink.clump.tsv"), path(harmonisation_qc), emit: aligned
    tuple val(meta), path("${meta.trait_id}.plink.reference_alignment_qc.tsv"), emit: qc
    path 'versions.yml', emit: versions

    script:
    """
    Rscript ${align_script} \
        --cojo '${cojo}' \
        --reference-pvar '${reference.path}.pvar' \
        --trait-id '${meta.trait_id}' \
        --prs-name '${meta.prs_name}' \
        --output-cojo '${meta.trait_id}.plink.cojo.ma' \
        --output-clump '${meta.trait_id}.plink.clump.tsv' \
        --output-qc '${meta.trait_id}.plink.reference_alignment_qc.tsv'

    cat > versions.yml <<-VERSIONS
    "${task.process}:${meta.trait_id}":
        R: \$(Rscript -e 'cat(as.character(getRversion()))')
        data.table: \$(Rscript -e "cat(as.character(packageVersion('data.table')))" 2>/dev/null)
    VERSIONS
    """

    stub:
    """
    cp ${cojo} ${meta.trait_id}.plink.cojo.ma
    cp ${clump_input} ${meta.trait_id}.plink.clump.tsv
    printf 'trait_id\tprs_name\tinput_variants\treference_aligned_variants\tfiltered_reference_missing\tfiltered_reference_ambiguous\tfiltered_reference_duplicate\tcomplemented_alleles\tstatus\n${meta.trait_id}\t${meta.prs_name}\t1\t1\t0\t0\t0\t0\tPASS\n' > ${meta.trait_id}.plink.reference_alignment_qc.tsv
    printf '"${task.process}:${meta.trait_id}":\n  R: stub\n  data.table: stub\n' > versions.yml
    """
}
