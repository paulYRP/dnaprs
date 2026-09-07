process ALIGN_PLINK_GWAS {
    tag "${target.cohort}:${meta.trait_id}"
    label 'process_low'
    label 'process_high_memory'

    container 'ghcr.io/paulyrp/dnaprs-analysis:1.0.1'

    input:
    tuple val(target), val(reference), path(compatibility_index), val(meta), path(cojo), path(clump_input), path(harmonisation_qc)
    path align_script

    output:
    tuple val(target), val(meta), path("${meta.trait_id}.plink.cojo.ma"), path("${meta.trait_id}.plink.clump.tsv"), path(harmonisation_qc), emit: aligned
    tuple val(target), val(meta), path("${meta.trait_id}.plink.reference_alignment_qc.tsv"), emit: qc
    tuple val("${task.process}"), val('R'), eval("Rscript -e 'cat(as.character(getRversion()))' 2>/dev/null || printf stub"), emit: versions_r, topic: versions
    tuple val("${task.process}"), val('data.table'), eval("Rscript -e 'cat(as.character(packageVersion(\"data.table\")))' 2>/dev/null || printf stub"), emit: versions_data_table, topic: versions
    script:
    """
    Rscript ${align_script} \
        --cojo '${cojo}' \
        --compatibility-index '${compatibility_index}' \
        --threads '${task.cpus}' \
        --cohort '${target.cohort}' \
        --trait-id '${meta.trait_id}' \
        --prs-name '${meta.prs_name}' \
        --output-cojo '${meta.trait_id}.plink.cojo.ma' \
        --output-clump '${meta.trait_id}.plink.clump.tsv' \
        --output-qc '${meta.trait_id}.plink.reference_alignment_qc.tsv'
    """

    stub:
    """
    cp ${cojo} ${meta.trait_id}.plink.cojo.ma
    cp ${clump_input} ${meta.trait_id}.plink.clump.tsv
    printf 'cohort\ttrait_id\tprs_name\tinput_variants\tfiltered_nonpositive_p\tfiltered_palindromic\ttarget_aligned_variants\tfiltered_target_missing_or_mismatched\tfiltered_target_ambiguous\treference_aligned_variants\tfiltered_reference_missing_or_mismatched\tfiltered_reference_ambiguous\tfiltered_reference_duplicate\tcomplemented_alleles\teffect_flips\tstatus\n${target.cohort}\t${meta.trait_id}\t${meta.prs_name}\t1\t0\t0\t1\t0\t0\t1\t0\t0\t0\t0\t0\tPASS\n' > ${meta.trait_id}.plink.reference_alignment_qc.tsv
    """
}
