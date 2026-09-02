process PREPARE_TARGET {
    tag "${meta.cohort}"
    label 'process_high'

    container 'ghcr.io/paulyrp/dnaprs-imputation:1.1.0'

    input:
    tuple val(meta), path(target_files)
    tuple val(dbsnp), path(dbsnp_files)
    tuple val(reference_fasta), path(reference_fasta_files)
    path prepare_script
    path adapter_script
    path marker_resolver

    output:
    tuple val(meta), path("${meta.cohort}"), path("${meta.cohort}.target_prep_summary.tsv"), emit: prepared
    tuple val(meta), path("${meta.cohort}.target_prep_summary.tsv"), emit: prep
    tuple val(meta), path("${meta.cohort}.marker_decisions.tsv"), emit: marker_decisions
    tuple val(meta), path("${meta.cohort}"), path("${meta.cohort}.corrected_target_manifest.tsv"), emit: checkpoint
    tuple val("${task.process}"), val('bcftools'), eval("command -v bcftools >/dev/null && bcftools --version 2>/dev/null | head -n 1 | cut -d ' ' -f 2 || printf stub"), emit: versions_bcftools, topic: versions
    tuple val("${task.process}"), val('plink2'), eval("command -v plink2 >/dev/null && plink2 --version 2>&1 | head -n 1 | cut -d ' ' -f 2 | sed 's/^v//' || printf stub"), emit: versions_plink2, topic: versions
    tuple val("${task.process}"), val('R'), eval("Rscript -e 'cat(as.character(getRversion()))' 2>/dev/null || printf stub"), emit: versions_r, topic: versions
    script:
    checkpoint_stage = meta.input_stage == 'raw' ? 'corrected' : meta.input_stage
    """
    bash ${prepare_script} \
        '${meta.cohort}' \
        '${meta.format}' \
        '${meta.genotype}' \
        '${meta.sample ?: ''}' \
        '${meta.keep ?: ''}' \
        '${meta.dosage ?: 'DS'}' \
        '${task.cpus}' \
        '${meta.input_stage ?: 'raw'}' \
        '${meta.assay_manifest ?: ''}' \
        '${meta.marker_map ?: ''}' \
        '${adapter_script}' \
        '${marker_resolver}' \
        '${dbsnp.path}' \
        '${reference_fasta.path}'

    printf 'cohort\trole\tsource_format\tgenotype\tsample\tkeep\tbuild\tancestry\tdosage\tinput_stage\tassay_manifest\tmarker_map\\n' > ${meta.cohort}.corrected_target_manifest.tsv
    printf '%s\t%s\tpgen\t%s\t\t\t%s\t%s\tDS\t%s\t\t\\n' \
        '${meta.cohort}' '${meta.role}' \
        'checkpoints/${checkpoint_stage}/${meta.cohort}/${meta.cohort}.pgen' \
        '${meta.build}' '${meta.ancestry}' '${checkpoint_stage}' \
        >> ${meta.cohort}.corrected_target_manifest.tsv
    """

    stub:
    checkpoint_stage = meta.input_stage == 'raw' ? 'corrected' : meta.input_stage
    """
    mkdir -p ${meta.cohort}
    printf 'stub\n' > ${meta.cohort}/${meta.cohort}.pgen
    printf '#CHROM\tPOS\tID\tREF\tALT\n1\t100\t1:100:A:G\tA\tG\n' > ${meta.cohort}/${meta.cohort}.pvar
    printf '#FID\tIID\nTEST01\tTEST01\nTEST02\tTEST02\n' > ${meta.cohort}/${meta.cohort}.psam
    cp ${meta.cohort}/${meta.cohort}.pgen ${meta.cohort}/${meta.cohort}_chr1.pgen
    cp ${meta.cohort}/${meta.cohort}.pvar ${meta.cohort}/${meta.cohort}_chr1.pvar
    cp ${meta.cohort}/${meta.cohort}.psam ${meta.cohort}/${meta.cohort}_chr1.psam
    printf 'cohort\tinput_stage\tstep\tparticipants\tvariants\tstatus\n${meta.cohort}\t${meta.input_stage ?: 'qc_completed'}\tNormalised to PGEN\t2\t1\tPASS\n' > ${meta.cohort}.target_prep_summary.tsv
    printf 'cohort\tinput_stage\tparticipants\tvariants\tchromosomes\tstatus\n${meta.cohort}\t${meta.input_stage ?: 'qc_completed'}\t2\t1\t1\tPASS\n' > ${meta.cohort}.target_qc.tsv
    printf 'source_id\tfinal_id\tsource_chr\tsource_pos\tfinal_chr\tfinal_pos\tsource_ref\tsource_alt\tfinal_ref\tfinal_alt\tdecision\treason\n1:100:A:G\t1:100:A:G\t1\t100\t1\t100\tA\tG\tA\tG\tINHERITED\tStub target\n' > ${meta.cohort}.marker_decisions.tsv
    printf 'cohort\trole\tsource_format\tgenotype\tsample\tkeep\tbuild\tancestry\tdosage\tinput_stage\tassay_manifest\tmarker_map\n' > ${meta.cohort}.corrected_target_manifest.tsv
    printf '${meta.cohort}\t${meta.role}\tpgen\tcheckpoints/${checkpoint_stage}/${meta.cohort}/${meta.cohort}.pgen\t\t\t${meta.build}\t${meta.ancestry}\tDS\t${checkpoint_stage}\t\t\n' >> ${meta.cohort}.corrected_target_manifest.tsv
    """
}
