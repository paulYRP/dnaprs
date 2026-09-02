process ASSEMBLE_TARGET_IMPUTATION {
    tag "${meta.cohort}"
    label 'process_medium'

    container 'ghcr.io/paulyrp/dnaprs-imputation:1.1.0'

    input:
    tuple val(meta), val(chromosomes), path(chromosome_dirs), path(chromosome_manifests), path(chromosome_qc), path(chromosome_dr2), path(chromosome_logs)
    path assembly_script

    output:
    tuple val(meta), path("${meta.cohort}.imputed"), path("${meta.cohort}.imputation_qc.tsv"), emit: prepared
    tuple val(meta), path("${meta.cohort}.imputation_manifest.tsv"), emit: manifest
    tuple val(meta), path("${meta.cohort}.imputation_qc.tsv"), emit: qc
    tuple val(meta), path("${meta.cohort}.imputation_dr2.tsv"), emit: dr2
    tuple val(meta), path("${meta.cohort}.imputed"), path("${meta.cohort}.imputed_target_manifest.tsv"), emit: checkpoint
    tuple val(meta), path("${meta.cohort}.target_imputation.log"), emit: logs
    tuple val("${task.process}"), val('plink2'), eval("command -v plink2 >/dev/null && plink2 --version 2>&1 | head -n 1 | cut -d ' ' -f 2 | sed 's/^v//' || printf stub"), emit: versions_plink2, topic: versions
    script:
    """
    bash ${assembly_script} '${meta.cohort}' '${meta.role}' '${meta.build}' '${meta.ancestry}' '${task.cpus}' '${chromosomes.join(',')}'
    """

    stub:
    """
    mkdir -p ${meta.cohort}.imputed
    cp ${chromosome_dirs[0]}/*.pgen ${meta.cohort}.imputed/${meta.cohort}.pgen
    cp ${chromosome_dirs[0]}/*.pvar ${meta.cohort}.imputed/${meta.cohort}.pvar
    cp ${chromosome_dirs[0]}/*.psam ${meta.cohort}.imputed/${meta.cohort}.psam
    cp -R ${chromosome_dirs[0]}/. ${meta.cohort}.imputed/
    cp ${chromosome_manifests[0]} ${meta.cohort}.imputation_manifest.tsv
    cp ${chromosome_qc[0]} ${meta.cohort}.imputation_qc.tsv
    cp ${chromosome_dr2[0]} ${meta.cohort}.imputation_dr2.tsv
    cp ${chromosome_logs[0]} ${meta.cohort}.target_imputation.log
    printf 'cohort\trole\tsource_format\tgenotype\tsample\tkeep\tbuild\tancestry\tdosage\tinput_stage\tassay_manifest\tmarker_map\n' > ${meta.cohort}.imputed_target_manifest.tsv
    printf '${meta.cohort}\t${meta.role}\tpgen\tcheckpoints/imputed/${meta.cohort}.imputed/${meta.cohort}.pgen\t\t\t${meta.build}\t${meta.ancestry}\tDS\timputed\t\t\n' >> ${meta.cohort}.imputed_target_manifest.tsv
    """
}
