process TARGET_QC {
    tag "${meta.cohort}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/plink2:2.0.0a.6.9--h9948957_0' :
        'quay.io/biocontainers/plink2:2.0.0a.6.9--h9948957_0' }"

    input:
    tuple val(meta), path(target_dir), path(prep_qc)
    path qc_script
    val imputation_geno
    val direct_geno
    val target_mind
    val target_maf
    val target_hwe

    output:
    tuple val(meta), path("${meta.cohort}.imputation_ready"), path("${meta.cohort}.target_qc.tsv"), emit: imputation_ready
    tuple val(meta), path("${meta.cohort}.direct_ready"), path("${meta.cohort}.target_qc.tsv"), emit: direct_ready
    tuple val(meta), path("${meta.cohort}.sample_decisions.tsv"), emit: sample_decisions
    tuple val(meta), path("${meta.cohort}.variant_decisions.tsv"), emit: variant_decisions
    tuple val(meta), path("${meta.cohort}.imputation_ready"), path("${meta.cohort}.imputation_ready_target.tsv"), emit: imputation_checkpoint
    tuple val(meta), path("${meta.cohort}.direct_ready"), path("${meta.cohort}.direct_ready_target.tsv"), emit: direct_checkpoint
    tuple val("${task.process}"), val('plink2'), eval("command -v plink2 >/dev/null && plink2 --version 2>&1 | head -n 1 | cut -d ' ' -f 2 | sed 's/^v//' || printf stub"), emit: versions_plink2, topic: versions
    script:
    """
    bash ${qc_script} \
        '${meta.cohort}' \
        '${target_dir}' \
        '${meta.input_stage ?: 'raw'}' \
        '${imputation_geno}' \
        '${direct_geno}' \
        '${target_mind}' \
        '${target_maf}' \
        '${target_hwe}' \
        '${task.cpus}'

    printf 'cohort\trole\tsource_format\tgenotype\tsample\tkeep\tbuild\tancestry\tdosage\tinput_stage\tassay_manifest\tmarker_map\\n' > ${meta.cohort}.imputation_ready_target.tsv
    printf '%s\t%s\tpgen\t%s\t\t\t%s\t%s\tDS\tqc_completed\t\t\\n' \
        '${meta.cohort}' '${meta.role}' \
        'data/target/prepared/${meta.cohort}/imputation_ready/${meta.cohort}.pgen' \
        '${meta.build}' '${meta.ancestry}' \
        >> ${meta.cohort}.imputation_ready_target.tsv

    printf 'cohort\trole\tsource_format\tgenotype\tsample\tkeep\tbuild\tancestry\tdosage\tinput_stage\tassay_manifest\tmarker_map\\n' > ${meta.cohort}.direct_ready_target.tsv
    printf '%s\t%s\tpgen\t%s\t\t\t%s\t%s\tDS\tqc_completed\t\t\\n' \
        '${meta.cohort}' '${meta.role}' \
        'data/target/prepared/${meta.cohort}/direct_ready/${meta.cohort}.pgen' \
        '${meta.build}' '${meta.ancestry}' \
        >> ${meta.cohort}.direct_ready_target.tsv
    """

    stub:
    """
    mkdir -p ${meta.cohort}.imputation_ready ${meta.cohort}.direct_ready
    for stage in imputation_ready direct_ready; do
        printf 'stub\n' > ${meta.cohort}.\${stage}/${meta.cohort}.pgen
        printf '#CHROM\tPOS\tID\tREF\tALT\n1\t100\t1:100:A:G\tA\tG\n' > ${meta.cohort}.\${stage}/${meta.cohort}.pvar
        printf '#FID\tIID\nTEST01\tTEST01\nTEST02\tTEST02\n' > ${meta.cohort}.\${stage}/${meta.cohort}.psam
        cp ${meta.cohort}.\${stage}/${meta.cohort}.pgen ${meta.cohort}.\${stage}/${meta.cohort}_chr1.pgen
        cp ${meta.cohort}.\${stage}/${meta.cohort}.pvar ${meta.cohort}.\${stage}/${meta.cohort}_chr1.pvar
        cp ${meta.cohort}.\${stage}/${meta.cohort}.psam ${meta.cohort}.\${stage}/${meta.cohort}_chr1.psam
    done
    printf 'cohort\tinput_stage\tqc_action\tsource_participants\tretained_participants\tfiltered_participants\tsource_variants\timputation_variants\tretained_variants\tfiltered_variants\timputation_variant_missingness\tdirect_variant_missingness\tsample_missingness\tmaf_filter\thwe_filter\tchromosomes\tstatus\n${meta.cohort}\t${meta.input_stage ?: 'raw'}\tFILTERED\t2\t2\t0\t1\t1\t1\t0\t${imputation_geno}\t${direct_geno}\t${target_mind}\t${target_maf}\t${target_hwe}\t1\tPASS\n' > ${meta.cohort}.target_qc.tsv
    printf 'cohort\tFID\tIID\tmissingness\tdecision\treason\n${meta.cohort}\tTEST01\tTEST01\t0\tRETAIN\tWithin threshold\n${meta.cohort}\tTEST02\tTEST02\t0\tRETAIN\tWithin threshold\n' > ${meta.cohort}.sample_decisions.tsv
    printf 'cohort\tID\tchromosome\tposition\tref\talt\tmissingness\tmaf\thwe_midp\timputation_decision\tdirect_decision\tdecision\treason\n${meta.cohort}\t1:100:A:G\t1\t100\tA\tG\t0\t0.5\t1\tRETAIN\tRETAIN\tRETAIN\tWithin thresholds\n' > ${meta.cohort}.variant_decisions.tsv
    printf 'cohort\trole\tsource_format\tgenotype\tsample\tkeep\tbuild\tancestry\tdosage\tinput_stage\tassay_manifest\tmarker_map\n' > ${meta.cohort}.imputation_ready_target.tsv
    printf '${meta.cohort}\t${meta.role}\tpgen\tdata/target/prepared/${meta.cohort}/imputation_ready/${meta.cohort}.pgen\t\t\t${meta.build}\t${meta.ancestry}\tDS\tqc_completed\t\t\n' >> ${meta.cohort}.imputation_ready_target.tsv
    printf 'cohort\trole\tsource_format\tgenotype\tsample\tkeep\tbuild\tancestry\tdosage\tinput_stage\tassay_manifest\tmarker_map\n' > ${meta.cohort}.direct_ready_target.tsv
    printf '${meta.cohort}\t${meta.role}\tpgen\tdata/target/prepared/${meta.cohort}/direct_ready/${meta.cohort}.pgen\t\t\t${meta.build}\t${meta.ancestry}\tDS\tqc_completed\t\t\n' >> ${meta.cohort}.direct_ready_target.tsv
    """
}
