process PREPARE_SBAYESRC_GENOTYPES {
    tag "${group_key.getGroupTarget().cohort}:chr${chromosome}"
    label 'process_low'

    container 'docker.io/zhiliz/sbayesrc:0.2.6'

    input:
    tuple val(group_key), val(chromosome), path(vcf), path(participant_keep)
    path preparation_script

    output:
    tuple val(group_key), val(chromosome), path('*.chr*.sbayesrc_genotypes'), path('*.chr*.sbayesrc.genotype_qc.tsv'), emit: chromosomes
    tuple val(group_key), path('*.chr*.sbayesrc_genotypes/*.log'), emit: logs
    tuple val("${task.process}"), val('plink2'), eval("command -v plink2 >/dev/null && plink2 --version 2>&1 | head -n 1 | cut -d ' ' -f 2 | sed 's/^v//' || printf stub"), emit: versions_plink2, topic: versions
    tuple val("${task.process}"), val('R'), eval("Rscript -e 'cat(as.character(getRversion()))' 2>/dev/null || printf stub"), emit: versions_r, topic: versions

    script:
    meta = group_key.getGroupTarget()
    memory_mb = Math.max(640, (task.memory.toMega() * 0.8) as int)
    """
    Rscript ${preparation_script} --action prepare \
        --cohort '${meta.cohort}' --chromosome '${chromosome}' \
        --vcf '${vcf}' --keep '${participant_keep}' \
        --threads '${task.cpus}' --memory '${memory_mb}'
    """

    stub:
    meta = group_key.getGroupTarget()
    """
    mkdir -p ${meta.cohort}.chr${chromosome}.sbayesrc_genotypes
    printf 'stub\n' > ${meta.cohort}.chr${chromosome}.sbayesrc_genotypes/${meta.cohort}_chr${chromosome}.pgen
    printf '#CHROM\tPOS\tID\tREF\tALT\n${chromosome}\t100\trs${chromosome}\tA\tG\n' > ${meta.cohort}.chr${chromosome}.sbayesrc_genotypes/${meta.cohort}_chr${chromosome}.pvar
    printf '#FID\tIID\n' > ${meta.cohort}.chr${chromosome}.sbayesrc_genotypes/${meta.cohort}_chr${chromosome}.psam
    cat '${participant_keep}' >> ${meta.cohort}.chr${chromosome}.sbayesrc_genotypes/${meta.cohort}_chr${chromosome}.psam
    printf 'SBayesRC genotype preparation stub\n' > ${meta.cohort}.chr${chromosome}.sbayesrc_genotypes/${meta.cohort}_chr${chromosome}.log
    printf 'cohort\tchromosome\tparticipants\tvariants\tpreserved_ids\tassigned_ids\tduplicate_ids\tparticipant_match\tstatus\n${meta.cohort}\t${chromosome}\t2\t1\t1\t0\t0\tTRUE\tPASS\n' > ${meta.cohort}.chr${chromosome}.sbayesrc.genotype_qc.tsv
    """
}
