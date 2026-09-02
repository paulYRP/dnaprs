process PREPARE_SBAYESRC_REFERENCE {
    tag "${ld_source.reference_id}"
    label 'process_high'

    container 'ghcr.io/paulyrp/dnaprs-imputation:1.1.0'

    input:
    tuple val(ld_source), path(ld_files)
    tuple val(annotation_source), path(annotation_files)
    path prepare_script

    output:
    tuple val(ld_source), path('prepared_sbayesrc/ld'), emit: ld
    tuple val(annotation_source), path('prepared_sbayesrc/annotation/annotation.txt'), emit: annotation
    path 'prepared_sbayesrc.summary.tsv', emit: summary
    tuple val("${task.process}"), val('unzip'), eval("command -v unzip >/dev/null && unzip -v 2>/dev/null | head -n 1 | cut -d ' ' -f 2 || printf stub"), emit: versions_unzip, topic: versions
    script:
    """
    bash ${prepare_script} '${ld_source.path}' '${annotation_source.path}' prepared_sbayesrc
    """

    stub:
    """
    mkdir -p prepared_sbayesrc/ld prepared_sbayesrc/annotation
    printf 'stub\n' > prepared_sbayesrc/ld/block1.eigen.bin
    printf 'SNP\tA1\tA2\n1:100:A:G\tG\tA\n' > prepared_sbayesrc/annotation/annotation.txt
    printf 'reference_type\tfiles_or_rows\tstatus\nsbayesrc_ld\t1\tPASS\nannotation\t1\tPASS\n' > prepared_sbayesrc.summary.tsv
    """
}
