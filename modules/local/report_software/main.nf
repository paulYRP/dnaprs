process REPORT_SOFTWARE {
    tag 'report software'
    label 'process_single'

    container 'ghcr.io/paulyrp/dnaprs-report:1.0.0'

    output:
    tuple val("${task.process}"), val('quarto'), eval("quarto --version 2>/dev/null || printf stub"), emit: versions_quarto, topic: versions
    tuple val("${task.process}"), val('R'), eval("Rscript -e 'cat(as.character(getRversion()))' 2>/dev/null || printf stub"), emit: versions_r, topic: versions
    tuple val("${task.process}"), val('ggplot2'), eval("Rscript -e 'cat(as.character(packageVersion(\"ggplot2\")))' 2>/dev/null || printf stub"), emit: versions_ggplot2, topic: versions
    tuple val("${task.process}"), val('jsonlite'), eval("Rscript -e 'cat(as.character(packageVersion(\"jsonlite\")))' 2>/dev/null || printf stub"), emit: versions_jsonlite, topic: versions
    tuple val("${task.process}"), val('knitr'), eval("Rscript -e 'cat(as.character(packageVersion(\"knitr\")))' 2>/dev/null || printf stub"), emit: versions_knitr, topic: versions
    tuple val("${task.process}"), val('openxlsx'), eval("Rscript -e 'cat(as.character(packageVersion(\"openxlsx\")))' 2>/dev/null || printf stub"), emit: versions_openxlsx, topic: versions
    tuple val("${task.process}"), val('rmarkdown'), eval("Rscript -e 'cat(as.character(packageVersion(\"rmarkdown\")))' 2>/dev/null || printf stub"), emit: versions_rmarkdown, topic: versions
    tuple val("${task.process}"), val('scales'), eval("Rscript -e 'cat(as.character(packageVersion(\"scales\")))' 2>/dev/null || printf stub"), emit: versions_scales, topic: versions

    script:
    """
    :
    """

    stub:
    """
    :
    """
}
