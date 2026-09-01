process PUBLIC_FIGURES {
    tag 'report figures'
    label 'process_single'

    input:
    path report_figures, stageAs: 'report_figures'

    output:
    path 'figures', emit: figures

    script:
    """
    mkdir -p figures
    cp -R ${report_figures}/. figures/
    """

    stub:
    """
    mkdir -p figures
    cp -R ${report_figures}/. figures/
    """
}
