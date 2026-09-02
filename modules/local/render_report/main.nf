process RENDER_REPORT {
    tag 'dnaprs HTML report'
    label 'process_low'

    container 'ghcr.io/paulyrp/dnaprs-report:1.0.0'

    input:
    path report_files, stageAs: 'report_inputs/*'
    path output_manifest
    path report_source

    output:
    path 'reports/*.html', emit: pages
    path 'reports/site_libs', emit: libraries
    path 'reports/assets', emit: assets
    path 'reports/downloads', optional: true, emit: downloads
    path 'reports/figures', optional: true, emit: figures
    path 'reports/provenance/*', emit: provenance

    script:
    """
    report_root=\$(pwd)
    mkdir -p report_project reports
    cp -R '${report_source}/.' report_project/

    export DNAPRS_REPORT_INPUTS="\$report_root/report_inputs"
    export DNAPRS_OUTPUT_MANIFEST="\$report_root/${output_manifest}"
    export QUARTO_VERSION=\$(quarto --version)
    export HOME="\$report_root"
    export XDG_CACHE_HOME="\$report_root/.cache"
    export QUARTO_CACHE_DIR="\$report_root/.cache/quarto"
    mkdir -p "\$QUARTO_CACHE_DIR"
    cd report_project
    Rscript prepare-report.R
    quarto render .
    cp -R _site/. ../reports/
    """

    stub:
    """
    mkdir -p reports/site_libs reports/assets reports/downloads reports/figures/tiff reports/figures/png reports/figures/jpeg reports/provenance
    printf '<html><body><h1>dnaprs stub report</h1></body></html>\n' > reports/index.html
    printf 'publish_path\tfile_name\n' > reports/provenance/output_files.tsv
    printf 'result_file\tcolumn\tdata_type\tmeaning\tmissing_value\n' > reports/provenance/data_dictionary.tsv
    printf 'figure_id\tpage\tsection\ttitle\tdescription\tinspection\twidth_in\theight_in\tdpi\ttiff\tpng\tjpeg\n' > reports/provenance/figure_manifest.tsv
    """
}
