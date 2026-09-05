process COLLECT_VERSIONS {
    tag 'software versions'
    label 'process_low'

    container 'ghcr.io/paulyrp/dnaprs-analysis:1.0.1'

    input:
    path version_files, stageAs: 'versions/versions????.yml'
    path collect_script

    output:
    path 'software_versions.yml', emit: versions

    script:
    """
    Rscript ${collect_script} --input-dir versions --output software_versions.yml
    """

    stub:
    """
    printf 'DNAPRS_STUB:\n  version: stub\n' > software_versions.yml
    """
}
