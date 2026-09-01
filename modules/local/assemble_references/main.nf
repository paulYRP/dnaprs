process ASSEMBLE_REFERENCES {
    tag "${reference_bundle}"
    label 'process_low'
    label 'process_r'

    container 'ghcr.io/paulyrp/dnaprs-analysis:1.0.0'

    input:
    val asset_rows
    path asset_files
    path provided_references, stageAs: 'provided_references.tsv'
    val reference_bundle
    val genome_build
    val cache_root
    path assemble_script

    output:
    path 'references.tsv', emit: references
    path 'reference_bundle', emit: bundle
    path 'reference_receipt.tsv', emit: receipt
    path 'versions.yml', emit: versions

    script:
    assets_json = groovy.json.JsonOutput.toJson(asset_rows).getBytes('UTF-8').encodeBase64().toString()
    """
    Rscript ${assemble_script} \
        --provided '${provided_references}' \
        --assets-json '${assets_json}' \
        --asset-dir . \
        --bundle '${reference_bundle}' \
        --genome-build '${genome_build}' \
        --cache-root '${cache_root}'
    printf '"${task.process}":\n  R: %s\n' "\$(Rscript -e 'cat(as.character(getRversion()))')" > versions.yml
    """

    stub:
    """
    mkdir -p reference_bundle
    cp ${provided_references} references.tsv
    printf 'asset_id\treference_type\tsource_url\tchecksum_algorithm\texpected_checksum\texpected_size\tcache_path\tstatus\n' > reference_receipt.tsv
    printf '"${task.process}":\n  R: stub\n' > versions.yml
    """
}
