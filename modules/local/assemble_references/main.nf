process ASSEMBLE_REFERENCES {
    tag "${reference_bundle}"
    label 'process_low'

    container 'ghcr.io/paulyrp/dnaprs-analysis:1.0.1'

    input:
    val asset_rows
    path asset_files
    path provided_references, stageAs: 'provided_references.tsv'
    tuple val(provided_paths), path(provided_assets, stageAs: 'provided_reference/asset??/*', arity: '1..*')
    val reference_base
    val reference_bundle
    val genome_build
    val cache_root
    path assemble_script

    output:
    path 'references.tsv', emit: references
    path 'reference_bundle', emit: bundle
    path 'reference_receipt.tsv', emit: receipt
    tuple val("${task.process}"), val('R'), eval("Rscript -e 'cat(as.character(getRversion()))' 2>/dev/null || printf stub"), emit: versions_r, topic: versions
    script:
    assets_json = groovy.json.JsonOutput.toJson(asset_rows).getBytes('UTF-8').encodeBase64().toString()
    provided_map = provided_paths.withIndex().collectEntries { source, index ->
        [(source): provided_assets[index].toString()]
    }
    provided_json = groovy.json.JsonOutput.toJson(provided_map).getBytes('UTF-8').encodeBase64().toString()
    """
    Rscript ${assemble_script} \
        --provided '${provided_references}' \
        --provided-map-json '${provided_json}' \
        --reference-base '${reference_base}' \
        --assets-json '${assets_json}' \
        --asset-dir . \
        --bundle '${reference_bundle}' \
        --genome-build '${genome_build}' \
        --cache-root '${cache_root}'
    """

    stub:
    """
    mkdir -p reference_bundle
    awk -F '\t' -v OFS='\t' -v base='${reference_base}' '
        function rooted(value) {
            return value == "" || substr(value, 1, 1) == "/" || value ~ /^[A-Za-z]:/ ? value : base "/" value
        }
        NR > 1 { \$5 = rooted(\$5); \$6 = rooted(\$6) }
        { print }
    ' ${provided_references} > references.tsv
    printf 'asset_id\treference_type\tsource_url\tchecksum_algorithm\texpected_checksum\tobserved_checksum\texpected_size\tcache_path\tstatus\n' > reference_receipt.tsv
    """
}
