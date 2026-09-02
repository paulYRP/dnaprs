process REFERENCE_ASSET {
    tag "${row.asset_id}"
    label 'process_low'

    container 'ghcr.io/paulyrp/dnaprs-imputation:1.1.0'

    input:
    tuple val(row), path(cached), val(cache_target)
    path embedded_fasta_index
    path asset_script

    output:
    tuple val(row), path("${row.asset_id}"), emit: asset
    tuple val("${task.process}"), val('curl'), eval("command -v curl >/dev/null && curl --version 2>/dev/null | head -n 1 | cut -d ' ' -f 2 || printf stub"), emit: versions_curl, topic: versions
    script:
    """
    bash ${asset_script} \
        '${row.asset_id}' '${row.url}' '${row.checksum_algorithm}' '${row.checksum}' \
        '${row.size}' '${row.etag}' '${cached}' '${cache_target}' '${embedded_fasta_index}' '${row.asset_id}'
    """

    stub:
    """
    if [[ '${row.asset_id}' == 'cache_complete' ]]; then
        printf 'stub\n' > ${row.asset_id}
    elif [[ "\$(basename '${cached}')" != 'empty_input' ]]; then
        cp ${cached} ${row.asset_id}
    else
        cp ${embedded_fasta_index} ${row.asset_id}
    fi
    """
}
