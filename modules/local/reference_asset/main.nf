process REFERENCE_ASSET {
    tag row.asset_id
    label 'process_low'

    container 'ghcr.io/paulyrp/dnaprs-imputation:1.1.0'

    publishDir "${params.reference_dir}/${params.reference_bundle}",
        mode: 'copy',
        overwrite: false,
        saveAs: { filename -> filename == row.asset_id ? (row.relative_path ?: null) : null }

    input:
    tuple val(row), path(cached)
    path embedded_fasta_index
    path asset_script

    output:
    tuple val(row), path("${row.asset_id}"), emit: asset
    path "${row.asset_id}.versions.yml", emit: versions

    script:
    """
    bash ${asset_script} \
        '${row.asset_id}' '${row.url}' '${row.checksum_algorithm}' '${row.checksum}' \
        '${row.size}' '${row.etag}' '${cached}' '${embedded_fasta_index}' '${row.asset_id}'
    printf '"${task.process}:${row.asset_id}":\n  curl: %s\n' "\$(curl --version | head -n 1 | awk '{print \$2}')" > ${row.asset_id}.versions.yml
    """

    stub:
    """
    if [[ '${row.asset_id}' == 'cache_complete' ]]; then
        printf 'stub\n' > ${row.asset_id}
    elif [[ '\$(basename ${cached})' != 'empty_input' ]]; then
        cp ${cached} ${row.asset_id}
    else
        cp ${embedded_fasta_index} ${row.asset_id}
    fi
    printf '"${task.process}:${row.asset_id}":\n  curl: stub\n' > ${row.asset_id}.versions.yml
    """
}
