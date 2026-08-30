process COLLECT_VERSIONS {
    tag 'software versions'
    label 'process_low'

    input:
    path version_files, stageAs: 'versions/versions????.yml'

    output:
    path 'software_versions.yml', emit: versions

    script:
    """
    awk '
        /^[^[:space:]].*:\$/ {
            if (seen[\$0]++) {
                print "Duplicate software-version key: " \$0 > "/dev/stderr"
                exit 1
            }
        }
        { print }
    ' versions/*.yml > software_versions.yml
    """

    stub:
    """
    printf 'DNAPRS_STUB:\n  version: stub\n' > software_versions.yml
    """
}
