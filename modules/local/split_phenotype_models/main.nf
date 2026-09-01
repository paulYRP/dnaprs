process SPLIT_PHENOTYPE_MODELS {
    tag 'phenotype models'
    label 'process_single'

    input:
    path phenotype_models

    output:
    path 'models/*.tsv', emit: models, optional: true

    script:
    """
    mkdir -p models
    awk -F '\t' '
        NR == 1 { header = \$0; next }
        NF > 1 {
            count++
            output = "models/" \$1 ".tsv"
            print header > output
            print \$0 >> output
            close(output)
        }
        END { if (count == 0) print header > "models/no_models.tsv" }
    ' ${phenotype_models}
    """

    stub:
    """
    mkdir -p models
    awk -F '\t' '
        NR == 1 { header = \$0; next }
        NF > 1 {
            count++
            output = "models/" \$1 ".tsv"
            print header > output
            print \$0 >> output
            close(output)
        }
        END { if (count == 0) print header > "models/no_models.tsv" }
    ' ${phenotype_models}
    """
}
