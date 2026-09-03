# Input record validation

This process checks the declared target, GWAS, reference, and phenotype inputs before
large PRS tasks start. It does not repeat genotype, GWAS, or phenotype QC that the input
provider has already completed. It checks only the data contract needed to run the
pipeline safely: required fields, supported formats, genome-build consistency, input
existence, genotype identifier uniqueness, repeated-phenotype timepoint rules, and
declared GWAS columns.

The validated target, GWAS, reference, and phenotype-model records are written as
separate TSV files. Source files remain read-only.
