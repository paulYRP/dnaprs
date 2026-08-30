#!/usr/bin/env Rscript

# Read a harmonised GWAS and its fixed PLINK clumping result.
argument <- commandArgs(trailingOnly = TRUE)
name <- sub("^--", "", argument[seq.int(1L, length(argument), 2L)])
value <- argument[seq.int(2L, length(argument), 2L)]
option <- stats::setNames(as.list(value), name)
if (!requireNamespace("data.table", quietly = TRUE)) stop("The data.table package is required.", call. = FALSE)

gwas <- data.table::fread(option[["cojo"]], select = c("SNP", "A1", "b"))
index <- data.table::fread(option[["clumps"]], select = "ID")
if (nrow(index) == 0L || anyDuplicated(index$ID)) stop("PLINK clumping returned no unique index variants.", call. = FALSE)

weight <- merge(index, gwas, by.x = "ID", by.y = "SNP", all.x = TRUE, sort = FALSE)
if (anyNA(weight$A1) || any(!is.finite(weight$b)) || nrow(weight) != nrow(index)) {
  stop("At least one clumped variant could not be mapped to its declared effect allele and effect.", call. = FALSE)
}
data.table::setnames(weight, c("ID", "A1", "b"), c("SNP", "A1", "BETA"))
data.table::fwrite(weight, paste0(option[["trait-id"]], ".plink_ct.weights.tsv"), sep = "\t", quote = FALSE)
data.table::fwrite(
  data.table::data.table(
    trait_id = option[["trait-id"]],
    prs_name = option[["prs-name"]],
    harmonised_variants = nrow(gwas),
    clumped_variants = nrow(weight),
    retained_percent = 100 * nrow(weight) / nrow(gwas),
    status = "PASS"
  ),
  paste0(option[["trait-id"]], ".plink_ct.weight_qc.tsv"),
  sep = "\t"
)
