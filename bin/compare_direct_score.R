#!/usr/bin/env Rscript

# Compare the primary imputed PLINK score with the direct-genotype sensitivity score.
argument <- commandArgs(trailingOnly = TRUE)
name <- sub("^--", "", argument[seq.int(1L, length(argument), 2L)])
value <- argument[seq.int(2L, length(argument), 2L)]
option <- stats::setNames(as.list(value), name)
if (!requireNamespace("data.table", quietly = TRUE)) stop("The data.table package is required.", call. = FALSE)

primary <- data.table::fread(option[["primary-score"]], colClasses = list(character = c("FID", "IID")))
direct <- data.table::fread(option[["direct-score"]], colClasses = list(character = c("FID", "IID")))
required <- c("cohort", "trait_id", "prs_name", "FID", "IID", "raw_prs")
if (!all(required %in% names(primary)) || !all(required %in% names(direct))) {
  stop("Both sensitivity inputs must follow the dnaprs score contract.", call. = FALSE)
}
key <- c("cohort", "trait_id", "prs_name", "FID", "IID")
comparison <- merge(
  primary[, c(key, "raw_prs"), with = FALSE],
  direct[, c(key, "raw_prs"), with = FALSE],
  by = key, all = TRUE, suffixes = c("_imputed", "_direct"), sort = TRUE
)
if (anyNA(comparison$raw_prs_imputed) || anyNA(comparison$raw_prs_direct)) {
  stop("Direct and imputed scores do not contain the same score-eligible participants.", call. = FALSE)
}
if (nrow(comparison) < 2L || stats::sd(comparison$raw_prs_imputed) <= 0 || stats::sd(comparison$raw_prs_direct) <= 0) {
  stop("At least two varying participant scores are required for sensitivity analysis.", call. = FALSE)
}
comparison[, `:=`(
  imputed_prs_z = as.numeric(scale(raw_prs_imputed)),
  direct_prs_z = as.numeric(scale(raw_prs_direct))
)]
comparison[, z_difference := imputed_prs_z - direct_prs_z]

readUSED <- function(path) unique(as.character(data.table::fread(path, header = FALSE)[[1L]]))
primaryUSED <- readUSED(option[["primary-used"]])
directUSED <- readUSED(option[["direct-used"]])
if (!length(primaryUSED) || !length(directUSED)) stop("A sensitivity score used no variants.", call. = FALSE)

summary <- unique(comparison[, .(cohort, trait_id, prs_name)])
summary[, `:=`(
  participants = nrow(comparison),
  imputed_scoring_variants = length(primaryUSED),
  direct_scoring_variants = length(directUSED),
  shared_variants = length(intersect(primaryUSED, directUSED)),
  imputed_only_variants = length(setdiff(primaryUSED, directUSED)),
  direct_only_variants = length(setdiff(directUSED, primaryUSED)),
  pearson_r = stats::cor(comparison$imputed_prs_z, comparison$direct_prs_z),
  spearman_r = stats::cor(comparison$imputed_prs_z, comparison$direct_prs_z, method = "spearman"),
  mean_z_difference = mean(comparison$z_difference),
  status = "PASS"
)]

prefix <- paste(summary$cohort[[1L]], summary$trait_id[[1L]], "plink_ct", sep = ".")
data.table::fwrite(comparison, paste0(prefix, ".sensitivity.tsv"), sep = "\t")
data.table::fwrite(summary, paste0(prefix, ".sensitivity_qc.tsv"), sep = "\t")
