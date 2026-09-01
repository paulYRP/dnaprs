#!/usr/bin/env Rscript

argument <- commandArgs(trailingOnly = TRUE)
name <- sub("^--", "", argument[seq.int(1L, length(argument), 2L)])
value <- argument[seq.int(2L, length(argument), 2L)]
option <- stats::setNames(as.list(value), name)
if (!requireNamespace("data.table", quietly = TRUE)) stop("The data.table package is required.", call. = FALSE)

inputDIR <- option[["input-dir"]]

combineTABLE <- function(suffix, output, orderCOLUMN) {
  paths <- sort(list.files(inputDIR, pattern = paste0("\\.", suffix, "$"), full.names = TRUE))
  if (length(paths) == 0L) stop(sprintf("No phenotype task output matched *.%s", suffix), call. = FALSE)
  tables <- lapply(paths, data.table::fread)
  combined <- data.table::rbindlist(tables, use.names = TRUE, fill = TRUE)
  ordering <- intersect(orderCOLUMN, names(combined))
  if (nrow(combined) > 0L && length(ordering) > 0L) data.table::setorderv(combined, ordering)
  data.table::fwrite(combined, output, sep = "\t", na = "NA")
}

combineTABLE(
  "phenotype_associations.tsv", "phenotype_associations.tsv",
  c("model_id", "cohort", "trait_id", "method")
)
combineTABLE(
  "phenotype_models_fitted.tsv", "phenotype_models_fitted.tsv",
  c("model_id", "cohort", "method")
)
combineTABLE(
  "phenotype_plot_data.tsv", "phenotype_plot_data.tsv",
  c("model_id", "cohort", "method", "IID")
)
combineTABLE(
  "phenotype_permutations.tsv", "phenotype_permutations.tsv",
  c("model_id", "cohort", "method", "permutation_id")
)
combineTABLE(
  "phenotype_influence.tsv", "phenotype_influence.tsv",
  c("model_id", "cohort", "method", "IID")
)

phenotypePATH <- sort(list.files(inputDIR, pattern = "\\.phenoPRS.csv$", full.names = TRUE))
if (length(phenotypePATH) == 0L) stop("No phenotype task produced phenoPRS.csv.", call. = FALSE)
digest <- unname(tools::md5sum(phenotypePATH))
if (length(unique(digest)) != 1L) {
  stop("Parallel phenotype tasks produced inconsistent phenotype-plus-PRS tables.", call. = FALSE)
}
if (!file.copy(phenotypePATH[[1L]], "phenoPRS.csv", overwrite = TRUE)) {
  stop("Could not copy the combined phenotype-plus-PRS table.", call. = FALSE)
}
