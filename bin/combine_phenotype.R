#!/usr/bin/env Rscript

argument <- commandArgs(trailingOnly = TRUE)
name <- sub("^--", "", argument[seq.int(1L, length(argument), 2L)])
value <- argument[seq.int(2L, length(argument), 2L)]
option <- stats::setNames(as.list(value), name)
if (!requireNamespace("data.table", quietly = TRUE)) stop("The data.table package is required.", call. = FALSE)

inputDIR <- option[["input-dir"]]

combineTABLE <- function(suffix, output, orderCOLUMN, uniqueROWS = FALSE) {
  paths <- sort(list.files(inputDIR, pattern = paste0("\\.", suffix, "$"), full.names = TRUE))
  if (length(paths) == 0L) stop(sprintf("No phenotype task output matched *.%s", suffix), call. = FALSE)
  tables <- lapply(paths, data.table::fread)
  combined <- data.table::rbindlist(tables, use.names = TRUE, fill = TRUE)
  if (uniqueROWS && nrow(combined) > 0L) combined <- unique(combined)
  ordering <- intersect(orderCOLUMN, names(combined))
  if (nrow(combined) > 0L && length(ordering) > 0L) data.table::setorderv(combined, ordering)
  data.table::fwrite(combined, output, sep = "\t", na = "NA")
  combined
}

association <- combineTABLE(
  "phenotype_associations.tsv", "phenotype_associations.tsv",
  c("model_id", "cohort", "trait_id", "method")
)
if ("permutation_holm" %in% names(association)) association[, permutation_holm := NULL]
association[, permutation_holm := NA_real_]
if (nrow(association) > 0L) {
  association[, permutation_holm := {
    adjusted <- rep(NA_real_, .N)
    selected <- primary %in% TRUE & is.finite(permutation_p)
    adjusted[selected] <- stats::p.adjust(permutation_p[selected], method = "holm")
    adjusted
  }, by = .(cohort, method)]
  data.table::setorder(association, model_id, cohort, trait_id, method)
  data.table::fwrite(association, "phenotype_associations.tsv", sep = "\t", na = "NA")
}
invisible(combineTABLE(
  "phenotype_models_fitted.tsv", "phenotype_models_fitted.tsv",
  c("model_id", "cohort", "method")
))
invisible(combineTABLE(
  "phenotype_plot_data.tsv", "phenotype_plot_data.tsv",
  c("model_id", "cohort", "method", "IID")
))
invisible(combineTABLE(
  "phenotype_permutations.tsv", "phenotype_permutations.tsv",
  c("model_id", "cohort", "method", "permutation_id")
))
invisible(combineTABLE(
  "phenotype_influence.tsv", "phenotype_influence.tsv",
  c("model_id", "cohort", "method", "IID")
))
invisible(combineTABLE(
  "phenotype_participant_level.tsv", "phenotype_participant_level.tsv",
  c("selection_model_id", "source_row"), uniqueROWS = TRUE
))
invisible(combineTABLE(
  "phenotype_timepoint_completeness.tsv", "phenotype_timepoint_completeness.tsv",
  c("model_id", "participant_id", "requested_timepoint"), uniqueROWS = TRUE
))

copyCONSISTENT <- function(pattern, output, label) {
  paths <- sort(list.files(inputDIR, pattern = pattern, full.names = TRUE))
  if (length(paths) == 0L) stop(sprintf("No phenotype task produced %s.", label), call. = FALSE)
  digest <- unname(tools::md5sum(paths))
  if (length(unique(digest)) != 1L) {
    stop(sprintf("Parallel phenotype tasks produced inconsistent %s tables.", label), call. = FALSE)
  }
  if (!file.copy(paths[[1L]], output, overwrite = TRUE)) {
    stop(sprintf("Could not copy the combined %s table.", label), call. = FALSE)
  }
}

copyCONSISTENT("\\.phenoPRS.csv$", "phenoPRS.csv", "legacy phenotype-plus-PRS")
copyCONSISTENT("\\.phenotype_with_prs.tsv$", "phenotype_with_prs.tsv", "row-level phenotype-plus-PRS")
