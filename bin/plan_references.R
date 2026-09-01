#!/usr/bin/env Rscript

argument <- commandArgs(trailingOnly = TRUE)
name <- sub("^--", "", argument[seq.int(1L, length(argument), 2L)])
value <- argument[seq.int(2L, length(argument), 2L)]
option <- stats::setNames(as.list(value), name)

readTSV <- function(path) utils::read.delim(
  path, sep = "\t", header = TRUE, quote = "", comment.char = "",
  check.names = FALSE, stringsAsFactors = FALSE
)
writeTSV <- function(value, path) utils::write.table(
  value, path, sep = "\t", quote = FALSE, row.names = FALSE, na = ""
)

catalogue <- readTSV(option[["catalogue"]])
provided <- readTSV(option[["provided"]])
method <- trimws(strsplit(option[["methods"]], ",", fixed = TRUE)[[1L]])
targetIMPUTATION <- tolower(option[["target-imputation"]]) == "true"
mode <- tolower(option[["mode"]])

required <- c(
  "dbsnp", "reference_fasta", "imputation_panel", "population_panel",
  "related_samples", "unbref3_jar"
)
if ("sbayesrc" %in% method) required <- c(required, "sbayesrc_ld_source", "annotation_source")
if (targetIMPUTATION) required <- c(required, "genetic_map", "imputation_panel", "beagle_jar")
required <- unique(required)

catalogueROLE <- unique(catalogue$reference_type)
missingCATALOGUE <- setdiff(required, catalogueROLE)
if (length(missingCATALOGUE) > 0L) {
  stop(sprintf("The pinned catalogue is missing required role(s): %s", paste(missingCATALOGUE, collapse = ", ")), call. = FALSE)
}

providedROLE <- if (mode == "download" || nrow(provided) == 0L) character() else unique(provided$reference_type)
downloadROLE <- setdiff(required, providedROLE)
planned <- catalogue[catalogue$reference_type %in% downloadROLE, , drop = FALSE]
if (nrow(planned) == 0L) {
  planned <- catalogue[1L, , drop = FALSE]
  for (column in names(planned)) planned[[column]] <- ""
  planned$asset_id <- "cache_complete"
  planned$reference_type <- "cache_complete"
}
writeTSV(planned, "reference_assets.tsv")

status <- data.frame(
  reference_type = required,
  source = ifelse(required %in% providedROLE, "provided", "pinned_download"),
  status = "PLANNED",
  stringsAsFactors = FALSE
)
writeTSV(status, "reference_plan.tsv")
