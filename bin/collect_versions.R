#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(yaml))

args <- commandArgs(trailingOnly = TRUE)
valueAfter <- function(flag) {
  position <- match(flag, args)
  if (is.na(position) || position == length(args)) stop("Missing argument: ", flag, call. = FALSE)
  args[[position + 1L]]
}

inputDIR <- valueAfter("--input-dir")
outputFILE <- valueAfter("--output")
versionFILES <- sort(list.files(inputDIR, pattern = "[.]ya?ml$", full.names = TRUE))
if (!length(versionFILES)) stop("No software-version records were supplied.", call. = FALSE)

records <- list()
for (versionFILE in versionFILES) {
  parsed <- yaml::read_yaml(versionFILE, eval.expr = FALSE)
  if (!is.list(parsed) || !length(parsed) || is.null(names(parsed)) || any(!nzchar(names(parsed)))) {
    stop("Invalid software-version YAML: ", basename(versionFILE), call. = FALSE)
  }
  duplicateKEYS <- intersect(names(records), names(parsed))
  if (length(duplicateKEYS)) {
    stop("Duplicate software-version key: ", duplicateKEYS[[1L]], call. = FALSE)
  }
  records[names(parsed)] <- parsed
}

records <- records[sort(names(records), method = "radix")]
yaml::write_yaml(records, outputFILE)
