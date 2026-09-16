#!/usr/bin/env Rscript
# Check task-local reference access without downloading production resources.
assembler <- parse("bin/assemble_references.R")
testROOT <- tempfile("reference-assembly-")
dir.create(testROOT)
testROOT <- normalizePath(testROOT, winslash = "/", mustWork = TRUE)

writeTSV <- function(value, path) utils::write.table(
  value, path, sep = "\t", quote = FALSE, row.names = FALSE, na = ""
)
encodeJSON <- function(value) jsonlite::base64_enc(charToRaw(jsonlite::toJSON(value, auto_unbox = TRUE)))
record <- function(role, path, companion = "", format = "directory") data.frame(
  bundle_id = "dnaprs_GRCh37", bundle_version = "test", reference_id = toupper(role),
  reference_type = role, path = path, companion = companion, build = "GRCh37",
  ancestry = "European", version = "test", checksum = "", source_format = format,
  reference_stage = "source", stringsAsFactors = FALSE
)
asset <- function(id, role, relative, reference, format = "directory") data.frame(
  asset_id = id, reference_type = role, relative_path = relative,
  reference_path = reference, companion = "", source_format = format,
  ancestry = "European", version = "test", url = "https://example.invalid/fixture",
  checksum_algorithm = "sha256", checksum = "fixture", size = "8", stringsAsFactors = FALSE
)
provided <- rbind(
  record("dbsnp", "rData/dbSNP157/source"),
  record("imputation_panel", "rData/imputation/reference", format = "bref3_directory"),
  record("reference_fasta", "rData/imputation/source/reference.fasta", "rData/imputation/source/reference.fasta.fai", "fasta"),
  record("genetic_map", "rData/imputation/map"),
  record("population_panel", "rData/population/source/metadata.tsv", format = "tsv"),
  record("related_samples", "rData/related/source/metadata.tsv", format = "tsv"),
  record("annotation_source", "rData/sbayesrc/annotation source.zip", format = "zip"),
  record("beagle_jar", "rData/imputation/source/beagle.jar", format = "jar"),
  record("unbref3_jar", "rData/imputation/source/unbref3.jar", format = "jar")
)
downloaded <- asset("sbayesrc_ld", "sbayesrc_ld_source", "sbayesrc/ld.zip", "sbayesrc/ld.zip", "zip")

runCASE <- function(name, mode = "relative", failure = "", download = TRUE) {
  task <- file.path(testROOT, name)
  dir.create(task)
  previous <- getwd()
  on.exit(setwd(previous))
  setwd(task)
  # Only staged copies exist. The declared source root is deliberately unavailable.
  base <- "/unmounted/reference-test/PRS"
  rows <- provided
  source <- unique(c(rows$path, rows$companion[nzchar(rows$companion)]))
  staged <- list()
  for (index in seq_along(source)) {
    destination <- file.path("provided_reference", sprintf("asset%02d", index), basename(source[[index]]))
    dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
    if (source[[index]] %in% rows$path[rows$source_format %in% c("directory", "bref3_directory")]) {
      dir.create(destination)
    } else {
      writeLines(source[[index]], destination)
    }
    staged[[paste0(base, "/", source[[index]])]] <- destination
  }
  dbsnp <- staged[[paste0(base, "/", rows$path[[1L]])]]
  panel <- staged[[paste0(base, "/", rows$path[[2L]])]]
  for (filename in c("assembly_report.txt", "GCF_000001405.25.gz", "GCF_000001405.25.gz.tbi")) {
    writeLines("fixture", file.path(dbsnp, filename))
  }
  for (chromosome in 1:22) writeLines("fixture", file.path(panel, sprintf("chr%s.bref3", chromosome)))
  if (mode == "absolute") {
    rows$path <- paste0(base, "/", rows$path)
    rows$companion[nzchar(rows$companion)] <- paste0(base, "/", rows$companion[nzchar(rows$companion)])
  } else {
    rows$path <- paste0("./unused/../", rows$path)
  }
  if (failure == "dbsnp") unlink(file.path(dbsnp, "assembly_report.txt"))
  if (failure == "panel") unlink(file.path(panel, "chr22.bref3"))
  if (failure == "companion") unlink(staged[[paste0(base, "/rData/imputation/source/reference.fasta.fai")]])
  writeTSV(rows, "provided.tsv")
  writeLines("fixture", "sbayesrc_ld")
  assets <- if (download) downloaded else transform(downloaded, asset_id = "cache_complete")
  arguments <- c(
    "--provided", "provided.tsv", "--provided-map-json", encodeJSON(staged),
    "--reference-base", base, "--assets-json", encodeJSON(assets), "--asset-dir", ".",
    "--bundle", "test", "--genome-build", "GRCh37", "--cache-root", "/cache/test"
  )
  environment <- new.env(parent = globalenv())
  environment$commandArgs <- function(trailingOnly = FALSE) arguments
  message <- tryCatch({ eval(assembler, environment); "" }, error = conditionMessage)
  if (nzchar(failure)) {
    expected <- switch(failure,
      dbsnp = "is missing required file(s): assembly_report.txt",
      panel = "must contain 22 BREF3 files; found 21",
      companion = "Reference 'REFERENCE_FASTA' companion cannot be read"
    )
    stopifnot(grepl(expected, message, fixed = TRUE), !file.exists("references.tsv"))
  } else {
    stopifnot(identical(message, ""))
    result <- read.delim("references.tsv", check.names = FALSE)
    local <- result[result$reference_type != "sbayesrc_ld_source", ]
    expected <- provided[match(local$reference_type, provided$reference_type), ]
    stopifnot(identical(local$path, paste0(base, "/", expected$path)))
    stopifnot(local$companion[local$reference_type == "reference_fasta"] == paste0(base, "/", expected$companion[expected$reference_type == "reference_fasta"]))
    stopifnot(nrow(read.delim("reference_receipt.tsv")) == as.integer(download))
    if (download) {
      added <- result[result$reference_type == "sbayesrc_ld_source", ]
      stopifnot(nrow(added) == 1L, file.exists(added$path), nchar(added$checksum) == 64L)
    }
  }
}

runCASE("mixed-relative")
runCASE("mixed-absolute", mode = "absolute")
runCASE("all-local", download = FALSE)
runCASE("missing-dbsnp", failure = "dbsnp")
runCASE("missing-panel", failure = "panel")
runCASE("missing-companion", failure = "companion")

runDOWNLOADED <- function() {
  task <- file.path(testROOT, "all-downloaded")
  dir.create(task)
  previous <- getwd()
  on.exit(setwd(previous))
  setwd(task)
  assets <- do.call(rbind, c(
    lapply(c("assembly_report.txt", "GCF_000001405.25.gz", "GCF_000001405.25.gz.tbi"), function(name) {
      asset(name, "dbsnp", paste0("dbsnp/", name), "dbsnp")
    }),
    lapply(1:22, function(chromosome) {
      name <- sprintf("chr%s.bref3", chromosome)
      asset(name, "imputation_panel", paste0("panel/", name), "panel", "bref3_directory")
    })
  ))
  for (name in assets$asset_id) writeLines("fixture", name)
  writeTSV(provided[0L, ], "provided.tsv")
  arguments <- c(
    "--provided", "provided.tsv", "--provided-map-json", encodeJSON(list()),
    "--reference-base", "/unmounted/reference-test/PRS", "--assets-json", encodeJSON(assets),
    "--asset-dir", ".", "--bundle", "test", "--genome-build", "GRCh37", "--cache-root", "/cache/test"
  )
  environment <- new.env(parent = globalenv())
  environment$commandArgs <- function(trailingOnly = FALSE) arguments
  eval(assembler, environment)
  result <- read.delim("references.tsv", check.names = FALSE)
  stopifnot(nrow(result) == 2L, all(dir.exists(result$path)))
  stopifnot(nrow(read.delim("reference_receipt.tsv")) == 25L)
}
runDOWNLOADED()
message("Reference assembly passed: mixed, local and downloaded sources, copied inputs, companions and incomplete references.")
