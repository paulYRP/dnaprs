#!/usr/bin/env Rscript

source("assets/report/report-inputs.R")

testREPORTINPUTS <- function() {
  inputROOT <- tempfile("report-inputs-")
  dir.create(inputROOT)
  on.exit(unlink(inputROOT, recursive = TRUE), add = TRUE)
  checks <- 0L

  manifestFOR <- function(name) {
    data.frame(publish_path = "phenotype", file_name = name, stringsAsFactors = FALSE)
  }
  expectERROR <- function(name, message) {
    error <- tryCatch({
      validateREPORTINPUTS(inputROOT, manifestFOR(name))
      NULL
    }, error = conditionMessage)
    stopifnot(
      is.character(error),
      grepl(message, error, fixed = TRUE),
      grepl(file.path(inputROOT, name), error, fixed = TRUE)
    )
    checks <<- checks + 1L
    cat("PASS:", name, "is rejected with a named error\n")
  }

  sourcePATH <- file.path(inputROOT, "scores with spaces.csv")
  writeLines(c("UID,MDD_PRS_RAW,MDD_PRS", "TEST01,0.5,1.0"), sourcePATH)
  stopifnot(identical(validateREPORTINPUTS(inputROOT, manifestFOR(basename(sourcePATH))), sourcePATH))
  checks <- checks + 1L
  cat("PASS: regular files with spaces are accepted\n")

  linkPATH <- file.path(inputROOT, "scores-link.csv")
  stopifnot(file.symlink(sourcePATH, linkPATH))
  stopifnot(identical(validateREPORTINPUTS(inputROOT, manifestFOR(basename(linkPATH))), linkPATH))
  checks <- checks + 1L
  cat("PASS: symbolic links to readable regular files are accepted\n")

  directoryPATH <- file.path(inputROOT, "TEST.sbayesrc_genotypes")
  dir.create(directoryPATH)
  expectERROR(basename(directoryPATH), "is a directory")
  stopifnot(file.symlink(directoryPATH, file.path(inputROOT, "genotype-link")))
  expectERROR("genotype-link", "is a directory")
  expectERROR("missing.csv", "is missing or is not a regular file")
  stopifnot(file.symlink(file.path(inputROOT, "missing.csv"), file.path(inputROOT, "broken.csv")))
  expectERROR("broken.csv", "is missing or is not a regular file")

  unreadablePATH <- file.path(inputROOT, "unreadable.csv")
  writeLines("UID,MDD_PRS", unreadablePATH)
  Sys.chmod(unreadablePATH, "0000")
  on.exit(Sys.chmod(unreadablePATH, "0600"), add = TRUE, after = FALSE)
  if (file.access(unreadablePATH, mode = 4L) == 0L) {
    stop("Run report-input tests as a non-root user to verify unreadable files.")
  }
  expectERROR(basename(unreadablePATH), "is not readable")
  cat(sprintf("All %d report-input checks passed.\n", checks))
}

testREPORTINPUTS()
