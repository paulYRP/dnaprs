#!/usr/bin/env Rscript

argument <- commandArgs(trailingOnly = TRUE)
if (length(argument) != 1L || !argument[[1L]] %in% c("analysis", "report")) {
  stop("Usage: check_r_dependencies.R <analysis|report>", call. = FALSE)
}

environment <- argument[[1L]]
manifestPATH <- switch(
  environment,
  analysis = "containers/analysis/R-packages.tsv",
  report = c("containers/analysis/R-packages.tsv", "containers/report/R-packages.tsv")
)
sourcePATH <- switch(
  environment,
  analysis = setdiff(list.files("bin", pattern = "[.]R$", full.names = TRUE), "bin/run_sbayesrc.R"),
  report = list.files("assets/report", pattern = "[.](R|qmd)$", full.names = TRUE, recursive = TRUE)
)

manifest <- do.call(
  rbind,
  lapply(manifestPATH, read.delim, check.names = FALSE, stringsAsFactors = FALSE)
)
manifest <- manifest[!duplicated(manifest$package, fromLast = TRUE), ]

installed <- vapply(manifest$package, requireNamespace, logical(1L), quietly = TRUE)
if (any(!installed)) {
  stop(
    sprintf(
      "%s container is missing declared R package(s): %s.",
      environment,
      paste(manifest$package[!installed], collapse = ", ")
    ),
    call. = FALSE
  )
}

observedVERSION <- vapply(
  manifest$package,
  function(package) as.character(packageVersion(package)),
  character(1L)
)
incorrect <- observedVERSION != manifest$expected_version
if (any(incorrect)) {
  detail <- sprintf(
    "%s=%s (expected %s)",
    manifest$package[incorrect],
    observedVERSION[incorrect],
    manifest$expected_version[incorrect]
  )
  stop(
    sprintf("%s container has incorrect R package version(s): %s.", environment, paste(detail, collapse = ", ")),
    call. = FALSE
  )
}

sourceTEXT <- paste(
  unlist(lapply(sourcePATH, readLines, warn = FALSE), use.names = FALSE),
  collapse = "\n"
)

extractMATCH <- function(pattern, value, transform = identity) {
  match <- regmatches(value, gregexpr(pattern, value, perl = TRUE))[[1L]]
  if (length(match) == 0L) character() else transform(match)
}

namespacePACKAGE <- extractMATCH(
  "[A-Za-z][A-Za-z0-9.]*::",
  sourceTEXT,
  function(value) sub("::$", "", value)
)
quotedPACKAGE <- extractMATCH(
  "(?:library|require|requireNamespace)\\s*\\(\\s*['\"][A-Za-z][A-Za-z0-9.]*['\"]",
  sourceTEXT,
  function(value) sub(".*['\"]([A-Za-z][A-Za-z0-9.]*)['\"]$", "\\1", value)
)
libraryPACKAGE <- extractMATCH(
  "(?:library|require)\\s*\\(\\s*[A-Za-z][A-Za-z0-9.]*",
  sourceTEXT,
  function(value) sub(".*\\(\\s*", "", value)
)

basePACKAGE <- c(
  "base", "compiler", "datasets", "graphics", "grDevices", "grid", "methods",
  "parallel", "splines", "stats", "stats4", "tcltk", "tools", "utils"
)
requiredPACKAGE <- sort(unique(c(namespacePACKAGE, quotedPACKAGE, libraryPACKAGE)))
missingDECLARATION <- setdiff(requiredPACKAGE, c(basePACKAGE, manifest$package))
if (length(missingDECLARATION) > 0L) {
  stop(
    sprintf(
      "%s source uses R package(s) absent from its manifests: %s.",
      environment,
      paste(missingDECLARATION, collapse = ", ")
    ),
    call. = FALSE
  )
}

message(
  sprintf(
    "Verified %d %s R package version(s) and source declarations.",
    nrow(manifest),
    environment
  )
)
