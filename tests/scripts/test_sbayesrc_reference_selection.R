#!/usr/bin/env Rscript
# Exercise reference discovery without downloading production resources.
resolver <- new.env(parent = globalenv())
for (expression in parse("bin/resolve_inputs.R")) {
  if (is.call(expression) && identical(expression[[1L]], as.name("<-")) &&
      is.call(expression[[3L]]) && identical(expression[[3L]][[1L]], as.name("function"))) {
    eval(expression, resolver)
  }
}
resolver$genomeBUILD <- "GRCh37"
resolver$option <- list("reference-bundle" = "test", "beagle-jar" = "", "unbref3-jar" = "")
directory <- tempfile("sbayesrc-discovery-")
dir.create(directory)
stopifnot(file.copy("tests/data/reference/sbayesrc_source.zip", file.path(directory, "ukbEUR_HM3.zip")))
required <- "sbayesrc_ld_source"
failure <- tryCatch({ resolver$discoverREFERENCES(directory, directory, "local", required); "" }, error = conditionMessage)
stopifnot(grepl("production uses ukbEUR_Imputed.zip", failure, fixed = TRUE))
automatic <- resolver$discoverREFERENCES(directory, directory, "auto", required)
stopifnot(nrow(automatic) == 0L)
stopifnot(file.copy("tests/data/reference/sbayesrc_source.zip", file.path(directory, "ukbEUR_Imputed.zip")))
selected <- resolver$discoverREFERENCES(directory, directory, "local", required)
stopifnot(basename(selected$path[selected$reference_type == required]) == "ukbEUR_Imputed.zip")
explicit <- resolver$explicitREFERENCES(list(root = directory, assets = list(sbayesrc_ld_source = "ukbEUR_HM3.zip")))
stopifnot(basename(explicit$path) == "ukbEUR_HM3.zip")
catalogue <- read.delim("assets/reference_catalogue.tsv", check.names = FALSE)
production <- catalogue[catalogue$reference_type == required, ]
stopifnot(nrow(production) == 1L, grepl("LD/Imputed/ukbEUR_Imputed.zip", production$url, fixed = TRUE))
message("SBayesRC reference selection passed: full European default, automatic replacement and explicit source override.")
