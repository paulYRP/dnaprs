#!/usr/bin/env Rscript
# Generate synthetic report inputs without running genotype or PRS analysis.
suppressPackageStartupMessages(library(data.table))
arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) < 1L) stop("Usage: test_report_scale.R <new-directory> [rows-per-trait,...]")
destination <- arguments[1L]
sizes <- if (length(arguments) > 1L) as.integer(strsplit(arguments[2L], ",", fixed = TRUE)[[1L]]) else c(6500000L, 7000000L, 7500000L)
stopifnot(all(is.finite(sizes)), all(sizes > 0), !dir.exists(destination))
dir.create(destination, recursive = TRUE)
destination <- normalizePath(destination, winslash = "/")
inputs <- file.path(destination, "inputs")
dir.create(inputs)
manifest <- data.table(publish_path = character(), file_name = character())
writeINPUT <- function(name, value, section) {
  fwrite(value, file.path(inputs, name), sep = "\t")
  manifest <<- rbind(manifest, data.table(publish_path = section, file_name = name))
}
for (trait in seq_along(sizes)) local({
  n <- sizes[trait]
  name <- paste0("TRAIT", trait)
  message(sprintf("Generate %s: %d rows", name, n))
  index <- seq_len(n)
  value <- data.table(
    SNP = paste0("rs", index), CHR = (index - 1L) %% 22L + 1L,
    BP = (index - 1L) %/% 22L * 100L + 1L, A1 = "A", A2 = "C",
    freq = .01 + (index %% 4900L) / 10000, b = ((index %% 2001L) - 1000) / 10000,
    se = .02, p = ifelse(index %% 10L == 0L, 1e-12, .01 + (index %% 9900L) / 10000), N = 100000L
  )
  writeINPUT(paste0(name, ".cojo.ma"), value, paste0("gwas/", name))
  writeINPUT(paste0(name, ".sbayesrc.txt"), value[, .(SNP, BETA = b / 10, PIP = freq)], "sbayesrc")
  writeINPUT(paste0(name, ".harmonisation_qc.tsv"), data.table(trait_id = name, source_variants = n, harmonised_variants = n), paste0("qc/gwas/", name))
})
writeINPUT("variant_flow.tsv", data.table(cohort = "SYNTHETIC", role = "target", trait_id = "TRAIT1",
  prs_name = "PRS", method = "plink_ct", stage = c("Source GWAS", "Scored"),
  variant_count = sizes[1L], stage_order = 1:2, percent_of_source = 100), "qc")
writeINPUT("prs_scores_long.tsv", data.table(method = character(), cohort = character(), IID = character(), trait_id = character()), "scores")
writeINPUT("score_qc.tsv", data.table(cohort = character(), prs_name = character(), method = character(),
  participants = integer(), used_variants = integer(), status = character()), "qc")
writeINPUT("input_checks.tsv", data.table(check = "synthetic_report", status = "PASS"), "run")
writeINPUT("gwas.tsv", data.table(trait_id = paste0("TRAIT", seq_along(sizes)), prs_name = paste0("PRS", seq_along(sizes))), "run")
writeINPUT("targets.tsv", data.table(cohort = "SYNTHETIC", role = "target"), "run")
writeINPUT("SYNTHETIC.participant_decisions.tsv", data.table(
  cohort = "SYNTHETIC", IID = paste0("S", 1:6),
  primary_analysis = c(TRUE, FALSE, FALSE, FALSE, NA, TRUE),
  score_eligible = c(TRUE, TRUE, TRUE, FALSE, TRUE, FALSE),
  reason = c("Pass", "Related", "Ancestry", "QC", "Missing", "Conflicting")
), "target_qc")
csv <- file.path(inputs, "phenoPRS.csv")
stopifnot(file.copy("tests/data/phenotype_combine/one.phenoPRS.csv", csv))
manifest <- rbind(manifest, data.table(publish_path = "phenotype", file_name = "phenoPRS.csv"))
fwrite(manifest, file.path(destination, "output_files.tsv"), sep = "\t")
stopifnot(all(file.copy(list.files("assets/report", full.names = TRUE), destination, recursive = TRUE)))
message(sprintf("Synthetic report ready: %s", destination))
