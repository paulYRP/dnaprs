#!/usr/bin/env Rscript
library(data.table)
root <- normalizePath(".")
fixture <- file.path(root, "tests/data/participant_decisions")
work <- tempfile("participant-inputs-")
dir.create(work)
setwd(work)
rscript <- file.path(R.home("bin"), "Rscript")
runSCRIPT <- function(script, args, error = NULL) {
  log <- suppressWarnings(system2(rscript, c(shQuote(file.path(root, "bin", script)), shQuote(args)), stdout = TRUE, stderr = TRUE))
  status <- attr(log, "status")
  if (is.null(error)) {
    if (!is.null(status) && status != 0L) stop(paste(log, collapse = "\n"))
  } else {
    stopifnot(!is.null(status), status != 0L, any(grepl(error, log, fixed = TRUE)))
  }
}
files <- c("sample_decisions", "relatedness", "heterozygosity", "sex_check", "target_ancestry", "genotype_eda_checks")
original <- setNames(lapply(files, function(name) fread(file.path(fixture, paste0("TEST.", name, ".tsv")))), files)
runQC <- function(tables = original, error = NULL) {
  for (name in names(tables)) fwrite(tables[[name]], paste0(name, ".tsv"), sep = "\t", na = "NA")
  runSCRIPT("participant_decisions.R", c(
    "--sample-decisions", "sample_decisions.tsv", "--relatedness", "relatedness.tsv",
    "--heterozygosity", "heterozygosity.tsv", "--sex-check", "sex_check.tsv",
    "--ancestry", "target_ancestry.tsv", "--qc-checks", "genotype_eda_checks.tsv",
    "--output", "decisions.tsv", "--keep", "eligible.keep"
  ), error)
  if (is.null(error)) fread("decisions.tsv")
}
baseline <- runQC()
stopifnot(nrow(baseline) == 4L, all(baseline$retained_after_qc), all(baseline$sex_check_status == "NOT_APPLICABLE"))
tables <- copy(original)
removed <- copy(tables$sample_decisions[1])
removed[, `:=`(FID = "REMOVED", IID = "REMOVED", retained_after_qc = FALSE,
  decision = "EXCLUDE", missingness = 0.2, reason = "Sample missingness exceeds the threshold")]
tables$sample_decisions <- rbind(tables$sample_decisions, removed)
result <- runQC(tables)
stopifnot(nrow(result) == 5L, !result[IID == "REMOVED", score_eligible],
  result[IID == "REMOVED", ancestry_flag] == "NOT_ASSESSED",
  !grepl("outside the European", result[IID == "REMOVED", reason], fixed = TRUE))
for (name in c("heterozygosity", "target_ancestry", "sex_check")) {
  tables <- copy(original)
  tables[[name]] <- tables[[name]][-1]
  runQC(tables, "each retained participant exactly once")
}
for (name in c("heterozygosity", "relatedness")) {
  tables <- copy(original)
  tables[[name]] <- tables[[name]][0]
  runQC(tables, if (name == "relatedness") "one finite PI_HAT" else "each retained participant exactly once")
  tables <- copy(original)
  tables$genotype_eda_checks[check == name, status := "FAIL"]
  runQC(tables, paste("Required", name, "QC did not complete successfully"))
}
tables <- copy(original)
tables$relatedness[1, pi_hat := NA_real_]
runQC(tables, "one finite PI_HAT")
tables <- copy(original)
tables$relatedness[, pi_hat := 0.01]
stopifnot(!any(runQC(tables)$related_flag))
cat("Retained participants require complete QC; excluded participants remain documented.\n")

phenotype <- fread(file.path(root, "tests/data/phenotype_repeated.tsv"))
scores <- fread(file.path(root, "tests/data/scores/test_timepoint.score.tsv"))
oldCURRENT <- seq_len(nrow(phenotype)) / 10
oldARCHIVE <- oldCURRENT + 1
oldARCHIVE[2] <- NA_real_
for (current in c(FALSE, TRUE)) for (archive in c(FALSE, TRUE)) {
  input <- copy(phenotype)
  input[, primary_analysis := "Original phenotype value"]
  if (current) input[, MDD_PRS := oldCURRENT]
  if (archive) input[, MDD_PRS_NIMP := oldARCHIVE]
  fwrite(input, "phenotype.tsv", sep = "\t", na = "NA")
  checksum <- tools::md5sum("phenotype.tsv")
  runSCRIPT("phenotype_association.R", c(
    "--scores", file.path(root, "tests/data/scores/test_timepoint.score.tsv"),
    "--phenotype", "phenotype.tsv", "--models", file.path(root, "tests/data/phenotype_models_timepoint.tsv"),
    "--model-id", "baseline", "--cohort", "TEST", "--trait-id", "MDD", "--method", "plink_ct", "--seed", "22"
  ))
  output <- fread("phenotype_with_prs.tsv")
  stopifnot(identical(checksum, tools::md5sum("phenotype.tsv")), nrow(output) == nrow(input),
    identical(output$IID, input$IID), identical(output$Visit, input$Visit),
    identical(output$Timepoint, input$Timepoint), identical(output$OUTCOME, input$OUTCOME),
    identical(output$Age, input$Age), identical(output$primary_analysis, input$primary_analysis),
    isTRUE(all.equal(output$MDD_PRS, scores$prs_z[match(input$IID, scores$IID)])))
  if (archive || current) stopifnot(isTRUE(all.equal(output$MDD_PRS_NIMP, if (archive) oldARCHIVE else oldCURRENT)))
  else stopifnot(!"MDD_PRS_NIMP" %in% names(output))
  stopifnot(nrow(fread("phenotype_participant_level.tsv")) == 6L)
}
cat("All four current-score/archive layouts preserve phenotype inputs and timepoints.\n")
