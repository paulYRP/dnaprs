#!/usr/bin/env Rscript

# Convert a PLINK .sscore file to the common dnaprs score contract.
argument <- commandArgs(trailingOnly = TRUE)
name <- sub("^--", "", argument[seq.int(1L, length(argument), 2L)])
value <- argument[seq.int(2L, length(argument), 2L)]
option <- stats::setNames(as.list(value), name)
if (!requireNamespace("data.table", quietly = TRUE)) stop("The data.table package is required.", call. = FALSE)

score <- data.table::fread(option[["score"]], colClasses = "character")
scoreCOLUMN <- grep("_SUM$", names(score), value = TRUE)
scoreCOLUMN <- setdiff(scoreCOLUMN, "NAMED_ALLELE_DOSAGE_SUM")
if (length(scoreCOLUMN) != 1L) stop("PLINK output must contain one score-sum column.", call. = FALSE)
alleleCOLUMN <- intersect(c("ALLELE_CT", "NAMED_ALLELE_DOSAGE_SUM"), names(score))
if (length(alleleCOLUMN) == 0L) alleleCOUNT <- rep(NA_real_, nrow(score)) else alleleCOUNT <- as.numeric(score[[alleleCOLUMN[1L]]])
rawPRS <- suppressWarnings(as.numeric(score[[scoreCOLUMN]]))
if (any(!is.finite(rawPRS))) stop("PLINK returned a non-finite score.", call. = FALSE)

usedVARIANT <- data.table::fread(option[["used"]], header = FALSE)[[1L]]
if (length(usedVARIANT) == 0L || anyDuplicated(usedVARIANT)) stop("PLINK did not report a unique used-variant set.", call. = FALSE)

result <- data.table::data.table(
  cohort = option[["cohort"]],
  role = option[["role"]],
  trait_id = option[["trait-id"]],
  prs_name = option[["prs-name"]],
  method = "plink_ct",
  FID = as.character(score[["#FID"]]),
  IID = as.character(score[["IID"]]),
  raw_prs = rawPRS,
  allele_count = alleleCOUNT,
  used_variants = length(usedVARIANT)
)
data.table::fwrite(result, paste0(option[["cohort"]], ".", option[["trait-id"]], ".plink_ct.score.tsv"), sep = "\t")
data.table::fwrite(
  unique(result[, .(cohort, role, trait_id, prs_name, method, participants = .N, used_variants)]),
  paste0(option[["cohort"]], ".", option[["trait-id"]], ".plink_ct.score_qc.tsv"),
  sep = "\t"
)
