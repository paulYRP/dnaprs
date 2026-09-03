#!/usr/bin/env Rscript

argument <- commandArgs(trailingOnly = TRUE)
name <- sub("^--", "", argument[seq.int(1L, length(argument), 2L)])
value <- argument[seq.int(2L, length(argument), 2L)]
option <- stats::setNames(as.list(value), name)
if (!requireNamespace("data.table", quietly = TRUE)) stop("The data.table package is required.", call. = FALSE)

sample <- data.table::fread(option[["sample-decisions"]], colClasses = list(character = c("FID", "IID")))
related <- data.table::fread(option[["relatedness"]], colClasses = "character")
heterozygosity <- data.table::fread(option[["heterozygosity"]], colClasses = list(character = c("FID", "IID")))
sex <- data.table::fread(option[["sex-check"]], colClasses = list(character = c("FID", "IID")))
ancestry <- data.table::fread(option[["ancestry"]], colClasses = list(character = c("FID", "IID")))
required <- c("cohort", "FID", "IID", "missingness", "decision", "reason")
if (!all(required %in% names(sample))) stop("The sample-decision table is invalid.", call. = FALSE)
if (anyDuplicated(sample[, .(cohort, FID, IID)])) stop("Sample decisions contain duplicate participants.", call. = FALSE)
cohortVALUE <- unique(sample$cohort)
if (length(cohortVALUE) != 1L) stop("Each participant-decision task requires one cohort.", call. = FALSE)
if ("cohort" %in% names(related) && nrow(related) > 0L) related[, cohort := as.character(cohortVALUE[[1L]])]
if ("cohort" %in% names(heterozygosity) && nrow(heterozygosity) > 0L) heterozygosity[, cohort := as.character(cohortVALUE[[1L]])]
if ("cohort" %in% names(sex) && nrow(sex) > 0L) sex[, cohort := as.character(cohortVALUE[[1L]])]

sample[, missingness := suppressWarnings(as.numeric(missingness))]
sample[, sample_missingness_pass := decision %in% c("RETAIN", "INHERITED")]
if (!all(c("cohort", "FID", "IID", "heterozygosity_z", "status") %in% names(heterozygosity))) {
  stop("The heterozygosity table is invalid.", call. = FALSE)
}
if (nrow(heterozygosity) > 0L) {
  if (anyDuplicated(heterozygosity[, .(cohort, FID, IID)])) stop("Heterozygosity results contain duplicate participants.", call. = FALSE)
  heterozygosity[, heterozygosity_z := suppressWarnings(as.numeric(heterozygosity_z))]
  heterozygosity[, heterozygosity_pass := status == "PASS"]
  sample <- merge(
    sample,
    heterozygosity[, .(cohort, FID, IID, heterozygosity_z, heterozygosity_pass)],
    by = c("cohort", "FID", "IID"), all.x = TRUE, sort = FALSE
  )
  if (any(is.na(sample$heterozygosity_pass) & sample$sample_missingness_pass)) {
    stop("At least one retained sample-QC participant has no heterozygosity result.", call. = FALSE)
  }
  sample[is.na(heterozygosity_pass), heterozygosity_pass := FALSE]
} else {
  sample[, `:=`(heterozygosity_z = NA_real_, heterozygosity_pass = TRUE)]
}

if (nrow(sex) > 0L) {
  if (!all(c("cohort", "FID", "IID", "status") %in% names(sex))) stop("The sex-check table is invalid.", call. = FALSE)
  if (anyDuplicated(sex[, .(cohort, FID, IID)])) stop("Sex-check results contain duplicate participants.", call. = FALSE)
  sex[, sex_check_pass := status == "PASS"]
  sample <- merge(
    sample,
    sex[, .(cohort, FID, IID, sex_check_pass)],
    by = c("cohort", "FID", "IID"), all.x = TRUE, sort = FALSE
  )
  sample[is.na(sex_check_pass), sex_check_pass := FALSE]
} else {
  sample[, sex_check_pass := TRUE]
}
sample[, technical_pass := sample_missingness_pass & heterozygosity_pass & sex_check_pass]
sample[, related_flag := FALSE]
if (!all(c("cohort", "FID", "IID", "ancestry_flag", "ancestry_distance") %in% names(ancestry))) {
  stop("The target ancestry table is invalid.", call. = FALSE)
}
if (anyDuplicated(ancestry[, .(cohort, FID, IID)])) stop("Target ancestry contains duplicate participants.", call. = FALSE)
ancestry[, ancestry_distance := suppressWarnings(as.numeric(ancestry_distance))]
sample <- merge(
  sample,
  ancestry[, .(cohort, FID, IID, ancestry_flag, ancestry_distance)],
  by = c("cohort", "FID", "IID"), all.x = TRUE, sort = FALSE
)
if (anyNA(sample$ancestry_flag)) stop("At least one technical-QC participant has no ancestry result.", call. = FALSE)

# Match the R Markdown analysis: every participant in a pair with PI_HAT >= 0.1875 is
# excluded from the primary unrelated-participant set.
if (nrow(related) > 0L) {
  relatedREQUIRED <- c("cohort", "FID1", "IID1", "FID2", "IID2", "pi_hat")
  if (!all(relatedREQUIRED %in% names(related))) stop("The relatedness table is invalid.", call. = FALSE)
  related[, pi_hat_numeric := suppressWarnings(as.numeric(pi_hat))]
  related <- related[is.finite(pi_hat_numeric) & pi_hat_numeric >= 0.1875]
  participant_key <- paste(sample$cohort, sample$FID, sample$IID, sep = "\r")
  for (row in seq_len(nrow(related))) {
    key1 <- paste(related$cohort[[row]], related$FID1[[row]], related$IID1[[row]], sep = "\r")
    key2 <- paste(related$cohort[[row]], related$FID2[[row]], related$IID2[[row]], sep = "\r")
    index1 <- match(key1, participant_key)
    index2 <- match(key2, participant_key)
    if (!is.na(index1)) sample$related_flag[[index1]] <- TRUE
    if (!is.na(index2)) sample$related_flag[[index2]] <- TRUE
  }
}

sample[, score_eligible := technical_pass]
sample[, primary_analysis := technical_pass & !related_flag & ancestry_flag == "PASS"]
sample[, reason := paste0(
  reason,
  ifelse(
    heterozygosity_pass,
    ifelse(is.na(heterozygosity_z), "; heterozygosity was not applicable", "; heterozygosity passed"),
    "; heterozygosity failed"
  ),
  ifelse(sex_check_pass, "; sex check passed or was not applicable", "; sex check failed"),
  ifelse(related_flag, "; excluded from the primary unrelated-participant set", "; no relatedness exclusion"),
  ifelse(ancestry_flag == "PASS", "; within the European reference distance", "; outside the European reference distance")
)]
result <- sample[, .(
  cohort, FID, IID, missingness, sample_missingness_pass, heterozygosity_z,
  heterozygosity_pass, sex_check_pass, technical_pass, score_eligible,
  related_flag, ancestry_flag, ancestry_distance, primary_analysis, reason
)]
data.table::setorder(result, cohort, FID, IID)
if (!any(result$score_eligible)) stop("No participant remains eligible for scoring.", call. = FALSE)
data.table::fwrite(result, option[["output"]], sep = "\t", na = "NA")
data.table::fwrite(
  result[score_eligible == TRUE, .(FID, IID)],
  option[["keep"]],
  sep = "\t", col.names = FALSE
)
