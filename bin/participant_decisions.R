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
checks <- data.table::fread(option[["qc-checks"]], colClasses = "character")
required <- c("cohort", "FID", "IID", "missingness", "retained_after_qc", "decision", "reason")
if (!all(required %in% names(sample))) stop("The sample-decision table is invalid.", call. = FALSE)
if (anyDuplicated(sample[, .(cohort, FID, IID)])) stop("Sample decisions contain duplicate participants.", call. = FALSE)
cohortVALUE <- unique(sample$cohort)
if (length(cohortVALUE) != 1L) stop("Each participant-decision task requires one cohort.", call. = FALSE)
if ("cohort" %in% names(related) && nrow(related) > 0L) related[, cohort := as.character(cohortVALUE[[1L]])]
if ("cohort" %in% names(heterozygosity) && nrow(heterozygosity) > 0L) heterozygosity[, cohort := as.character(cohortVALUE[[1L]])]
if ("cohort" %in% names(sex) && nrow(sex) > 0L) sex[, cohort := as.character(cohortVALUE[[1L]])]

sample[, missingness := suppressWarnings(as.numeric(missingness))]
sample[, retained_after_qc := as.logical(retained_after_qc)]
if (anyNA(sample$retained_after_qc)) stop("Sample retention must be TRUE or FALSE for every participant.", call. = FALSE)
retainedKEY <- paste(sample[retained_after_qc == TRUE, FID], sample[retained_after_qc == TRUE, IID], sep = "\r")
if (!length(retainedKEY)) stop("No participant remains after sample QC.", call. = FALSE)
if (!all(c("check", "status", "reason") %in% names(checks))) stop("The diagnostic-status table is invalid.", call. = FALSE)
for (checkNAME in c("heterozygosity", "relatedness")) {
  checkROW <- checks[check == checkNAME]
  if (nrow(checkROW) != 1L || !checkROW$status %in% c("PASS", "REVIEW")) {
    stop(sprintf("Required %s QC did not complete successfully. Inspect %s and its genotype_eda log. %s",
      checkNAME, option[["qc-checks"]], paste(checkROW$reason, collapse = "; ")), call. = FALSE)
  }
}
validatePARTICIPANTS <- function(value, label) {
  if (!all(c("FID", "IID") %in% names(value))) stop(sprintf("%s results require FID and IID.", label), call. = FALSE)
  if (!"cohort" %in% names(value) || anyNA(value$cohort) || any(value$cohort != cohortVALUE)) {
    stop(sprintf("%s results must belong to cohort %s.", label, cohortVALUE), call. = FALSE)
  }
  key <- paste(value$FID, value$IID, sep = "\r")
  if (anyDuplicated(key) || !setequal(key, retainedKEY)) {
    stop(sprintf("%s results must contain each retained participant exactly once, with no extra participants.", label), call. = FALSE)
  }
}
sample[, sample_missingness_pass := decision %in% c("RETAIN", "INHERITED")]
if (!all(c("cohort", "FID", "IID", "heterozygosity_z", "status") %in% names(heterozygosity))) {
  stop("The heterozygosity table is invalid.", call. = FALSE)
}
validatePARTICIPANTS(heterozygosity, "Heterozygosity")
heterozygosity[, heterozygosity_z := suppressWarnings(as.numeric(heterozygosity_z))]
if (any(!is.finite(heterozygosity$heterozygosity_z)) || any(!heterozygosity$status %in% c("PASS", "REVIEW"))) {
  stop("Required heterozygosity QC contains invalid values or unsuccessful results. Inspect the genotype_eda log.", call. = FALSE)
}
heterozygosity[, heterozygosity_pass := status == "PASS"]
sample <- merge(
  sample,
  heterozygosity[, .(cohort, FID, IID, heterozygosity_z, heterozygosity_pass)],
  by = c("cohort", "FID", "IID"), all.x = TRUE, sort = FALSE
)
sample[is.na(heterozygosity_pass), heterozygosity_pass := FALSE]

if (!all(c("cohort", "FID", "IID", "status") %in% names(sex))) stop("The sex-check table is invalid.", call. = FALSE)
validatePARTICIPANTS(sex, "Sex check")
if (any(!sex$status %in% c("PASS", "REVIEW", "NOT_APPLICABLE"))) stop("Required sex QC contains unsuccessful results.", call. = FALSE)
sex[, `:=`(sex_check_pass = status %in% c("PASS", "NOT_APPLICABLE"), sex_check_status = status)]
sample <- merge(
  sample,
  sex[, .(cohort, FID, IID, sex_check_pass, sex_check_status)],
  by = c("cohort", "FID", "IID"), all.x = TRUE, sort = FALSE
)
sample[is.na(sex_check_pass), sex_check_pass := FALSE]
sample[is.na(sex_check_status), sex_check_status := "NOT_ASSESSED"]
sample[, technical_pass := retained_after_qc & sample_missingness_pass & heterozygosity_pass & sex_check_pass]
sample[, related_flag := FALSE]
if (!all(c("cohort", "FID", "IID", "ancestry_flag", "ancestry_distance") %in% names(ancestry))) {
  stop("The target ancestry table is invalid.", call. = FALSE)
}
validatePARTICIPANTS(ancestry, "Ancestry")
ancestry[, ancestry_distance := suppressWarnings(as.numeric(ancestry_distance))]
if (any(!is.finite(ancestry$ancestry_distance)) || any(!ancestry$ancestry_flag %in% c("PASS", "OUTLIER"))) {
  stop("Ancestry contains invalid distances or flags for retained participants.", call. = FALSE)
}
sample <- merge(
  sample,
  ancestry[, .(cohort, FID, IID, ancestry_flag, ancestry_distance)],
  by = c("cohort", "FID", "IID"), all.x = TRUE, sort = FALSE
)
sample[retained_after_qc == FALSE, ancestry_flag := "NOT_ASSESSED"]

# Both members of each flagged pair are excluded from the primary unrelated set.
relatedREQUIRED <- c("cohort", "FID1", "IID1", "FID2", "IID2", "pi_hat")
if (!all(relatedREQUIRED %in% names(related))) stop("The relatedness table is invalid.", call. = FALSE)
related[, pi_hat_numeric := suppressWarnings(as.numeric(pi_hat))]
firstKEY <- paste(related$FID1, related$IID1, sep = "\r")
secondKEY <- paste(related$FID2, related$IID2, sep = "\r")
pairKEY <- paste(pmin(firstKEY, secondKEY), pmax(firstKEY, secondKEY), sep = "\n")
if (nrow(related) != length(retainedKEY) * (length(retainedKEY) - 1) / 2 ||
    anyDuplicated(pairKEY) || any(firstKEY == secondKEY) ||
    any(!firstKEY %in% retainedKEY | !secondKEY %in% retainedKEY) || any(!is.finite(related$pi_hat_numeric))) {
  stop("Required relatedness QC must contain one finite PI_HAT per retained participant pair. Inspect the genotype_eda log.", call. = FALSE)
}
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

sample[, score_eligible := technical_pass]
sample[, primary_analysis := technical_pass & !related_flag & ancestry_flag == "PASS"]
sample[retained_after_qc == TRUE, reason := paste0(
  reason,
  ifelse(
    heterozygosity_pass,
    "; heterozygosity passed",
    "; heterozygosity failed"
  ),
  paste0("; sex check ", tolower(sex_check_status)),
  ifelse(related_flag, "; excluded from the primary unrelated-participant set", "; no relatedness exclusion"),
  ifelse(ancestry_flag == "PASS", "; within the European reference distance", "; outside the European reference distance")
)]
sample[retained_after_qc == FALSE, reason := paste0(reason, "; downstream diagnostics not assessed after sample exclusion")]
result <- sample[, .(
  cohort, FID, IID, missingness, retained_after_qc, sample_missingness_pass, heterozygosity_z,
  heterozygosity_pass, sex_check_pass, sex_check_status, technical_pass, score_eligible,
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
