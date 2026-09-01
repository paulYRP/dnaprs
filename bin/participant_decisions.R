#!/usr/bin/env Rscript

argument <- commandArgs(trailingOnly = TRUE)
name <- sub("^--", "", argument[seq.int(1L, length(argument), 2L)])
value <- argument[seq.int(2L, length(argument), 2L)]
option <- stats::setNames(as.list(value), name)
if (!requireNamespace("data.table", quietly = TRUE)) stop("The data.table package is required.", call. = FALSE)

sample <- data.table::fread(option[["sample-decisions"]], colClasses = list(character = c("FID", "IID")))
related <- data.table::fread(option[["relatedness"]], colClasses = "character")
ancestry <- data.table::fread(option[["ancestry"]], colClasses = list(character = c("FID", "IID")))
required <- c("cohort", "FID", "IID", "missingness", "decision", "reason")
if (!all(required %in% names(sample))) stop("The sample-decision table is invalid.", call. = FALSE)
if (anyDuplicated(sample[, .(cohort, FID, IID)])) stop("Sample decisions contain duplicate participants.", call. = FALSE)

sample[, missingness := suppressWarnings(as.numeric(missingness))]
sample[, technical_pass := decision %in% c("RETAIN", "INHERITED")]
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

# Select one participant from each reviewable related pair deterministically. Prefer the
# participant with lower missingness; use cohort/FID/IID order as the stable tie-break.
if (nrow(related) > 0L && all(c("cohort", "FID1", "IID1", "FID2", "IID2", "kinship", "review_status") %in% names(related))) {
  related[, kinship_numeric := suppressWarnings(as.numeric(kinship))]
  related <- related[review_status != "PASS"]
  data.table::setorder(related, -kinship_numeric, cohort, FID1, IID1, FID2, IID2, na.last = TRUE)
  participant_key <- paste(sample$cohort, sample$FID, sample$IID, sep = "\r")
  for (row in seq_len(nrow(related))) {
    key1 <- paste(related$cohort[[row]], related$FID1[[row]], related$IID1[[row]], sep = "\r")
    key2 <- paste(related$cohort[[row]], related$FID2[[row]], related$IID2[[row]], sep = "\r")
    index1 <- match(key1, participant_key)
    index2 <- match(key2, participant_key)
    if (is.na(index1) || is.na(index2) || !sample$technical_pass[[index1]] || !sample$technical_pass[[index2]]) next
    if (sample$related_flag[[index1]] || sample$related_flag[[index2]]) next
    miss1 <- sample$missingness[[index1]]
    miss2 <- sample$missingness[[index2]]
    if (!is.finite(miss1)) miss1 <- Inf
    if (!is.finite(miss2)) miss2 <- Inf
    drop_index <- if (miss1 > miss2) index1 else if (miss2 > miss1) index2 else {
      if (key1 > key2) index1 else index2
    }
    sample$related_flag[[drop_index]] <- TRUE
  }
}

sample[, score_eligible := technical_pass]
sample[, primary_analysis := technical_pass & !related_flag & ancestry_flag == "PASS"]
sample[, reason := paste0(
  reason,
  ifelse(related_flag, "; excluded from the primary unrelated-participant set", "; no relatedness exclusion"),
  ifelse(ancestry_flag == "PASS", "; within the European reference distance", "; outside the European reference distance")
)]
result <- sample[, .(
  cohort, FID, IID, missingness, technical_pass, score_eligible,
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
