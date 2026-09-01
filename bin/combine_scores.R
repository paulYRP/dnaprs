#!/usr/bin/env Rscript

# Combine method-specific scores and standardise only within each target cohort.
argument <- commandArgs(trailingOnly = TRUE)
name <- sub("^--", "", argument[seq.int(1L, length(argument), 2L)])
value <- argument[seq.int(2L, length(argument), 2L)]
option <- stats::setNames(as.list(value), name)
if (!requireNamespace("data.table", quietly = TRUE)) stop("The data.table package is required.", call. = FALSE)

scoreFILE <- strsplit(option[["scores"]], ",", fixed = TRUE)[[1L]]
score <- data.table::rbindlist(lapply(scoreFILE, data.table::fread), use.names = TRUE, fill = TRUE)
required <- c("cohort", "role", "trait_id", "prs_name", "method", "FID", "IID", "raw_prs", "used_variants")
if (!all(required %in% names(score))) stop("A method score does not follow the dnaprs score contract.", call. = FALSE)
if (any(!is.finite(score$raw_prs))) stop("At least one raw PRS is not finite.", call. = FALSE)
if (anyDuplicated(score[, .(cohort, trait_id, method, IID)])) stop("A participant has a duplicated method score.", call. = FALSE)

decisionFILE <- strsplit(option[["participant-decisions"]], ",", fixed = TRUE)[[1L]]
decision <- data.table::rbindlist(
  lapply(decisionFILE, data.table::fread), use.names = TRUE, fill = TRUE
)
decisionREQUIRED <- c(
  "cohort", "FID", "IID", "technical_pass", "score_eligible", "related_flag",
  "ancestry_flag", "ancestry_distance", "primary_analysis"
)
if (!all(decisionREQUIRED %in% names(decision))) stop("A participant-decision table is invalid.", call. = FALSE)
if (anyDuplicated(decision[, .(cohort, FID, IID)])) stop("Participant decisions contain duplicate IDs.", call. = FALSE)
score <- merge(
  score,
  decision[, ..decisionREQUIRED],
  by = c("cohort", "FID", "IID"), all.x = TRUE, sort = FALSE
)
if (any(is.na(score$score_eligible)) || any(!score$score_eligible)) {
  stop("A generated score has no matching score-eligible participant decision.", call. = FALSE)
}

score[, `:=`(
  cohort_mean = mean(raw_prs),
  cohort_sd = stats::sd(raw_prs)
), by = .(cohort, role, trait_id, prs_name, method)]
if (any(!is.finite(score$cohort_sd) | score$cohort_sd <= 0)) {
  stop("At least one cohort-method PRS has no participant variation.", call. = FALSE)
}
score[, prs_z := (raw_prs - cohort_mean) / cohort_sd]
score[, score_name := ifelse(
  method == "plink_ct",
  prs_name,
  paste0(prs_name, "_", toupper(method))
)]
data.table::setorder(score, cohort, trait_id, method, IID)
data.table::fwrite(score, "prs_scores_long.tsv", sep = "\t", na = "NA")

wideZ <- data.table::dcast(score, cohort + role + FID + IID ~ score_name, value.var = "prs_z")
score[, raw_name := paste0(score_name, "_RAW")]
wideRAW <- data.table::dcast(score, cohort + role + FID + IID ~ raw_name, value.var = "raw_prs")
wide <- merge(wideZ, wideRAW, by = c("cohort", "role", "FID", "IID"), all = TRUE, sort = FALSE)
wide <- merge(
  wide,
  unique(score[, .(
    cohort, FID, IID, technical_pass, score_eligible, related_flag,
    ancestry_flag, ancestry_distance, primary_analysis
  )]),
  by = c("cohort", "FID", "IID"), all.x = TRUE, sort = FALSE
)
data.table::fwrite(wide, "prs_scores_wide.tsv", sep = "\t", na = "NA")

scoreQC <- score[, .(
  participants = .N,
  finite_scores = sum(is.finite(raw_prs)),
  raw_mean = unique(cohort_mean),
  raw_sd = unique(cohort_sd),
  z_mean = mean(prs_z),
  z_sd = stats::sd(prs_z),
  used_variants = if (all(is.na(used_variants))) NA_integer_ else unique(used_variants)
), by = .(cohort, role, trait_id, prs_name, method)]
scoreQC[, status := ifelse(finite_scores == participants & abs(z_mean) < 1e-10 & abs(z_sd - 1) < 1e-10, "PASS", "FAIL")]
if (any(scoreQC$status != "PASS")) stop("The combined score checks failed.", call. = FALSE)
data.table::fwrite(scoreQC, "score_qc.tsv", sep = "\t", na = "NA")

methodCOUNT <- score[, .(methods = data.table::uniqueN(method)), by = .(cohort, trait_id, prs_name)]
if (any(methodCOUNT$methods > 1L)) {
  concordance <- methodCOUNT[methods > 1L, {
    value <- data.table::dcast(
      score[cohort == .BY$cohort & trait_id == .BY$trait_id & prs_name == .BY$prs_name],
      IID ~ method,
      value.var = "prs_z"
    )
    methodNAME <- setdiff(names(value), "IID")
    pair <- utils::combn(methodNAME, 2L)
    data.table::data.table(
      method_1 = pair[1L, ],
      method_2 = pair[2L, ],
      participants = nrow(value),
      pearson_r = vapply(seq_len(ncol(pair)), function(column) stats::cor(value[[pair[1L, column]]], value[[pair[2L, column]]]), numeric(1L)),
      spearman_r = vapply(seq_len(ncol(pair)), function(column) stats::cor(value[[pair[1L, column]]], value[[pair[2L, column]]], method = "spearman"), numeric(1L))
    )
  }, by = .(cohort, trait_id, prs_name)]
} else {
  concordance <- data.table::data.table(
    cohort = character(), trait_id = character(), prs_name = character(), method_1 = character(),
    method_2 = character(), participants = integer(), pearson_r = numeric(), spearman_r = numeric()
  )
}
data.table::fwrite(concordance, "method_concordance.tsv", sep = "\t", na = "NA")
