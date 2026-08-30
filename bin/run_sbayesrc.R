#!/usr/bin/env Rscript

# Parse the selected official SBayesRC stage.
argument <- commandArgs(trailingOnly = TRUE)
name <- sub("^--", "", argument[seq.int(1L, length(argument), 2L)])
value <- argument[seq.int(2L, length(argument), 2L)]
option <- stats::setNames(as.list(value), name)
if (!requireNamespace("SBayesRC", quietly = TRUE)) stop("SBayesRC is required.", call. = FALSE)
if (!requireNamespace("data.table", quietly = TRUE)) stop("The data.table package is required.", call. = FALSE)

action <- option[["action"]]
traitID <- option[["trait-id"]]

# Copy the package log to one predictable stage-specific name.
copyLOG <- function(candidate, destination) {
  source <- candidate[file.exists(candidate)]
  if (length(source) == 0L) {
    stop(sprintf("SBayesRC did not create its expected log for %s.", action), call. = FALSE)
  }
  if (!file.copy(source[1L], destination, overwrite = TRUE)) {
    stop(sprintf("Could not save the SBayesRC log for %s.", action), call. = FALSE)
  }
}

if (action == "tidy") {
  source <- data.table::fread(option[["input"]], select = c("SNP", "freq"))
  if (any(!is.finite(source$freq) | source$freq <= 0 | source$freq >= 1)) {
    stop("SBayesRC tidy requires a valid effect-allele frequency for every GWAS variant.", call. = FALSE)
  }
  output <- paste0(traitID, ".tidy.ma")
  SBayesRC::tidy(mafile = option[["input"]], LDdir = option[["ld-dir"]], output = output, log2file = TRUE)
  copyLOG(paste0(output, ".log"), paste0(traitID, ".sbayesrc.tidy.log"))
  retained <- data.table::fread(output, select = "SNP")
  data.table::fwrite(
    data.table::data.table(
      trait_id = traitID,
      input_variants = nrow(source),
      ld_aligned_variants = nrow(retained),
      retained_percent = 100 * nrow(retained) / nrow(source),
      review_below_70_percent = nrow(retained) / nrow(source) < 0.70
    ),
    paste0(traitID, ".tidy_qc.tsv"),
    sep = "\t"
  )
} else if (action == "impute") {
  output <- paste0(traitID, ".imputed.ma")
  SBayesRC::impute(mafile = option[["input"]], LDdir = option[["ld-dir"]], output = output, log2file = TRUE)
  copyLOG(paste0(output, ".log"), paste0(traitID, ".sbayesrc.impute.log"))
  inputN <- nrow(data.table::fread(option[["input"]], select = "SNP"))
  outputN <- nrow(data.table::fread(output, select = "SNP"))
  data.table::fwrite(
    data.table::data.table(trait_id = traitID, ld_aligned_variants = inputN, summary_imputed_variants = outputN, status = "PASS"),
    paste0(traitID, ".impute_qc.tsv"),
    sep = "\t"
  )
} else if (action == "model") {
  output <- paste0(traitID, ".sbayesrc")
  set.seed(as.integer(option[["seed"]]))
  SBayesRC::sbayesrc(
    mafile = option[["input"]],
    LDdir = option[["ld-dir"]],
    outPrefix = output,
    annot = option[["annotation"]],
    log2file = TRUE
  )
  copyLOG(paste0(output, ".log"), paste0(traitID, ".sbayesrc.model.log"))
  weight <- data.table::fread(paste0(output, ".txt"), select = c("SNP", "BETA", "PIP"))
  if (nrow(weight) == 0L || any(!is.finite(weight$BETA)) || any(!is.finite(weight$PIP))) {
    stop("SBayesRC returned invalid posterior weights.", call. = FALSE)
  }
  data.table::fwrite(
    data.table::data.table(
      trait_id = traitID,
      weights = nrow(weight),
      nonzero_weights = sum(weight$BETA != 0),
      maximum_pip = max(weight$PIP),
      status = "PASS"
    ),
    paste0(traitID, ".sbayesrc.model_qc.tsv"),
    sep = "\t"
  )
} else if (action == "score") {
  targetPREFIX <- file.path(option[["target-dir"]], paste0(option[["cohort"]], "_chr{CHR}"))
  targetPSAM <- data.table::fread(
    file.path(option[["target-dir"]], paste0(option[["cohort"]], ".psam")),
    colClasses = "character"
  )
  iidCOLUMN <- if ("IID" %in% names(targetPSAM)) "IID" else "#IID"
  fidCOLUMN <- if ("#FID" %in% names(targetPSAM)) "#FID" else iidCOLUMN
  keep <- targetPSAM[, .(FID = get(fidCOLUMN), IID = get(iidCOLUMN))]
  keepPATH <- paste0(option[["cohort"]], ".keep.tsv")
  data.table::fwrite(keep, keepPATH, sep = "\t", col.names = FALSE)

  output <- paste0(option[["cohort"]], ".", traitID, ".sbayesrc")
  SBayesRC::prs(
    weight = option[["input"]],
    genoPrefix = targetPREFIX,
    outPrefix = output,
    genoCHR = "1-22",
    keepid = keepPATH,
    scoreFlag = "1 2 3 header no-mean-imputation",
    tool = option[["plink"]],
    log2file = TRUE
  )
  copyLOG(
    c(paste0(output, ".log.log"), paste0(output, ".log")),
    paste0(option[["cohort"]], ".", traitID, ".sbayesrc.score.log")
  )
  score <- data.table::fread(paste0(output, ".score.txt"), colClasses = "character")
  rawPRS <- suppressWarnings(as.numeric(score$SCORE))
  alleleCOUNT <- suppressWarnings(as.numeric(score$ALLELE_CT))
  if (nrow(score) != nrow(keep) || any(!is.finite(rawPRS)) || length(unique(alleleCOUNT)) != 1L) {
    stop("SBayesRC scoring did not return one finite score per prepared participant.", call. = FALSE)
  }
  usedN <- unique(alleleCOUNT) / 2
  weightN <- nrow(data.table::fread(option[["input"]], select = 1L))
  result <- data.table::data.table(
    cohort = option[["cohort"]],
    role = option[["role"]],
    trait_id = traitID,
    prs_name = option[["prs-name"]],
    method = "sbayesrc",
    FID = score$FID,
    IID = score$IID,
    raw_prs = rawPRS,
    allele_count = alleleCOUNT,
    used_variants = usedN
  )
  data.table::fwrite(result, paste0(option[["cohort"]], ".", traitID, ".sbayesrc.score.tsv"), sep = "\t", na = "NA")
  data.table::fwrite(
    data.table::data.table(
      cohort = option[["cohort"]], role = option[["role"]], trait_id = traitID,
      prs_name = option[["prs-name"]], method = "sbayesrc", participants = nrow(result),
      model_weights = weightN, finite_scores = sum(is.finite(result$raw_prs)), status = "PASS"
    ),
    paste0(option[["cohort"]], ".", traitID, ".sbayesrc.score_qc.tsv"),
    sep = "\t"
  )
} else {
  stop(sprintf("Unknown SBayesRC action: %s", action), call. = FALSE)
}
