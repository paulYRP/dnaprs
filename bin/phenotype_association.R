#!/usr/bin/env Rscript

# Evaluate prespecified matching phenotype-PRS models without tuning score construction.
argument <- commandArgs(trailingOnly = TRUE)
name <- sub("^--", "", argument[seq.int(1L, length(argument), 2L)])
value <- argument[seq.int(2L, length(argument), 2L)]
option <- stats::setNames(as.list(value), name)
if (!requireNamespace("data.table", quietly = TRUE)) stop("The data.table package is required.", call. = FALSE)
for (item in c("cohort", "trait-id", "method")) if (is.null(option[[item]])) option[[item]] <- ""

score <- data.table::fread(option[["scores"]], colClasses = list(character = c("FID", "IID")))
phenotype <- data.table::fread(option[["phenotype"]])
if ("source_row" %in% names(phenotype)) {
  stop("Phenotype input cannot contain the reserved column 'source_row'.", call. = FALSE)
}
phenotype[, source_row := .I]
modelSPEC <- data.table::fread(option[["models"]])
if (is.null(option[["model-id"]]) || !nzchar(option[["model-id"]])) stop("--model-id is required.", call. = FALSE)
modelSPEC <- modelSPEC[model_id == option[["model-id"]]]
if (nrow(modelSPEC) != 1L) stop("--model-id must select exactly one resolved phenotype model.", call. = FALSE)
for (column in setdiff(c("control_value", "case_value"), names(modelSPEC))) {
  modelSPEC[[column]] <- ""
}
if (!"model_type" %in% names(modelSPEC)) {
  modelSPEC[, model_type := data.table::fifelse(!is.na(group_id) & group_id != "", "mixed", family)]
}
for (column in setdiff(c("timepoint_column", "timepoint_values"), names(modelSPEC))) {
  modelSPEC[[column]] <- ""
}

# Add the completed PRSs to the original phenotype rows without changing row order.
participantCOLUMN <- unique(modelSPEC$participant_id)
participantCOLUMN <- participantCOLUMN[!is.na(participantCOLUMN) & participantCOLUMN != ""]
if (length(participantCOLUMN) == 0L && "IID" %in% names(phenotype)) {
  participantCOLUMN <- "IID"
}
if (length(participantCOLUMN) > 1L) {
  stop("All phenotype models must use one participant identifier column.", call. = FALSE)
}
if (length(participantCOLUMN) == 1L) {
  participantCOLUMN <- participantCOLUMN[[1L]]
  phenotype[[participantCOLUMN]] <- as.character(phenotype[[participantCOLUMN]])
  scoreWIDEZ <- data.table::dcast(
    score,
    cohort + IID ~ score_name,
    value.var = "prs_z"
  )
  score[, raw_score_name := paste0(score_name, "_RAW")]
  scoreWIDERAW <- data.table::dcast(
    score,
    cohort + IID ~ raw_score_name,
    value.var = "raw_prs"
  )
  decisionCOLUMN <- intersect(
    c(
      "cohort", "FID", "IID", "technical_pass", "score_eligible", "related_flag",
      "sample_missingness_pass", "heterozygosity_z", "heterozygosity_pass", "sex_check_pass",
      "ancestry_flag", "ancestry_distance", "primary_analysis"
    ),
    names(score)
  )
  decision <- unique(score[, ..decisionCOLUMN])
  if (anyDuplicated(decision[, .(cohort, IID)])) {
    stop("Participant decisions are inconsistent across score records.", call. = FALSE)
  }
  if ("FID" %in% names(decision)) data.table::setnames(decision, "FID", "genotype_fid")
  scoreWIDE <- merge(scoreWIDEZ, scoreWIDERAW, by = c("cohort", "IID"), all = TRUE, sort = FALSE)
  scoreWIDE <- merge(scoreWIDE, decision, by = c("cohort", "IID"), all.x = TRUE, sort = FALSE)
  if (data.table::uniqueN(scoreWIDE$cohort) > 1L && !"cohort" %in% names(phenotype)) {
    stop("Phenotype data must contain a cohort column when scores contain multiple cohorts.", call. = FALSE)
  }
  scoreCOLUMN <- setdiff(names(scoreWIDE), c("cohort", "IID"))
  plinkCOLUMN <- unique(score[method == "plink_ct", score_name])
  for (column in intersect(plinkCOLUMN, names(phenotype))) {
    archivedCOLUMN <- paste0(column, "_NIMP")
    if (archivedCOLUMN %in% names(phenotype)) {
      stop(sprintf("Phenotype data already contain both '%s' and '%s'.", column, archivedCOLUMN), call. = FALSE)
    }
    data.table::setnames(phenotype, column, archivedCOLUMN)
  }
  if (length(scoreCOLUMN) > 0L) {
    collision <- intersect(scoreCOLUMN, names(phenotype))
    for (column in collision) {
      archivedCOLUMN <- paste0(column, "_INPUT")
      if (archivedCOLUMN %in% names(phenotype)) {
        stop(sprintf("Phenotype data already contain both '%s' and '%s'.", column, archivedCOLUMN), call. = FALSE)
      }
      data.table::setnames(phenotype, column, archivedCOLUMN)
    }
    if ("cohort" %in% names(phenotype)) {
      scoreKEY <- paste(scoreWIDE$cohort, scoreWIDE$IID, sep = "\r")
      phenotypeKEY <- paste(phenotype$cohort, phenotype[[participantCOLUMN]], sep = "\r")
    } else {
      if (anyDuplicated(scoreWIDE$IID)) {
        stop("Participant identifiers are duplicated across score cohorts.", call. = FALSE)
      }
      scoreKEY <- scoreWIDE$IID
      phenotypeKEY <- phenotype[[participantCOLUMN]]
    }
    scoreROW <- match(phenotypeKEY, scoreKEY)
    for (column in scoreCOLUMN) {
      phenotype[[column]] <- scoreWIDE[[column]][scoreROW]
    }
  }
}
data.table::fwrite(phenotype, "phenoPRS.csv", sep = ",", na = "NA")
data.table::fwrite(phenotype, "phenotype_with_prs.tsv", sep = "\t", na = "NA")

splitVALUES <- function(value, delimiter = ",") {
  value <- as.character(value)
  if (length(value) == 0L || is.na(value) || !nzchar(trimws(value))) return(character())
  result <- trimws(strsplit(value, delimiter, fixed = TRUE)[[1L]])
  result[nzchar(result)]
}

specification <- modelSPEC[1L]
groupedMODEL <- identical(tolower(as.character(specification$model_type)), "mixed") ||
  (!is.na(specification$group_id) && nzchar(as.character(specification$group_id)))
timepointCOLUMN <- as.character(specification$timepoint_column)
if (is.na(timepointCOLUMN)) timepointCOLUMN <- ""
timepointVALUE <- splitVALUES(specification$timepoint_values, "|")
covariateVALUE <- as.character(specification$covariates)
if (is.na(covariateVALUE)) covariateVALUE <- ""
agreementCOLUMN <- unique(c(
  as.character(specification$outcome),
  splitVALUES(covariateVALUE),
  if (option[["method"]] == "plink_ct") as.character(specification$prs_name) else
    paste0(as.character(specification$prs_name), "_", toupper(option[["method"]]))
))
agreementCOLUMN <- intersect(agreementCOLUMN, names(phenotype))

if (length(timepointVALUE) > 0L) {
  if (!nzchar(timepointCOLUMN) || !timepointCOLUMN %in% names(phenotype)) {
    stop(sprintf("Model '%s' requires phenotype timepoint column '%s'.", specification$model_id, timepointCOLUMN), call. = FALSE)
  }
  availableTIMEPOINT <- unique(as.character(phenotype[[timepointCOLUMN]]))
  availableTIMEPOINT <- availableTIMEPOINT[!is.na(availableTIMEPOINT) & nzchar(availableTIMEPOINT)]
  missingTIMEPOINT <- setdiff(timepointVALUE, availableTIMEPOINT)
  if (length(missingTIMEPOINT) > 0L) {
    stop(
      sprintf(
        "Model '%s' selects timepoint value(s) %s, but column '%s' contains %s. Select available values in timepoint_values.",
        specification$model_id,
        paste(missingTIMEPOINT, collapse = ", "),
        timepointCOLUMN,
        if (length(availableTIMEPOINT) > 0L) paste(availableTIMEPOINT, collapse = ", ") else "no non-missing values"
      ),
      call. = FALSE
    )
  }
  selectedCANDIDATE <- phenotype[as.character(get(timepointCOLUMN)) %in% timepointVALUE]
  recordKEY <- paste(selectedCANDIDATE[[participantCOLUMN]], selectedCANDIDATE[[timepointCOLUMN]], sep = "\r")
  duplicateKEY <- unique(recordKEY[duplicated(recordKEY) | duplicated(recordKEY, fromLast = TRUE)])
  for (key in duplicateKEY) {
    index <- which(recordKEY == key)
    for (column in agreementCOLUMN) {
      observedVALUE <- unique(as.character(selectedCANDIDATE[[column]][index]))
      observedVALUE <- observedVALUE[!is.na(observedVALUE)]
      if (length(observedVALUE) > 1L) {
        stop(
          sprintf(
            "Model '%s' has conflicting values for participant '%s', field '%s', timepoint '%s': %s. Correct the technical records or timepoint selection.",
            specification$model_id, selectedCANDIDATE[[participantCOLUMN]][index[[1L]]], column,
            selectedCANDIDATE[[timepointCOLUMN]][index[[1L]]], paste(observedVALUE, collapse = ", ")
          ),
          call. = FALSE
        )
      }
    }
  }
  selectedCANDIDATE[, technical_record_count := .N, by = c(participantCOLUMN, timepointCOLUMN)]
  data.table::setorderv(selectedCANDIDATE, "source_row")
  selectedPhenotype <- selectedCANDIDATE[!duplicated(recordKEY)]
  selectedPhenotype[, `:=`(
    selection_model_id = as.character(specification$model_id),
    requested_timepoint = paste(timepointVALUE, collapse = ","),
    selected_timepoint = as.character(get(timepointCOLUMN)),
    selection_reason = data.table::fifelse(
      technical_record_count > 1L,
      "selected_first_of_identical_technical_records",
      "selected_requested_timepoint"
    )
  )]

  participantORDER <- unique(phenotype[[participantCOLUMN]])
  completeness <- data.table::CJ(
    participant_index = seq_along(participantORDER),
    timepoint_index = seq_along(timepointVALUE),
    sorted = FALSE
  )
  completeness[, `:=`(
    model_id = as.character(specification$model_id),
    participant_id = participantORDER[participant_index],
    requested_timepoint = timepointVALUE[timepoint_index]
  )]
  selectedKEY <- paste(selectedPhenotype[[participantCOLUMN]], selectedPhenotype[[timepointCOLUMN]], sep = "\r")
  completenessKEY <- paste(completeness$participant_id, completeness$requested_timepoint, sep = "\r")
  selectedROW <- match(completenessKEY, selectedKEY)
  completeness[, `:=`(
    source_row = selectedPhenotype$source_row[selectedROW],
    status = data.table::fifelse(is.na(selectedROW), "MISSING", "SELECTED"),
    reason = data.table::fifelse(
      is.na(selectedROW),
      "Participant has no record at the requested timepoint.",
      "The first agreeing source record was selected."
    )
  )]
  completeness[, c("participant_index", "timepoint_index") := NULL]
} else {
  if (!groupedMODEL && anyDuplicated(phenotype[[participantCOLUMN]])) {
    stop(
      sprintf(
        "Fixed model '%s' has repeated participant IDs. Supply timepoint_column and one timepoint_values entry.",
        specification$model_id
      ),
      call. = FALSE
    )
  }
  selectedPhenotype <- data.table::copy(phenotype)
  selectedPhenotype[, `:=`(
    technical_record_count = 1L,
    selection_model_id = as.character(specification$model_id),
    requested_timepoint = "",
    selected_timepoint = "",
    selection_reason = "selected_without_timepoint_filter"
  )]
  completeness <- data.table::data.table(
    model_id = character(), participant_id = character(), requested_timepoint = character(),
    source_row = integer(), status = character(), reason = character()
  )
}
data.table::setorderv(selectedPhenotype, "source_row")
data.table::fwrite(selectedPhenotype, "phenotype_participant_level.tsv", sep = "\t", na = "NA")
data.table::fwrite(completeness, "phenotype_timepoint_completeness.tsv", sep = "\t", na = "NA")

emptyRESULT <- data.table::data.table(
  model_id = character(), cohort = character(), role = character(), trait_id = character(),
  prs_name = character(), outcome = character(), method = character(), family = character(),
  primary = logical(), participants = integer(), n = integer(), beta = numeric(),
  std_error = numeric(), ci_low = numeric(), ci_high = numeric(), parametric_p = numeric(),
  p_value = numeric(), r2_base = numeric(), r2_full = numeric(), delta_r2 = numeric(),
  partial_r2 = numeric(), null_fit = numeric(), full_fit = numeric(), incremental_fit = numeric(),
  fit_metric = character(), permutation_scheme = character(), permutations = integer(),
  permutation_p = numeric(), permutation_holm = numeric(), empirical_p = numeric(),
  expected_direction = character(), direction_match = character(), input_rows = integer(),
  primary_rows = integer(), complete_cases = integer(), excluded_non_primary = integer(),
  excluded_missing = integer(), convergence = character(), singular = character(),
  separation = character(), diagnostics = character(), status = character()
)
emptyMODEL <- data.table::data.table(
  model_id = character(), cohort = character(), method = character(), formula = character(),
  null_formula = character(), estimator = character(), expected_direction = character(), primary = logical(),
  convergence = character(), singular = character(), separation = character(), diagnostics = character(),
  status = character()
)
emptyPLOT <- data.table::data.table(
  model_id = character(), outcome = character(), cohort = character(), role = character(),
  trait_id = character(), prs_name = character(), method = character(), family = character(),
  estimator = character(), IID = character(), observed = numeric(), fitted_null = numeric(),
  fitted_full = numeric(), residual_full = numeric(), adjusted_outcome = numeric(), adjusted_prs = numeric()
)
emptyPERMUTATION <- data.table::data.table(
  model_id = character(), cohort = character(), role = character(), trait_id = character(),
  prs_name = character(), method = character(), family = character(), estimator = character(),
  permutation_id = integer(), permuted_beta = numeric(), observed_beta = numeric(),
  permutation_scheme = character(), status = character(), reason = character()
)
emptyINFLUENCE <- data.table::data.table(
  model_id = character(), cohort = character(), role = character(), trait_id = character(),
  prs_name = character(), method = character(), family = character(), estimator = character(),
  IID = character(), full_beta = numeric(), beta_without = numeric(), beta_change = numeric(),
  status = character(), reason = character()
)

if (nrow(modelSPEC) == 0L) {
  data.table::fwrite(emptyRESULT, "phenotype_associations.tsv", sep = "\t")
  data.table::fwrite(emptyMODEL, "phenotype_models_fitted.tsv", sep = "\t")
  data.table::fwrite(emptyPLOT, "phenotype_plot_data.tsv", sep = "\t")
  data.table::fwrite(emptyPERMUTATION, "phenotype_permutations.tsv", sep = "\t")
  data.table::fwrite(emptyINFLUENCE, "phenotype_influence.tsv", sep = "\t")
  quit(save = "no", status = 0L)
}

quoteNAME <- function(value) paste0("`", gsub("`", "", value, fixed = TRUE), "`")

# Preserve the legacy Freedman-Lane diagnostic while bounding work for larger cohorts.
# Exact permutations are used through n=8; larger models use a reproducible sample.
allORDERS <- function(value) {
  if (length(value) == 1L) return(matrix(value, nrow = 1L))
  do.call(rbind, lapply(seq_along(value), function(index) {
    cbind(value[index], allORDERS(value[-index]))
  }))
}

baseSEED <- suppressWarnings(as.integer(option[["seed"]]))
if (length(baseSEED) != 1L || is.na(baseSEED)) baseSEED <- 20260829L
diagnosticSEED <- function(...) {
  code <- utf8ToInt(paste(..., collapse = "|"))
  value <- (as.double(baseSEED) + sum(as.double(code) * seq_along(code))) %% 2147483646
  as.integer(value + 1)
}

resultLIST <- list()
fittedLIST <- list()
plotLIST <- list()
permutationLIST <- list()
influenceLIST <- list()
resultN <- 0L

for (modelROW in seq_len(nrow(modelSPEC))) {
  specification <- modelSPEC[modelROW]
  familyVALUE <- tolower(specification$family)
  covariateVALUE <- as.character(specification$covariates)
  if (length(covariateVALUE) != 1L || is.na(covariateVALUE)) covariateVALUE <- ""
  covariate <- trimws(strsplit(covariateVALUE, ",", fixed = TRUE)[[1L]])
  covariate <- covariate[covariate != ""]
  scoreSUBSET <- score[
    prs_name == specification$prs_name &
      (!nzchar(option[["cohort"]]) | cohort == option[["cohort"]]) &
      (!nzchar(option[["trait-id"]]) | trait_id == option[["trait-id"]]) &
      (!nzchar(option[["method"]]) | method == option[["method"]])
  ]
  selectedPhenotype[[specification$participant_id]] <- as.character(selectedPhenotype[[specification$participant_id]])

  for (cohortVALUE in unique(scoreSUBSET$cohort)) {
    for (methodVALUE in unique(scoreSUBSET[cohort == cohortVALUE, method])) {
      genetic <- scoreSUBSET[cohort == cohortVALUE & method == methodVALUE]
      genetic$IID <- as.character(genetic$IID)
      if (anyDuplicated(genetic$IID)) stop("A participant has more than one matching score for a model task.", call. = FALSE)
      analysis <- data.table::copy(selectedPhenotype)
      scoreROW <- match(analysis[[specification$participant_id]], genetic$IID)
      analysis[, `:=`(
        .PRS_Z = genetic$prs_z[scoreROW],
        .PRIMARY_ANALYSIS = genetic$primary_analysis[scoreROW],
        trait_id = genetic$trait_id[scoreROW],
        role = genetic$role[scoreROW]
      )]
      analysis <- analysis[!is.na(.PRS_Z)]
      inputROWS <- nrow(analysis)
      primaryROWS <- sum(analysis$.PRIMARY_ANALYSIS %in% TRUE)
      analysis <- analysis[.PRIMARY_ANALYSIS %in% TRUE]
      if (nrow(analysis) == 0L) stop("No participant is eligible for a primary phenotype analysis.", call. = FALSE)
      if (familyVALUE == "binomial") {
        observed <- as.character(analysis[[specification$outcome]])
        control <- as.character(specification$control_value)
        case <- as.character(specification$case_value)
        if (is.na(control)) control <- ""
        if (is.na(case)) case <- ""
        if (!nzchar(control) && !nzchar(case) && setequal(unique(observed[!is.na(observed)]), c("0", "1"))) {
          control <- "0"
          case <- "1"
        }
        analysis[[specification$outcome]] <- ifelse(
          observed == control,
          0,
          ifelse(observed == case, 1, NA_real_)
        )
      }
      needed <- unique(c(specification$outcome, ".PRS_Z", covariate, specification$group_id))
      needed <- needed[!is.na(needed) & needed != ""]
      analysis <- analysis[stats::complete.cases(analysis[, ..needed])]
      modelCOVARIATE <- covariate
      for (index in seq_along(covariate)) {
        column <- covariate[[index]]
        if (is.numeric(analysis[[column]])) {
          centredCOLUMN <- paste0(column, "_centered")
          if (centredCOLUMN %in% names(analysis) && centredCOLUMN != column) {
            stop(sprintf("Cannot create centred covariate '%s' because that column already exists.", centredCOLUMN), call. = FALSE)
          }
          analysis[[centredCOLUMN]] <- analysis[[column]] - mean(analysis[[column]])
          modelCOVARIATE[[index]] <- centredCOLUMN
        }
      }
      parameterN <- length(covariate) + 2L
      resultN <- resultN + 1L

      if (nrow(analysis) <= parameterN) {
        resultLIST[[resultN]] <- data.table::data.table(
          model_id = specification$model_id, cohort = cohortVALUE, role = unique(genetic$role),
          trait_id = unique(genetic$trait_id), prs_name = specification$prs_name,
          outcome = specification$outcome, method = methodVALUE, family = specification$family,
          primary = as.logical(specification$primary),
          participants = data.table::uniqueN(analysis[[specification$participant_id]]), n = nrow(analysis),
          beta = NA_real_, std_error = NA_real_, ci_low = NA_real_, ci_high = NA_real_,
          parametric_p = NA_real_, p_value = NA_real_, r2_base = NA_real_, r2_full = NA_real_,
          delta_r2 = NA_real_, partial_r2 = NA_real_, null_fit = NA_real_, full_fit = NA_real_,
          incremental_fit = NA_real_, fit_metric = NA_character_, permutation_scheme = "NOT_RUN",
          permutations = 0L, permutation_p = NA_real_, permutation_holm = NA_real_, empirical_p = NA_real_,
          expected_direction = specification$expected_direction, direction_match = "NOT_ESTIMATED",
          input_rows = inputROWS, primary_rows = primaryROWS, complete_cases = nrow(analysis),
          excluded_non_primary = inputROWS - primaryROWS, excluded_missing = primaryROWS - nrow(analysis),
          convergence = "NOT_ESTIMATED",
          singular = "NOT_ESTIMATED", separation = "NOT_ESTIMATED",
          diagnostics = "Too few complete cases for the declared model.",
          status = "INSUFFICIENT_COMPLETE_CASES"
        )
        next
      }

      rhsNULL <- if (length(modelCOVARIATE) == 0L) "1" else paste(quoteNAME(modelCOVARIATE), collapse = " + ")
      rhsFULL <- paste(c(if (length(modelCOVARIATE) > 0L) quoteNAME(modelCOVARIATE), ".PRS_Z"), collapse = " + ")
      grouped <- !is.na(specification$group_id) && specification$group_id != ""
      if (grouped) {
        rhsNULL <- paste(rhsNULL, paste0("(1 | ", quoteNAME(specification$group_id), ")"), sep = " + ")
        rhsFULL <- paste(rhsFULL, paste0("(1 | ", quoteNAME(specification$group_id), ")"), sep = " + ")
      }
      formulaNULL <- stats::as.formula(paste(quoteNAME(specification$outcome), "~", rhsNULL))
      formulaFULL <- stats::as.formula(paste(quoteNAME(specification$outcome), "~", rhsFULL))
      if (grouped && data.table::uniqueN(analysis[[specification$group_id]]) < 2L) {
        resultLIST[[resultN]] <- data.table::data.table(
          model_id = specification$model_id, cohort = cohortVALUE, role = unique(genetic$role),
          trait_id = unique(genetic$trait_id), prs_name = specification$prs_name,
          outcome = specification$outcome, method = methodVALUE, family = specification$family,
          primary = as.logical(specification$primary),
          participants = data.table::uniqueN(analysis[[specification$participant_id]]), n = nrow(analysis),
          beta = NA_real_, std_error = NA_real_, ci_low = NA_real_, ci_high = NA_real_,
          parametric_p = NA_real_, p_value = NA_real_, r2_base = NA_real_, r2_full = NA_real_,
          delta_r2 = NA_real_, partial_r2 = NA_real_, null_fit = NA_real_, full_fit = NA_real_,
          incremental_fit = NA_real_, fit_metric = NA_character_, permutation_scheme = "NOT_RUN",
          permutations = 0L, permutation_p = NA_real_, permutation_holm = NA_real_, empirical_p = NA_real_,
          expected_direction = specification$expected_direction, direction_match = "NOT_ESTIMATED",
          input_rows = inputROWS, primary_rows = primaryROWS, complete_cases = nrow(analysis),
          excluded_non_primary = inputROWS - primaryROWS, excluded_missing = primaryROWS - nrow(analysis),
          convergence = "NOT_ESTIMATED",
          singular = "NOT_ESTIMATED", separation = "NOT_ESTIMATED",
          diagnostics = "A grouped model requires at least two observed groups.",
          status = "INSUFFICIENT_GROUPS"
        )
        next
      }

      fitWARNINGS <- character()
      confidenceLOW <- NA_real_
      confidenceHIGH <- NA_real_
      withCallingHandlers({
      if (grouped) {
        if (!requireNamespace("lme4", quietly = TRUE)) stop("lme4 is required for declared grouped models.", call. = FALSE)
        if (familyVALUE == "gaussian") {
          nullFIT <- lme4::lmer(formulaNULL, data = analysis, REML = FALSE)
          fullFIT <- lme4::lmer(formulaFULL, data = analysis, REML = FALSE)
        } else {
          familyFUNCTION <- get(familyVALUE, envir = asNamespace("stats"))
          nullFIT <- lme4::glmer(formulaNULL, data = analysis, family = familyFUNCTION())
          fullFIT <- lme4::glmer(formulaFULL, data = analysis, family = familyFUNCTION())
        }
        beta <- unname(lme4::fixef(fullFIT)[".PRS_Z"])
        standardERROR <- unname(sqrt(diag(as.matrix(stats::vcov(fullFIT))))[".PRS_Z"])
        confidenceLOW <- beta - 1.96 * standardERROR
        confidenceHIGH <- beta + 1.96 * standardERROR
        pVALUE <- 2 * stats::pnorm(abs(beta / standardERROR), lower.tail = FALSE)
        nullVALUE <- stats::AIC(nullFIT)
        fullVALUE <- stats::AIC(fullFIT)
        incremental <- nullVALUE - fullVALUE
        estimator <- if (familyVALUE == "gaussian") "lme4::lmer" else "lme4::glmer"
        fitMETRIC <- "AIC reduction"
      } else if (familyVALUE == "gaussian") {
        nullFIT <- stats::lm(formulaNULL, data = analysis)
        fullFIT <- stats::lm(formulaFULL, data = analysis)
        coefficient <- summary(fullFIT)$coefficients[".PRS_Z", ]
        beta <- unname(coefficient["Estimate"])
        standardERROR <- unname(coefficient["Std. Error"])
        pVALUE <- unname(coefficient["Pr(>|t|)"])
        interval <- stats::confint(fullFIT, ".PRS_Z", level = 0.95)
        confidenceLOW <- unname(interval[[1L]])
        confidenceHIGH <- unname(interval[[2L]])
        nullVALUE <- summary(nullFIT)$r.squared
        fullVALUE <- summary(fullFIT)$r.squared
        incremental <- fullVALUE - nullVALUE
        estimator <- "stats::lm"
        fitMETRIC <- "Delta R-squared"
      } else {
        if (!requireNamespace("glm2", quietly = TRUE)) stop("glm2 is required for declared generalised models.", call. = FALSE)
        familyFUNCTION <- get(familyVALUE, envir = asNamespace("stats"))
        nullFIT <- glm2::glm2(formulaNULL, data = analysis, family = familyFUNCTION())
        fullFIT <- glm2::glm2(formulaFULL, data = analysis, family = familyFUNCTION())
        coefficient <- summary(fullFIT)$coefficients[".PRS_Z", ]
        beta <- unname(coefficient["Estimate"])
        standardERROR <- unname(coefficient["Std. Error"])
        pVALUE <- unname(coefficient[ncol(summary(fullFIT)$coefficients)])
        confidenceLOW <- beta - 1.96 * standardERROR
        confidenceHIGH <- beta + 1.96 * standardERROR
        nullVALUE <- 1 - nullFIT$deviance / nullFIT$null.deviance
        fullVALUE <- 1 - fullFIT$deviance / fullFIT$null.deviance
        incremental <- fullVALUE - nullVALUE
        estimator <- "glm2::glm2"
        fitMETRIC <- "Delta pseudo-R-squared"
      }
      }, warning = function(condition) {
        fitWARNINGS <<- unique(c(fitWARNINGS, conditionMessage(condition)))
        invokeRestart("muffleWarning")
      })
      partialR2 <- if (!grouped && familyVALUE == "gaussian" && is.finite(nullVALUE) && nullVALUE < 1) {
        incremental / (1 - nullVALUE)
      } else {
        NA_real_
      }

      convergenceSTATUS <- "PASS"
      singularSTATUS <- "NOT_APPLICABLE"
      separationSTATUS <- "NOT_APPLICABLE"
      diagnosticMESSAGE <- fitWARNINGS
      if (grouped) {
        convergenceMESSAGE <- unlist(fullFIT@optinfo$conv$lme4$messages, use.names = FALSE)
        convergenceMESSAGE <- convergenceMESSAGE[!is.na(convergenceMESSAGE) & nzchar(convergenceMESSAGE)]
        if (length(convergenceMESSAGE) > 0L) {
          convergenceSTATUS <- "REVIEW"
          diagnosticMESSAGE <- c(diagnosticMESSAGE, paste(convergenceMESSAGE, collapse = "; "))
        }
        singularSTATUS <- if (lme4::isSingular(fullFIT, tol = 1e-4)) "REVIEW" else "PASS"
        if (singularSTATUS == "REVIEW") diagnosticMESSAGE <- c(diagnosticMESSAGE, "The fitted mixed model is singular.")
      } else if (familyVALUE == "gaussian") {
        singularSTATUS <- if (fullFIT$rank < ncol(stats::model.matrix(fullFIT))) "REVIEW" else "PASS"
        if (singularSTATUS == "REVIEW") diagnosticMESSAGE <- c(diagnosticMESSAGE, "The fixed-effects design is rank deficient.")
      } else if (!isTRUE(fullFIT$converged)) {
        convergenceSTATUS <- "REVIEW"
        diagnosticMESSAGE <- c(diagnosticMESSAGE, "The generalised model did not converge.")
      }
      if (familyVALUE == "binomial") {
        fittedPROBABILITY <- as.numeric(stats::fitted(fullFIT))
        coefficientVALUE <- if (grouped) lme4::fixef(fullFIT) else stats::coef(fullFIT)
        separationSTATUS <- if (
          any(abs(coefficientVALUE) > 10, na.rm = TRUE) ||
            any(fittedPROBABILITY < 1e-8 | fittedPROBABILITY > 1 - 1e-8, na.rm = TRUE)
        ) "REVIEW" else "PASS"
        if (separationSTATUS == "REVIEW") {
          diagnosticMESSAGE <- c(diagnosticMESSAGE, "Extreme coefficients or fitted probabilities indicate possible separation.")
        }
      }
      finiteESTIMATE <- all(is.finite(c(beta, standardERROR, pVALUE, nullVALUE, fullVALUE, incremental)))
      resultSTATUS <- if (!finiteESTIMATE) {
        diagnosticMESSAGE <- c(diagnosticMESSAGE, "One or more reported model estimates are non-finite.")
        "NON_FINITE"
      } else if (any(c(convergenceSTATUS, singularSTATUS, separationSTATUS) == "REVIEW") || length(fitWARNINGS) > 0L) {
        "REVIEW_DIAGNOSTICS"
      } else {
        "ESTIMATED"
      }
      diagnosticMESSAGE <- paste(unique(diagnosticMESSAGE[nzchar(diagnosticMESSAGE)]), collapse = " | ")
      if (!nzchar(diagnosticMESSAGE)) diagnosticMESSAGE <- "No model-fit warnings were detected."

      permutationSCHEME <- "NOT_RUN"
      permutationN <- 0L
      empiricalP <- NA_real_
      diagnosticREASON <- "Residual-permutation and case-deletion diagnostics are available for fixed-effects Gaussian models only."
      diagnosticKEY <- list(
        model_id = specification$model_id,
        cohort = cohortVALUE,
        role = unique(genetic$role),
        trait_id = unique(genetic$trait_id),
        prs_name = specification$prs_name,
        method = methodVALUE,
        family = familyVALUE,
        estimator = estimator
      )
      diagnosticROW <- data.table::as.data.table(diagnosticKEY)

      if (!grouped && familyVALUE == "gaussian") {
        observationN <- nrow(analysis)
        if (observationN <= 8L) {
          orderMATRIX <- allORDERS(seq_len(observationN))
          permutationSCHEME <- "Freedman-Lane exact"
        } else {
          permutationN <- 10000L
          set.seed(diagnosticSEED(specification$model_id, cohortVALUE, methodVALUE))
          orderMATRIX <- t(replicate(permutationN, sample.int(observationN), simplify = "matrix"))
          permutationSCHEME <- "Freedman-Lane Monte Carlo"
        }
        permutationN <- nrow(orderMATRIX)
        nullRESIDUAL <- as.numeric(stats::residuals(nullFIT))
        nullFITTED <- as.numeric(stats::fitted(nullFIT))
        designMATRIX <- stats::model.matrix(fullFIT)
        designQR <- qr(designMATRIX)
        coefficientINDEX <- match(".PRS_Z", colnames(designMATRIX))
        permutedBETA <- vapply(seq_len(permutationN), function(permutationID) {
          coefficient <- qr.coef(
            designQR,
            nullFITTED + nullRESIDUAL[orderMATRIX[permutationID, ]]
          )
          unname(coefficient[coefficientINDEX])
        }, numeric(1L))
        finitePERMUTATION <- is.finite(permutedBETA)
        validN <- sum(finitePERMUTATION)
        if (validN > 0L) {
          exceedN <- sum(abs(permutedBETA[finitePERMUTATION]) >= abs(beta) - 1e-12)
          empiricalP <- if (identical(permutationSCHEME, "Freedman-Lane exact")) {
            exceedN / validN
          } else {
            (exceedN + 1) / (validN + 1)
          }
        }
        permutationVALUE <- diagnosticROW[rep(1L, permutationN)]
        permutationVALUE[, `:=`(
          permutation_id = seq_len(permutationN),
          permuted_beta = permutedBETA,
          observed_beta = beta,
          permutation_scheme = permutationSCHEME,
          status = data.table::fifelse(finitePERMUTATION, "ESTIMATED", "NON_FINITE"),
          reason = "Null-model residuals were permuted and added to null fitted values before refitting the full design."
        )]
        permutationLIST[[length(permutationLIST) + 1L]] <- permutationVALUE

        deletionCHANGE <- as.numeric(stats::dfbeta(fullFIT)[, ".PRS_Z"])
        deletionBETA <- beta - deletionCHANGE
        influenceVALUE <- diagnosticROW[rep(1L, nrow(analysis))]
        influenceVALUE[, `:=`(
          IID = as.character(analysis[[specification$participant_id]]),
          full_beta = beta,
          beta_without = deletionBETA,
          beta_change = deletionBETA - beta,
          status = data.table::fifelse(is.finite(deletionBETA), "ESTIMATED", "NON_FINITE"),
          reason = "Coefficient after deleting the named participant, calculated from the exact OLS case-deletion identity."
        )]
        influenceLIST[[length(influenceLIST) + 1L]] <- influenceVALUE
      } else {
        permutationLIST[[length(permutationLIST) + 1L]] <- data.table::copy(diagnosticROW)[, `:=`(
          permutation_id = NA_integer_, permuted_beta = NA_real_, observed_beta = beta,
          permutation_scheme = permutationSCHEME, status = "NOT_RUN", reason = diagnosticREASON
        )]
        influenceLIST[[length(influenceLIST) + 1L]] <- data.table::copy(diagnosticROW)[, `:=`(
          IID = NA_character_, full_beta = beta, beta_without = NA_real_, beta_change = NA_real_,
          status = "NOT_RUN", reason = diagnosticREASON
        )]
      }

      expectedDIRECTION <- trimws(as.character(specification$expected_direction))
      if (is.na(expectedDIRECTION) || !nzchar(expectedDIRECTION)) expectedDIRECTION <- "not_declared"
      directionMATCH <- if (expectedDIRECTION == "not_declared") {
        "NOT_DECLARED"
      } else if ((expectedDIRECTION == "positive" && beta > 0) || (expectedDIRECTION == "negative" && beta < 0)) {
        "MATCH"
      } else {
        "MISMATCH"
      }
      resultLIST[[resultN]] <- data.table::data.table(
        model_id = specification$model_id, cohort = cohortVALUE, role = unique(genetic$role),
        trait_id = unique(genetic$trait_id), prs_name = specification$prs_name,
        outcome = specification$outcome, method = methodVALUE, family = familyVALUE,
        primary = as.logical(specification$primary),
        participants = data.table::uniqueN(analysis[[specification$participant_id]]), n = nrow(analysis),
        beta, std_error = standardERROR, ci_low = confidenceLOW, ci_high = confidenceHIGH,
        parametric_p = pVALUE, p_value = pVALUE, r2_base = nullVALUE, r2_full = fullVALUE,
        delta_r2 = incremental, partial_r2 = partialR2, null_fit = nullVALUE,
        full_fit = fullVALUE, incremental_fit = incremental, fit_metric = fitMETRIC,
        permutation_scheme = permutationSCHEME, permutations = permutationN,
        permutation_p = empiricalP, permutation_holm = NA_real_, empirical_p = empiricalP,
        expected_direction = expectedDIRECTION, direction_match = directionMATCH,
        input_rows = inputROWS, primary_rows = primaryROWS, complete_cases = nrow(analysis),
        excluded_non_primary = inputROWS - primaryROWS, excluded_missing = primaryROWS - nrow(analysis),
        convergence = convergenceSTATUS, singular = singularSTATUS,
        separation = separationSTATUS, diagnostics = diagnosticMESSAGE, status = resultSTATUS
      )
      fittedLIST[[length(fittedLIST) + 1L]] <- data.table::data.table(
        model_id = specification$model_id, cohort = cohortVALUE, method = methodVALUE,
        formula = paste(deparse(formulaFULL), collapse = ""),
        null_formula = paste(deparse(formulaNULL), collapse = ""), estimator,
        expected_direction = specification$expected_direction,
        primary = as.logical(specification$primary), convergence = convergenceSTATUS,
        singular = singularSTATUS, separation = separationSTATUS,
        diagnostics = diagnosticMESSAGE, status = resultSTATUS
      )

      adjustedOUTCOME <- rep(NA_real_, nrow(analysis))
      adjustedPRS <- rep(NA_real_, nrow(analysis))
      if (!grouped && familyVALUE == "gaussian") {
        adjustedOUTCOME <- stats::residuals(nullFIT)
        prsFORMULA <- stats::as.formula(paste(".PRS_Z ~", rhsNULL))
        adjustedPRS <- stats::residuals(stats::lm(prsFORMULA, data = analysis))
      }
      plotLIST[[length(plotLIST) + 1L]] <- data.table::data.table(
        model_id = specification$model_id,
        outcome = specification$outcome,
        cohort = cohortVALUE,
        role = unique(genetic$role),
        trait_id = unique(genetic$trait_id),
        prs_name = specification$prs_name,
        method = methodVALUE,
        family = familyVALUE,
        estimator = estimator,
        IID = as.character(analysis[[specification$participant_id]]),
        observed = as.numeric(stats::model.response(stats::model.frame(fullFIT))),
        fitted_null = as.numeric(stats::fitted(nullFIT)),
        fitted_full = as.numeric(stats::fitted(fullFIT)),
        residual_full = as.numeric(stats::residuals(fullFIT, type = "response")),
        adjusted_outcome = as.numeric(adjustedOUTCOME),
        adjusted_prs = as.numeric(adjustedPRS)
      )
    }
  }
}

result <- if (length(resultLIST) == 0L) emptyRESULT else data.table::rbindlist(resultLIST, use.names = TRUE, fill = TRUE)
fitted <- if (length(fittedLIST) == 0L) emptyMODEL else data.table::rbindlist(fittedLIST, use.names = TRUE, fill = TRUE)
plotDATA <- if (length(plotLIST) == 0L) emptyPLOT else data.table::rbindlist(plotLIST, use.names = TRUE, fill = TRUE)
permutationDATA <- if (length(permutationLIST) == 0L) emptyPERMUTATION else data.table::rbindlist(permutationLIST, use.names = TRUE, fill = TRUE)
influenceDATA <- if (length(influenceLIST) == 0L) emptyINFLUENCE else data.table::rbindlist(influenceLIST, use.names = TRUE, fill = TRUE)
data.table::fwrite(result, "phenotype_associations.tsv", sep = "\t", na = "NA")
data.table::fwrite(fitted, "phenotype_models_fitted.tsv", sep = "\t", na = "NA")
data.table::fwrite(plotDATA, "phenotype_plot_data.tsv", sep = "\t", na = "NA")
data.table::fwrite(permutationDATA, "phenotype_permutations.tsv", sep = "\t", na = "NA")
data.table::fwrite(influenceDATA, "phenotype_influence.tsv", sep = "\t", na = "NA")
