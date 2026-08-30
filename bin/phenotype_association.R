#!/usr/bin/env Rscript

# Evaluate prespecified matching phenotype-PRS models without tuning score construction.
argument <- commandArgs(trailingOnly = TRUE)
name <- sub("^--", "", argument[seq.int(1L, length(argument), 2L)])
value <- argument[seq.int(2L, length(argument), 2L)]
option <- stats::setNames(as.list(value), name)
if (!requireNamespace("data.table", quietly = TRUE)) stop("The data.table package is required.", call. = FALSE)

score <- data.table::fread(option[["scores"]], colClasses = list(character = c("FID", "IID")))
phenotype <- data.table::fread(option[["phenotype"]])
modelSPEC <- data.table::fread(option[["models"]])

emptyRESULT <- data.table::data.table(
  model_id = character(), cohort = character(), role = character(), trait_id = character(),
  prs_name = character(), method = character(), family = character(), n = integer(), beta = numeric(),
  std_error = numeric(), ci_low = numeric(), ci_high = numeric(), p_value = numeric(),
  null_fit = numeric(), full_fit = numeric(), incremental_fit = numeric(), fit_metric = character(),
  status = character()
)
emptyMODEL <- data.table::data.table(
  model_id = character(), cohort = character(), method = character(), formula = character(),
  null_formula = character(), estimator = character(), expected_direction = character(), primary = logical()
)
emptyPLOT <- data.table::data.table(
  model_id = character(), outcome = character(), cohort = character(), role = character(),
  trait_id = character(), prs_name = character(), method = character(), family = character(),
  estimator = character(), IID = character(), observed = numeric(), fitted_null = numeric(),
  fitted_full = numeric(), residual_full = numeric(), adjusted_outcome = numeric(), adjusted_prs = numeric()
)

if (nrow(modelSPEC) == 0L) {
  data.table::fwrite(emptyRESULT, "phenotype_associations.tsv", sep = "\t")
  data.table::fwrite(emptyMODEL, "phenotype_models_fitted.tsv", sep = "\t")
  data.table::fwrite(emptyPLOT, "phenotype_plot_data.tsv", sep = "\t")
  quit(save = "no", status = 0L)
}

quoteNAME <- function(value) paste0("`", gsub("`", "", value, fixed = TRUE), "`")
resultLIST <- list()
fittedLIST <- list()
plotLIST <- list()
resultN <- 0L

for (modelROW in seq_len(nrow(modelSPEC))) {
  specification <- modelSPEC[modelROW]
  covariate <- trimws(strsplit(specification$covariates, ",", fixed = TRUE)[[1L]])
  covariate <- covariate[covariate != ""]
  scoreSUBSET <- score[prs_name == specification$prs_name]
  phenotype[[specification$participant_id]] <- as.character(phenotype[[specification$participant_id]])

  for (cohortVALUE in unique(scoreSUBSET$cohort)) {
    for (methodVALUE in unique(scoreSUBSET[cohort == cohortVALUE, method])) {
      genetic <- scoreSUBSET[cohort == cohortVALUE & method == methodVALUE]
      genetic$IID <- as.character(genetic$IID)
      analysis <- merge(
        phenotype,
        genetic[, .(IID, .PRS_Z = prs_z, trait_id, role)],
        by.x = specification$participant_id,
        by.y = "IID",
        all = FALSE
      )
      needed <- unique(c(specification$outcome, ".PRS_Z", covariate, specification$group_id))
      needed <- needed[!is.na(needed) & needed != ""]
      analysis <- analysis[stats::complete.cases(analysis[, ..needed])]
      parameterN <- length(covariate) + 2L
      resultN <- resultN + 1L

      if (nrow(analysis) <= parameterN) {
        resultLIST[[resultN]] <- data.table::data.table(
          model_id = specification$model_id, cohort = cohortVALUE, role = unique(genetic$role),
          trait_id = unique(genetic$trait_id), prs_name = specification$prs_name, method = methodVALUE,
          family = specification$family, n = nrow(analysis), beta = NA_real_, std_error = NA_real_,
          ci_low = NA_real_, ci_high = NA_real_, p_value = NA_real_, null_fit = NA_real_,
          full_fit = NA_real_, incremental_fit = NA_real_, fit_metric = NA_character_,
          status = "INSUFFICIENT_COMPLETE_CASES"
        )
        next
      }

      rhsNULL <- if (length(covariate) == 0L) "1" else paste(quoteNAME(covariate), collapse = " + ")
      rhsFULL <- paste(c(".PRS_Z", if (length(covariate) > 0L) quoteNAME(covariate)), collapse = " + ")
      grouped <- !is.na(specification$group_id) && specification$group_id != ""
      if (grouped) {
        rhsNULL <- paste(rhsNULL, paste0("(1 | ", quoteNAME(specification$group_id), ")"), sep = " + ")
        rhsFULL <- paste(rhsFULL, paste0("(1 | ", quoteNAME(specification$group_id), ")"), sep = " + ")
      }
      formulaNULL <- stats::as.formula(paste(quoteNAME(specification$outcome), "~", rhsNULL))
      formulaFULL <- stats::as.formula(paste(quoteNAME(specification$outcome), "~", rhsFULL))
      familyVALUE <- tolower(specification$family)

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
        nullVALUE <- 1 - nullFIT$deviance / nullFIT$null.deviance
        fullVALUE <- 1 - fullFIT$deviance / fullFIT$null.deviance
        incremental <- fullVALUE - nullVALUE
        estimator <- "glm2::glm2"
        fitMETRIC <- "Delta pseudo-R-squared"
      }

      resultLIST[[resultN]] <- data.table::data.table(
        model_id = specification$model_id, cohort = cohortVALUE, role = unique(genetic$role),
        trait_id = unique(genetic$trait_id), prs_name = specification$prs_name, method = methodVALUE,
        family = familyVALUE, n = nrow(analysis), beta, std_error = standardERROR,
        ci_low = beta - 1.96 * standardERROR, ci_high = beta + 1.96 * standardERROR,
        p_value = pVALUE, null_fit = nullVALUE, full_fit = fullVALUE,
        incremental_fit = incremental, fit_metric = fitMETRIC, status = "ESTIMATED"
      )
      fittedLIST[[length(fittedLIST) + 1L]] <- data.table::data.table(
        model_id = specification$model_id, cohort = cohortVALUE, method = methodVALUE,
        formula = paste(deparse(formulaFULL), collapse = ""),
        null_formula = paste(deparse(formulaNULL), collapse = ""), estimator,
        expected_direction = specification$expected_direction,
        primary = as.logical(specification$primary)
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
data.table::fwrite(result, "phenotype_associations.tsv", sep = "\t", na = "NA")
data.table::fwrite(fitted, "phenotype_models_fitted.tsv", sep = "\t", na = "NA")
data.table::fwrite(plotDATA, "phenotype_plot_data.tsv", sep = "\t", na = "NA")
