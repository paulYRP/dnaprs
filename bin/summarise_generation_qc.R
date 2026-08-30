#!/usr/bin/env Rscript

# Combine method-specific variant counts into one ordered PRS-generation record.
argument <- commandArgs(trailingOnly = TRUE)
name <- sub("^--", "", argument[seq.int(1L, length(argument), 2L)])
value <- argument[seq.int(2L, length(argument), 2L)]
option <- stats::setNames(as.list(value), name)
if (!requireNamespace("data.table", quietly = TRUE)) stop("The data.table package is required.", call. = FALSE)

qcFILE <- strsplit(option[["qc-files"]], ",", fixed = TRUE)[[1L]]
qcTABLE <- lapply(qcFILE, data.table::fread)
scoreQC <- data.table::fread(option[["score-qc"]])

selectQC <- function(required) {
  selected <- qcTABLE[vapply(qcTABLE, function(value) all(required %in% names(value)), logical(1L))]
  if (length(selected) == 0L) return(data.table::data.table())
  unique(data.table::rbindlist(selected, use.names = TRUE, fill = TRUE))
}

harmonisedQC <- selectQC(c("trait_id", "source_variants", "harmonised_variants"))
plinkQC <- selectQC(c("trait_id", "clumped_variants"))
tidyQC <- selectQC(c("trait_id", "input_variants", "ld_aligned_variants"))
imputeQC <- selectQC(c("trait_id", "summary_imputed_variants"))
modelQC <- selectQC(c("trait_id", "weights", "nonzero_weights"))

if (nrow(harmonisedQC) == 0L) stop("No GWAS harmonisation QC table was supplied.", call. = FALSE)
if (nrow(scoreQC) == 0L) stop("No final score QC row was supplied.", call. = FALSE)

oneROW <- function(value, traitID, label) {
  selected <- value[trait_id == traitID]
  if (nrow(selected) != 1L) {
    stop(sprintf("Expected one %s QC row for trait '%s'.", label, traitID), call. = FALSE)
  }
  selected
}

flowLIST <- vector("list", nrow(scoreQC))
for (row in seq_len(nrow(scoreQC))) {
  scoreROW <- scoreQC[row]
  traitID <- scoreROW$trait_id
  harmonisedROW <- oneROW(harmonisedQC, traitID, "harmonisation")

  if (scoreROW$method == "plink_ct") {
    methodROW <- oneROW(plinkQC, traitID, "PLINK C+T")
    stage <- c("Source GWAS", "Harmonised", "LD-clumped", "Scored")
    count <- c(
      harmonisedROW$source_variants,
      harmonisedROW$harmonised_variants,
      methodROW$clumped_variants,
      scoreROW$used_variants
    )
  } else if (scoreROW$method == "sbayesrc") {
    tidyROW <- oneROW(tidyQC, traitID, "SBayesRC tidy")
    imputeROW <- oneROW(imputeQC, traitID, "SBayesRC imputation")
    modelROW <- oneROW(modelQC, traitID, "SBayesRC model")
    stage <- c(
      "Source GWAS", "Harmonised", "LD-aligned", "Summary imputed",
      "Model weights", "Non-zero effects", "Scored"
    )
    count <- c(
      harmonisedROW$source_variants,
      harmonisedROW$harmonised_variants,
      tidyROW$ld_aligned_variants,
      imputeROW$summary_imputed_variants,
      modelROW$weights,
      modelROW$nonzero_weights,
      scoreROW$used_variants
    )
  } else {
    stop(sprintf("Unsupported PRS method in score QC: %s", scoreROW$method), call. = FALSE)
  }

  count <- as.numeric(count)
  if (any(!is.finite(count) | count <= 0)) {
    stop(sprintf("Trait '%s' has a missing or non-positive variant count.", traitID), call. = FALSE)
  }
  flowLIST[[row]] <- data.table::data.table(
    cohort = scoreROW$cohort,
    role = scoreROW$role,
    trait_id = traitID,
    prs_name = scoreROW$prs_name,
    method = scoreROW$method,
    stage = stage,
    stage_order = seq_along(stage),
    variant_count = count,
    percent_of_source = 100 * count / count[1L],
    percent_of_previous = 100 * count / c(count[1L], count[-length(count)])
  )
}

variantFLOW <- data.table::rbindlist(flowLIST)
data.table::setorder(variantFLOW, cohort, trait_id, method, stage_order)
data.table::fwrite(variantFLOW, "variant_flow.tsv", sep = "\t", na = "NA")
