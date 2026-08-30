#!/usr/bin/env Rscript

# Parse the declared GWAS column mapping.
argument <- commandArgs(trailingOnly = TRUE)
name <- sub("^--", "", argument[seq.int(1L, length(argument), 2L)])
value <- argument[seq.int(2L, length(argument), 2L)]
option <- stats::setNames(as.list(value), name)

if (!requireNamespace("data.table", quietly = TRUE)) {
  stop("The data.table package is required.", call. = FALSE)
}

requiredOPTION <- c(
  "input", "trait-id", "prs-name", "effect-type", "sample-size", "snp-col", "chr-col",
  "bp-col", "effect-allele-col", "other-allele-col", "beta-col", "se-col", "p-col",
  "freq-col", "n-col"
)
missingOPTION <- setdiff(requiredOPTION, names(option))
if (length(missingOPTION) > 0L) {
  stop(sprintf("Missing option(s): %s", paste(missingOPTION, collapse = ", ")), call. = FALSE)
}

# Read only the columns explicitly declared in the manifest.
sourceCOLUMN <- unique(c(
  option[["snp-col"]], option[["chr-col"]], option[["bp-col"]],
  option[["effect-allele-col"]], option[["other-allele-col"]], option[["beta-col"]],
  option[["se-col"]], option[["p-col"]], option[["freq-col"]], option[["n-col"]]
))
sourceCOLUMN <- sourceCOLUMN[!is.na(sourceCOLUMN) & sourceCOLUMN != ""]
gwas <- data.table::fread(option[["input"]], select = sourceCOLUMN, showProgress = FALSE)

extractCOLUMN <- function(column, default = NA) {
  if (is.null(column) || is.na(column) || column == "") rep(default, nrow(gwas)) else gwas[[column]]
}

effect <- suppressWarnings(as.numeric(extractCOLUMN(option[["beta-col"]])))
effectTYPE <- tolower(option[["effect-type"]])
if (effectTYPE == "or") {
  if (any(!is.finite(effect) | effect <= 0)) {
    stop("Odds ratios must be finite and greater than zero.", call. = FALSE)
  }
  effect <- log(effect)
} else if (!effectTYPE %in% c("beta", "log_or")) {
  stop("effect_type must be beta, log_or, or or.", call. = FALSE)
}

cojo <- data.table::data.table(
  SNP = as.character(extractCOLUMN(option[["snp-col"]])),
  CHR = sub("^chr", "", as.character(extractCOLUMN(option[["chr-col"]])), ignore.case = TRUE),
  BP = suppressWarnings(as.integer(extractCOLUMN(option[["bp-col"]]))),
  A1 = toupper(as.character(extractCOLUMN(option[["effect-allele-col"]]))),
  A2 = toupper(as.character(extractCOLUMN(option[["other-allele-col"]]))),
  freq = suppressWarnings(as.numeric(extractCOLUMN(option[["freq-col"]]))),
  b = effect,
  se = suppressWarnings(as.numeric(extractCOLUMN(option[["se-col"]]))),
  p = suppressWarnings(as.numeric(extractCOLUMN(option[["p-col"]]))),
  N = suppressWarnings(as.numeric(extractCOLUMN(option[["n-col"]], option[["sample-size"]])))
)

# Stop on structural errors instead of silently redoing completed GWAS QC.
validALLELE <- grepl("^[ACGT]+$", cojo$A1) & grepl("^[ACGT]+$", cojo$A2) & cojo$A1 != cojo$A2
validROW <- !is.na(cojo$SNP) & cojo$SNP != "" & cojo$CHR %in% as.character(1:22) &
  is.finite(cojo$BP) & cojo$BP > 0 & validALLELE & is.finite(cojo$b) &
  is.finite(cojo$se) & cojo$se > 0 & is.finite(cojo$p) & cojo$p > 0 & cojo$p <= 1 &
  is.finite(cojo$N) & cojo$N > 0
if (!all(validROW)) {
  stop(sprintf("GWAS '%s' contains %s structurally invalid row(s).", option[["trait-id"]], sum(!validROW)), call. = FALSE)
}
if (anyDuplicated(cojo$SNP)) {
  stop(sprintf("GWAS '%s' contains duplicated SNP identifiers.", option[["trait-id"]]), call. = FALSE)
}
if (all(is.na(cojo$freq))) {
  cojo$freq <- NA_real_
} else if (any(!is.finite(cojo$freq) | cojo$freq <= 0 | cojo$freq >= 1)) {
  stop(sprintf("GWAS '%s' contains invalid effect-allele frequencies.", option[["trait-id"]]), call. = FALSE)
}

data.table::setorder(cojo, CHR, BP)
data.table::fwrite(cojo, paste0(option[["trait-id"]], ".cojo.ma"), sep = "\t", quote = FALSE, na = "NA")
data.table::fwrite(
  cojo[, .(ID = SNP, CHR, POS = BP, A1, P = p)],
  paste0(option[["trait-id"]], ".clump.tsv"),
  sep = "\t",
  quote = FALSE
)
data.table::fwrite(
  data.table::data.table(
    trait_id = option[["trait-id"]],
    prs_name = option[["prs-name"]],
    source_variants = nrow(gwas),
    harmonised_variants = nrow(cojo),
    duplicated_snp = 0L,
    structural_status = "PASS"
  ),
  paste0(option[["trait-id"]], ".harmonisation_qc.tsv"),
  sep = "\t"
)
