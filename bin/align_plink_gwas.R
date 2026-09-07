#!/usr/bin/env Rscript

argument <- commandArgs(trailingOnly = TRUE)
if (length(argument) %% 2L != 0L) stop("Arguments must be --name value pairs.", call. = FALSE)
option <- stats::setNames(
  as.list(argument[seq.int(2L, length(argument), 2L)]),
  sub("^--", "", argument[seq.int(1L, length(argument), 2L)])
)
if (!requireNamespace("data.table", quietly = TRUE)) stop("The data.table package is required.", call. = FALSE)

complement <- function(value) chartr("ACGT", "TGCA", toupper(value))
canonicalPAIR <- function(first, second) {
  direct <- paste(pmin(first, second), pmax(first, second), sep = "/")
  reverse <- paste(pmin(complement(first), complement(second)), pmax(complement(first), complement(second)), sep = "/")
  pmin(direct, reverse)
}

data.table::setDTthreads(as.integer(option[["threads"]]))
gwas <- data.table::fread(option[["cojo"]])
index <- data.table::fread(option[["compatibility-index"]])
requiredGWAS <- c("SNP", "CHR", "BP", "A1", "A2", "freq", "b", "se", "p", "N")
if (!all(requiredGWAS %in% names(gwas))) stop("Harmonised GWAS is missing required COJO columns.", call. = FALSE)

gwas[, `:=`(
  CHR = sub("^chr", "", as.character(CHR), ignore.case = TRUE),
  BP = as.integer(BP), A1 = toupper(A1), A2 = toupper(A2), source_row = .I
)]
gwas[, canonical_pair := canonicalPAIR(A1, A2)]
gwas[, canonical_key := paste(CHR, BP, canonical_pair, sep = ":")]

palindromic <- paste0(pmin(gwas$A1, gwas$A2), pmax(gwas$A1, gwas$A2)) %in% c("AT", "CG")
positiveP <- is.finite(gwas$p) & gwas$p > 0
eligible <- gwas[!palindromic & positiveP]
eligibleROWS <- eligible$source_row

targetDUPLICATEKEY <- index[target_key_count != 1L, canonical_key]
targetUNIQUE <- index[target_key_count == 1L]
selected <- merge(
  eligible,
  targetUNIQUE[, .(canonical_key, target_id, target_ref, target_alt)],
  by = "canonical_key", all = FALSE, sort = FALSE
)
selected[, orientation := data.table::fcase(
  A1 == target_alt & A2 == target_ref, "DIRECT_ALT",
  A1 == target_ref & A2 == target_alt, "DIRECT_REF",
  complement(A1) == target_alt & complement(A2) == target_ref, "COMPLEMENT_ALT",
  complement(A1) == target_ref & complement(A2) == target_alt, "COMPLEMENT_REF",
  default = "INCOMPATIBLE"
)]
selected <- selected[orientation != "INCOMPATIBLE"]
selected[, `:=`(
  b = ifelse(orientation %in% c("DIRECT_REF", "COMPLEMENT_REF"), -b, b),
  freq = ifelse(orientation %in% c("DIRECT_REF", "COMPLEMENT_REF"), 1 - freq, freq),
  A1 = target_alt,
  A2 = target_ref,
  SNP = target_id
)]

targetALIGNED <- unique(selected$source_row)
targetAMBIGUOUS <- eligible[canonical_key %in% targetDUPLICATEKEY, source_row]
selected <- selected[!source_row %in% targetAMBIGUOUS]
selected <- selected[, .SD[1L], by = source_row]

referenceDUPLICATEKEY <- index[reference_key_count > 1L, canonical_key]
referenceUNIQUE <- index[reference_key_count == 1L]
referenceAMBIGUOUS <- selected[canonical_key %in% referenceDUPLICATEKEY, source_row]
selected <- selected[!source_row %in% referenceAMBIGUOUS]
selected <- merge(
  selected,
  referenceUNIQUE[, .(canonical_key, reference_id, reference_ref, reference_alt)],
  by = "canonical_key", all = FALSE, sort = FALSE
)
selected <- selected[target_ref == reference_ref & target_alt == reference_alt]
referenceALIGNED <- unique(selected$source_row)
selected <- selected[, .SD[1L], by = source_row]
selected[, SNP := reference_id]

duplicate <- duplicated(selected$SNP) | duplicated(selected$SNP, fromLast = TRUE)
duplicateCOUNT <- sum(duplicate)
selected <- selected[!duplicate]
if (nrow(selected) == 0L) {
  stop(sprintf("GWAS '%s' has no non-palindromic variants aligned to both target '%s' and the European LD reference.", option[["trait-id"]], option[["cohort"]]), call. = FALSE)
}

aligned <- selected[, ..requiredGWAS]
data.table::setorder(aligned, CHR, BP)
data.table::fwrite(aligned, option[["output-cojo"]], sep = "\t", quote = FALSE, na = "NA")
data.table::fwrite(
  aligned[, .(ID = SNP, CHR, POS = BP, A1, P = p)],
  option[["output-clump"]], sep = "\t", quote = FALSE
)
data.table::fwrite(
  data.table::data.table(
    cohort = option[["cohort"]],
    trait_id = option[["trait-id"]],
    prs_name = option[["prs-name"]],
    input_variants = nrow(gwas),
    filtered_nonpositive_p = sum(!positiveP),
    filtered_palindromic = sum(palindromic & positiveP),
    target_aligned_variants = length(targetALIGNED),
    filtered_target_missing_or_mismatched = length(eligibleROWS) - length(targetALIGNED),
    filtered_target_ambiguous = length(targetAMBIGUOUS),
    reference_aligned_variants = nrow(aligned),
    filtered_reference_missing_or_mismatched = length(setdiff(targetALIGNED, referenceALIGNED)),
    filtered_reference_ambiguous = length(referenceAMBIGUOUS),
    filtered_reference_duplicate = duplicateCOUNT,
    complemented_alleles = sum(selected$orientation %in% c("COMPLEMENT_ALT", "COMPLEMENT_REF")),
    effect_flips = sum(selected$orientation %in% c("DIRECT_REF", "COMPLEMENT_REF")),
    status = "PASS"
  ),
  option[["output-qc"]], sep = "\t", quote = FALSE
)
