#!/usr/bin/env Rscript

argument <- commandArgs(trailingOnly = TRUE)
option <- stats::setNames(as.list(argument[seq.int(2L, length(argument), 2L)]), sub("^--", "", argument[seq.int(1L, length(argument), 2L)]))
if (!requireNamespace("data.table", quietly = TRUE)) stop("The data.table package is required.", call. = FALSE)

weight <- data.table::fread(option[["weights"]], colClasses = "character")
pvar <- data.table::fread(option[["pvar"]], colClasses = "character")
used <- unique(data.table::fread(option[["used"]], header = FALSE, colClasses = "character")[[1L]])
names(pvar)[names(pvar) == "#CHROM"] <- "CHROM"
if (!all(c("SNP", "A1") %in% names(weight)) || !all(c("ID", "REF", "ALT") %in% names(pvar))) {
  stop("Scoring compatibility requires SNP/A1 weights and ID/REF/ALT PVAR columns.", call. = FALSE)
}
if (anyDuplicated(weight$SNP) || anyDuplicated(pvar$ID) || anyDuplicated(used)) {
  stop("Weights, target IDs, and used-variant IDs must each be unique.", call. = FALSE)
}
complement <- function(value) chartr("ACGT", "TGCA", toupper(value))
target <- pvar[match(weight$SNP, pvar$ID)]
effect <- toupper(weight$A1)
ref <- toupper(target$REF)
alt <- toupper(target$ALT)
state <- ifelse(
  is.na(target$ID), "MISSING_TARGET",
  ifelse(effect == ref | effect == alt, "DIRECT", ifelse(complement(effect) == ref | complement(effect) == alt, "COMPLEMENT", "INCOMPATIBLE"))
)
audit <- data.table::data.table(
  cohort = option[["cohort"]], trait_id = option[["trait-id"]], method = option[["method"]],
  scoring_stage = option[["scoring-stage"]], variant_id = weight$SNP, effect_allele = effect,
  target_ref = ref, target_alt = alt, allele_state = state, used = weight$SNP %in% used
)
audit[, reason := data.table::fcase(
  used, "Used by PLINK scoring",
  state == "MISSING_TARGET", "Weight variant is absent from target",
  state == "INCOMPATIBLE", "Effect allele is incompatible with target alleles",
  default = "Compatible variant was not reported in the PLINK used-variant set"
)]
if (any(audit$used & audit$allele_state %in% c("MISSING_TARGET", "INCOMPATIBLE"))) {
  stop("PLINK reported an absent or allele-incompatible weight as used.", call. = FALSE)
}
requested <- nrow(audit)
usedN <- sum(audit$used)
if (usedN == 0L) stop("No requested weight variant was used for scoring.", call. = FALSE)
summary <- data.table::data.table(
  cohort = option[["cohort"]], trait_id = option[["trait-id"]], method = option[["method"]],
  scoring_stage = option[["scoring-stage"]], requested_variants = requested,
  target_compatible_variants = sum(audit$allele_state %in% c("DIRECT", "COMPLEMENT")),
  used_variants = usedN, used_fraction = usedN / requested,
  review_required = usedN < requested,
  status = if (usedN < requested) "REVIEW" else "PASS"
)
prefix <- paste(option[["cohort"]], option[["trait-id"]], option[["method"]], sep = ".")
data.table::fwrite(audit, paste0(prefix, ".variant_compatibility.tsv"), sep = "\t", na = "")
data.table::fwrite(summary, paste0(prefix, ".variant_coverage.tsv"), sep = "\t")
