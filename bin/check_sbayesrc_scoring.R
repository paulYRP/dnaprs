#!/usr/bin/env Rscript

if (!requireNamespace("data.table", quietly = TRUE)) stop("The data.table package is required.", call. = FALSE)
argument <- commandArgs(trailingOnly = TRUE)
option <- stats::setNames(as.list(argument[seq.int(2L, length(argument), 2L)]), sub("^--", "", argument[seq.int(1L, length(argument), 2L)]))
data.table::setDTthreads(1L)
cohort <- option[["cohort"]]
trait <- option[["trait-id"]]
fail <- function(message) stop(sprintf("SBayesRC scoring for cohort '%s', trait '%s': %s", cohort, trait, message), call. = FALSE)
weightPATH <- option[["weight"]]
weight <- data.table::fread(weightPATH, select = c("SNP", "A1", "BETA"), colClasses = c(SNP = "character", A1 = "character"), showProgress = FALSE)
if (nrow(weight) == 0L || anyNA(weight) || any(!is.finite(weight$BETA)) || anyDuplicated(weight$SNP)) {
  fail(paste("Weights must have unique, non-missing SNP IDs, effect alleles and finite effects:", weightPATH))
}
keep <- data.table::fread(option[["keep"]], header = FALSE, colClasses = "character", showProgress = FALSE)
if (ncol(keep) != 2L || nrow(keep) == 0L || anyNA(keep)) fail("The participant keep file must contain FID/IID pairs.")
keepKEY <- paste(keep[[1L]], keep[[2L]], sep = "\t")
if (anyDuplicated(keepKEY)) fail("The participant keep file contains duplicate pairs.")
matched <- rep(FALSE, nrow(weight))
summaries <- vector("list", 22L)
firstKEY <- NULL
for (chromosome in 1:22) {
  prefix <- file.path(option[["target-dir"]], paste0(cohort, "_chr", chromosome))
  companions <- paste0(prefix, c(".pgen", ".pvar", ".psam"))
  if (any(!file.exists(companions)) || any(file.info(companions)$size == 0)) {
    fail(sprintf("Missing chromosome %s scoring files at '%s'. Prepare a complete chromosome 1-22 genotype set.", chromosome, prefix))
  }
  sample <- data.table::fread(paste0(prefix, ".psam"), colClasses = "character", showProgress = FALSE)
  data.table::setnames(sample, sub("^#", "", names(sample)))
  if (!all(c("FID", "IID") %in% names(sample))) fail(paste("Missing FID/IID columns in", paste0(prefix, ".psam")))
  sampleKEY <- paste(sample$FID, sample$IID, sep = "\t")
  if (anyDuplicated(sampleKEY) || any(!keepKEY %in% sampleKEY)) fail(sprintf("Chromosome %s has duplicate participants or is missing eligible participants.", chromosome))
  selectedKEY <- sampleKEY[sampleKEY %in% keepKEY]
  if (is.null(firstKEY)) firstKEY <- selectedKEY
  if (!identical(firstKEY, selectedKEY)) fail(sprintf("Chromosome %s has a different eligible-participant order.", chromosome))
  variant <- data.table::fread(paste0(prefix, ".pvar"), skip = "#CHROM", select = c("#CHROM", "ID", "REF", "ALT"), colClasses = "character", showProgress = FALSE)
  if (nrow(variant) == 0L || anyNA(variant) || anyDuplicated(variant$ID) || any(variant[["#CHROM"]] != as.character(chromosome))) {
    fail(sprintf("Chromosome %s has empty, duplicate or inconsistent variant records in '%s.pvar'.", chromosome, prefix))
  }
  # Match each genotype ID once; do not search the complete weights table per variant.
  index <- match(variant$ID, weight$SNP)
  found <- !is.na(index)
  compatible <- found & (weight$A1[index] == variant$REF | weight$A1[index] == variant$ALT)
  matched[index[found]] <- TRUE
  summaries[[chromosome]] <- data.table::data.table(
    cohort = cohort, trait_id = trait, chromosome = as.character(chromosome),
    target_variants = nrow(variant), matched_ids = sum(found),
    incompatible_alleles = sum(found & !compatible), usable_variants = sum(compatible),
    status = if (sum(compatible) == 0L) "FAIL" else if (any(found & !compatible)) "REVIEW" else "PASS"
  )
}
summary <- data.table::rbindlist(summaries)
summary[, requested_weights := nrow(weight)]
summary[, unmatched_weights := sum(!matched)]
output <- paste0(cohort, ".", trait, ".sbayesrc.match_qc.tsv")
data.table::fwrite(summary, output, sep = "\t", quote = FALSE)
failed <- summary[usable_variants == 0L]
if (nrow(failed)) {
  chromosome <- failed$chromosome[1L]
  fail(sprintf(
    "Chromosome %s has no usable variants: %s matched IDs and %s incompatible alleles. Weights: '%s'; genotypes: '%s/%s_chr%s.pvar'. Review variant IDs and effect alleles; no chromosome was skipped. Match counts: %s.",
    chromosome, failed$matched_ids[1L], failed$incompatible_alleles[1L], weightPATH, option[["target-dir"]], cohort, chromosome, output
  ))
}
message(sprintf("SBayesRC match check: %s/%s weights have genotype IDs; %s usable chromosome-variant records. See %s.", sum(matched), nrow(weight), sum(summary$usable_variants), output))
