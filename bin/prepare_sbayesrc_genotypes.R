#!/usr/bin/env Rscript

if (!requireNamespace("data.table", quietly = TRUE)) stop("The data.table package is required.", call. = FALSE)
argument <- commandArgs(trailingOnly = TRUE)
option <- stats::setNames(as.list(argument[seq.int(2L, length(argument), 2L)]), sub("^--", "", argument[seq.int(1L, length(argument), 2L)]))
data.table::setDTthreads(as.integer(if (is.null(option[["threads"]])) 1L else option[["threads"]]))
cohort <- option[["cohort"]]
fail <- function(message) stop(sprintf("SBayesRC genotype preparation for cohort '%s': %s", cohort, message), call. = FALSE)
readTABLE <- function(path, ...) data.table::fread(path, colClasses = "character", showProgress = FALSE, ...)
sampleKEY <- function(sample) paste(sample[["FID"]], sample[["IID"]], sep = "\t")
readSAMPLES <- function(prefix) {
  sample <- readTABLE(paste0(prefix, ".psam"))
  data.table::setnames(sample, sub("^#", "", names(sample)))
  if (!all(c("FID", "IID") %in% names(sample))) fail(paste("Missing FID/IID columns in", paste0(prefix, ".psam")))
  sample[, c("FID", "IID"), with = FALSE]
}
checkSAMPLES <- function(sample, keep, chromosome) {
  if (nrow(sample) != nrow(keep) || anyDuplicated(sampleKEY(sample)) || !setequal(sampleKEY(sample), sampleKEY(keep))) {
    fail(sprintf("Chromosome %s does not contain exactly the eligible participants. Review the keep file and PSAM.", chromosome))
  }
}
writeTABLE <- function(table, path) data.table::fwrite(table, path, sep = "\t", quote = FALSE)
keep <- readTABLE(option[["keep"]], header = FALSE)
if (ncol(keep) != 2L || nrow(keep) == 0L || anyNA(keep)) fail("The participant keep file must contain non-empty FID/IID pairs.")
data.table::setnames(keep, c("FID", "IID"))
if (anyDuplicated(sampleKEY(keep))) fail("The participant keep file contains duplicate FID/IID pairs.")

if (option[["action"]] == "prepare") {
  chromosome <- as.integer(option[["chromosome"]])
  if (length(chromosome) != 1L || is.na(chromosome) || !chromosome %in% 1:22) fail("Chromosome must be an integer from 1 to 22.")
  vcf <- option[["vcf"]]
  outputDIR <- paste0(cohort, ".chr", chromosome, ".sbayesrc_genotypes")
  if (dir.exists(outputDIR)) fail(paste("Output directory already exists:", outputDIR))
  dir.create(outputDIR)
  prefix <- file.path(outputDIR, paste0(cohort, "_chr", chromosome))
  # Read variant metadata only; gzip is available in the scoring container.
  source <- data.table::fread(
    cmd = sprintf("gzip -cd -- %s | cut -f1-5", shQuote(vcf)),
    skip = "#CHROM", select = c("#CHROM", "POS", "ID", "REF", "ALT"),
    colClasses = "character", showProgress = FALSE
  )
  sourceCHR <- suppressWarnings(as.integer(sub("^chr", "", source[["#CHROM"]])))
  if (nrow(source) == 0L || anyNA(source) || anyNA(sourceCHR) || any(sourceCHR != chromosome)) {
    fail(sprintf("VCF '%s' must contain non-empty chromosome %s records.", vcf, chromosome))
  }
  expectedID <- source$ID
  missingID <- expectedID == "."
  expectedID[missingID] <- paste(chromosome, source$POS[missingID], source$REF[missingID], source$ALT[missingID], sep = ":")
  duplicatedID <- unique(expectedID[duplicated(expectedID)])
  if (length(duplicatedID)) {
    fail(sprintf("Chromosome %s has duplicate variant IDs: %s. Review the VCF; no records were removed.", chromosome, paste(head(duplicatedID, 5L), collapse = ", ")))
  }
  status <- system2("plink2", c(
    "--vcf", shQuote(vcf), "dosage=DS", "--double-id", "--keep", shQuote(option[["keep"]]),
    "--set-missing-var-ids", shQuote("@:#:$r:$a"), "--make-pgen",
    "--threads", option[["threads"]], "--memory", option[["memory"]], "--out", shQuote(prefix)
  ))
  if (status != 0L) fail(sprintf("PLINK conversion failed for chromosome %s. Review %s.log.", chromosome, prefix))
  pvar <- readTABLE(paste0(prefix, ".pvar"), skip = "#CHROM", select = c("#CHROM", "POS", "ID", "REF", "ALT"))
  sample <- readSAMPLES(prefix)
  checkSAMPLES(sample, keep, chromosome)
  if (nrow(pvar) != nrow(source) || !identical(pvar$ID, expectedID) ||
      !identical(pvar$POS, source$POS) || !identical(pvar$REF, source$REF) || !identical(pvar$ALT, source$ALT) ||
      any(pvar[["#CHROM"]] != as.character(chromosome))) {
    fail(sprintf("Chromosome %s conversion changed variant records, alleles or expected IDs. Review %s.", chromosome, prefix))
  }
  writeTABLE(data.table::data.table(
    cohort = cohort, chromosome = chromosome, participants = nrow(sample), variants = nrow(pvar),
    preserved_ids = sum(!missingID), assigned_ids = sum(missingID), duplicate_ids = 0L,
    participant_match = TRUE, status = "PASS"
  ), paste0(cohort, ".chr", chromosome, ".sbayesrc.genotype_qc.tsv"))
} else if (option[["action"]] == "assemble") {
  manifest <- data.table::fread(option[["manifest"]], colClasses = list(integer = "chromosome"), showProgress = FALSE)
  if (!all(c("chromosome", "directory", "qc") %in% names(manifest)) || nrow(manifest) != 22L ||
      anyDuplicated(manifest$chromosome) || !setequal(manifest$chromosome, 1:22)) {
    fail("Expected exactly one prepared bundle for each chromosome 1-22. Review missing or duplicate chromosome inputs.")
  }
  data.table::setorder(manifest, chromosome)
  outputDIR <- paste0(cohort, ".sbayesrc_genotypes")
  if (dir.exists(outputDIR)) fail(paste("Output directory already exists:", outputDIR))
  dir.create(outputDIR)
  firstKEY <- NULL
  summaries <- vector("list", nrow(manifest))
  for (i in seq_len(nrow(manifest))) {
    chromosome <- manifest$chromosome[i]
    prefix <- file.path(manifest$directory[i], paste0(cohort, "_chr", chromosome))
    companions <- paste0(prefix, c(".pgen", ".pvar", ".psam"))
    if (any(!file.exists(companions)) || any(file.info(companions)$size == 0)) fail(paste("Incomplete chromosome bundle:", prefix))
    sample <- readSAMPLES(prefix)
    checkSAMPLES(sample, keep, chromosome)
    if (is.null(firstKEY)) firstKEY <- sampleKEY(sample)
    if (!identical(sampleKEY(sample), firstKEY)) fail(sprintf("Chromosome %s has a different participant order.", chromosome))
    summaries[[i]] <- data.table::fread(manifest$qc[i], showProgress = FALSE)
    if (nrow(summaries[[i]]) != 1L || summaries[[i]]$cohort != cohort || summaries[[i]]$chromosome != chromosome || summaries[[i]]$status != "PASS") {
      fail(paste("QC does not match the prepared chromosome:", prefix))
    }
    if (!all(file.copy(companions, outputDIR, overwrite = FALSE))) fail(paste("Could not collect prepared chromosome:", prefix))
  }
  writeTABLE(data.table::rbindlist(summaries), paste0(cohort, ".sbayesrc.genotype_qc.tsv"))
} else {
  fail(paste("Unknown action:", option[["action"]]))
}
