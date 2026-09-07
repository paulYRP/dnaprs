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
readPVAR <- function(path, prefix) {
  value <- data.table::fread(path, skip = "#CHROM", select = c("#CHROM", "POS", "ID", "REF", "ALT"), colClasses = "character")
  data.table::setnames(value, "#CHROM", "CHR", skip_absent = TRUE)
  required <- c("CHR", "POS", "ID", "REF", "ALT")
  if (!all(required %in% names(value))) stop(sprintf("%s PVAR is missing CHR, POS, ID, REF, or ALT.", prefix), call. = FALSE)
  value <- value[
    ID != "" & ID != "." & grepl("^[ACGT]$", toupper(REF)) & grepl("^[ACGT]$", toupper(ALT)) & toupper(REF) != toupper(ALT),
    .(
      CHR = sub("^chr", "", as.character(CHR), ignore.case = TRUE),
      BP = as.integer(POS),
      ID = as.character(ID),
      REF = toupper(REF),
      ALT = toupper(ALT)
    )
  ]
  value[, canonical_pair := canonicalPAIR(REF, ALT)]
  value[, canonical_key := paste(CHR, BP, canonical_pair, sep = ":")]
  value
}

data.table::setDTthreads(as.integer(option[["threads"]]))

target <- readPVAR(option[["target-pvar"]], "Target")
reference <- readPVAR(option[["reference-pvar"]], "PLINK reference")
if (anyDuplicated(target$ID)) stop("Target PVAR variant IDs must be unique before PLINK alignment.", call. = FALSE)
if (anyDuplicated(reference$ID)) stop("PLINK reference PVAR variant IDs must be unique.", call. = FALSE)
target[, target_key_count := .N, by = canonical_key]
reference <- reference[canonical_key %in% target$canonical_key]
reference[, reference_key_count := .N, by = canonical_key]
index <- merge(
  target[!duplicated(canonical_key), .(canonical_key, CHR, BP, canonical_pair, target_id = ID, target_ref = REF, target_alt = ALT, target_key_count)],
  reference[!duplicated(canonical_key), .(canonical_key, reference_id = ID, reference_ref = REF, reference_alt = ALT, reference_key_count)],
  by = "canonical_key", all.x = TRUE, sort = FALSE
)
index[is.na(reference_key_count), reference_key_count := 0L]
index[, unique_compatible := target_key_count == 1L & reference_key_count == 1L &
  target_ref == reference_ref & target_alt == reference_alt]
index[is.na(unique_compatible), unique_compatible := FALSE]
data.table::setorder(index, CHR, BP, canonical_pair)
data.table::fwrite(index, option[["output-index"]], sep = "\t", quote = FALSE, na = "")
