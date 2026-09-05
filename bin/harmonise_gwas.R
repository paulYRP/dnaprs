#!/usr/bin/env Rscript

# Parse the declared GWAS column mapping.
argument <- commandArgs(trailingOnly = TRUE)
name <- sub("^--", "", argument[seq.int(1L, length(argument), 2L)])
value <- argument[seq.int(2L, length(argument), 2L)]
option <- stats::setNames(as.list(value), name)

if (!requireNamespace("data.table", quietly = TRUE)) {
  stop("The data.table package is required.", call. = FALSE)
}

optionalOPTION <- c(
  "freq-col", "case-freq-col", "control-freq-col", "case-n-col", "control-n-col",
  "n-col", "info-col", "info-min", "maf-min", "source-format"
)
for (item in optionalOPTION) if (is.null(option[[item]])) option[[item]] <- ""
if (option[["maf-min"]] == "") option[["maf-min"]] <- "0.01"
if (option[["source-format"]] == "") option[["source-format"]] <- "auto"

requiredOPTION <- c(
  "input", "trait-id", "prs-name", "effect-type", "sample-size", "snp-col", "chr-col",
  "bp-col", "effect-allele-col", "other-allele-col", "beta-col", "se-col", "p-col"
)
missingOPTION <- setdiff(requiredOPTION, names(option))
if (length(missingOPTION) > 0L) {
  stop(sprintf("Missing option(s): %s", paste(missingOPTION, collapse = ", ")), call. = FALSE)
}

compressedINPUT <- grepl("\\.(gz|bgz)$", option[["input"]], ignore.case = TRUE)
if (compressedINPUT && !requireNamespace("R.utils", quietly = TRUE)) {
  stop(
    sprintf(
      "Compressed GWAS input '%s' requires the R.utils package. Use the pipeline analysis container.",
      basename(option[["input"]])
    ),
    call. = FALSE
  )
}

# Locate the real table header after consortium metadata, then read only the columns
# explicitly declared in the manifest. This supports metadata-prefixed TSV/VCF-style
# releases and whitespace-delimited tables without study-specific code in the scorer.
sourceCOLUMN <- unique(c(
  option[["snp-col"]], option[["chr-col"]], option[["bp-col"]],
  option[["effect-allele-col"]], option[["other-allele-col"]], option[["beta-col"]],
  option[["se-col"]], option[["p-col"]], option[["freq-col"]], option[["n-col"]],
  option[["case-freq-col"]], option[["control-freq-col"]], option[["case-n-col"]],
  option[["control-n-col"]], option[["info-col"]]
))
sourceCOLUMN <- sourceCOLUMN[!is.na(sourceCOLUMN) & sourceCOLUMN != ""]
connection <- if (compressedINPUT) gzfile(option[["input"]], "rt") else file(option[["input"]], "rt")
headerCANDIDATE <- readLines(connection, n = 500L, warn = FALSE)
close(connection)
splitHEADER <- function(value) {
  if (grepl("\t", value, fixed = TRUE)) strsplit(value, "\t", fixed = TRUE)[[1L]] else strsplit(trimws(value), "[[:space:]]+")[[1L]]
}
headerLIST <- lapply(headerCANDIDATE, splitHEADER)
normalisedHEADER <- lapply(headerLIST, function(value) sub("^#", "", value))
headerINDEX <- which(vapply(normalisedHEADER, function(value) all(sourceCOLUMN %in% value), logical(1L)))[1L]
if (is.na(headerINDEX)) stop("Could not locate a GWAS header containing every declared source column.", call. = FALSE)
actualHEADER <- headerLIST[[headerINDEX]]
normalHEADER <- normalisedHEADER[[headerINDEX]]
actualCOLUMN <- actualHEADER[match(sourceCOLUMN, normalHEADER)]
gwas <- data.table::fread(
  option[["input"]],
  skip = headerINDEX - 1L,
  select = actualCOLUMN,
  showProgress = FALSE,
  check.names = FALSE
)
data.table::setnames(gwas, names(gwas), sub("^#", "", names(gwas)))

extractCOLUMN <- function(column, default = NA) {
  if (is.null(column) || is.na(column) || column == "") rep(default, nrow(gwas)) else gwas[[column]]
}

# Some case-control releases provide frequencies separately for cases and controls.
# Derive the combined effect-allele frequency with the row-specific sample counts so
# the downstream SBayesRC input has the same meaning as a directly supplied EAF.
directFREQUENCY <- suppressWarnings(as.numeric(extractCOLUMN(option[["freq-col"]])))
caseFREQUENCY <- suppressWarnings(as.numeric(extractCOLUMN(option[["case-freq-col"]])))
controlFREQUENCY <- suppressWarnings(as.numeric(extractCOLUMN(option[["control-freq-col"]])))
caseN <- suppressWarnings(as.numeric(extractCOLUMN(option[["case-n-col"]])))
controlN <- suppressWarnings(as.numeric(extractCOLUMN(option[["control-n-col"]])))
weightedFREQUENCY <- (caseFREQUENCY * caseN + controlFREQUENCY * controlN) / (caseN + controlN)
frequency <- ifelse(is.finite(directFREQUENCY), directFREQUENCY, weightedFREQUENCY)

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
  freq = frequency,
  b = effect,
  se = suppressWarnings(as.numeric(extractCOLUMN(option[["se-col"]]))),
  p = suppressWarnings(as.numeric(extractCOLUMN(option[["p-col"]]))),
  N = suppressWarnings(as.numeric(extractCOLUMN(option[["n-col"]], option[["sample-size"]])))
)
information <- suppressWarnings(as.numeric(extractCOLUMN(option[["info-col"]])))

# Apply one common PRS-specific QC contract after format adaptation.
validALLELE <- grepl("^[ACGT]+$", cojo$A1) & grepl("^[ACGT]+$", cojo$A2) & cojo$A1 != cojo$A2
structuralROW <- !is.na(cojo$CHR) & cojo$CHR %in% as.character(1:22) &
  is.finite(cojo$BP) & cojo$BP > 0 & validALLELE & is.finite(cojo$b) &
  is.finite(cojo$se) & cojo$se > 0 & is.finite(cojo$p) & cojo$p > 0 & cojo$p <= 1 &
  is.finite(cojo$N) & cojo$N > 0
missingID <- is.na(cojo$SNP) | cojo$SNP == "" | cojo$SNP == "."
cojo$SNP[missingID & structuralROW] <- paste(cojo$CHR[missingID & structuralROW], cojo$BP[missingID & structuralROW], cojo$A1[missingID & structuralROW], cojo$A2[missingID & structuralROW], sep = ":")
frequencyROW <- is.na(cojo$freq) | (is.finite(cojo$freq) & cojo$freq > 0 & cojo$freq < 1)
mafMIN <- suppressWarnings(as.numeric(option[["maf-min"]]))
if (!is.finite(mafMIN)) mafMIN <- 0.01
mafROW <- is.na(cojo$freq) | pmin(cojo$freq, 1 - cojo$freq) >= mafMIN
infoMIN <- suppressWarnings(as.numeric(option[["info-min"]]))
infoROW <- if (is.finite(infoMIN)) is.finite(information) & information >= infoMIN else rep(TRUE, nrow(cojo))
palindromic <- paste0(cojo$A1, cojo$A2) %in% c("AT", "TA", "CG", "GC")
ambiguousROW <- !(palindromic & is.finite(cojo$freq) & cojo$freq >= 0.4 & cojo$freq <= 0.6)
retainROW <- structuralROW & frequencyROW & mafROW & infoROW & ambiguousROW
filteredSTRUCTURAL <- sum(!structuralROW)
filteredFREQUENCY <- sum(structuralROW & !frequencyROW)
filteredMAF <- sum(structuralROW & frequencyROW & !mafROW)
filteredINFO <- sum(structuralROW & frequencyROW & mafROW & !infoROW)
filteredAMBIGUOUS <- sum(structuralROW & frequencyROW & mafROW & infoROW & !ambiguousROW)
cojo <- cojo[retainROW]
if (nrow(cojo) == 0L) stop(sprintf("GWAS '%s' retained no variants after common QC.", option[["trait-id"]]), call. = FALSE)
duplicateROW <- duplicated(cojo$SNP) | duplicated(cojo$SNP, fromLast = TRUE)
filteredDUPLICATE <- sum(duplicateROW)
cojo <- cojo[!duplicateROW]
if (nrow(cojo) == 0L) stop(sprintf("GWAS '%s' retained no unique variants.", option[["trait-id"]]), call. = FALSE)

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
    source_format = option[["source-format"]],
    source_variants = nrow(gwas),
    harmonised_variants = nrow(cojo),
    filtered_structural = filteredSTRUCTURAL,
    filtered_frequency = filteredFREQUENCY,
    filtered_maf = filteredMAF,
    filtered_info = filteredINFO,
    filtered_ambiguous = filteredAMBIGUOUS,
    filtered_duplicate = filteredDUPLICATE,
    maf_min = mafMIN,
    info_min = if (is.finite(infoMIN)) infoMIN else NA_real_,
    duplicated_snp = filteredDUPLICATE,
    structural_status = "PASS"
  ),
  paste0(option[["trait-id"]], ".harmonisation_qc.tsv"),
  sep = "\t"
)
