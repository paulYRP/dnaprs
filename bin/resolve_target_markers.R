#!/usr/bin/env Rscript

argument <- commandArgs(trailingOnly = TRUE)
if (length(argument) %% 2L != 0L) stop("Arguments must be --name value pairs.", call. = FALSE)
option <- stats::setNames(as.list(argument[seq.int(2L, length(argument), 2L)]), sub("^--", "", argument[seq.int(1L, length(argument), 2L)]))

readTSV <- function(path, ...) utils::read.delim(
  path, sep = "\t", header = TRUE, quote = "", comment.char = "", check.names = FALSE,
  stringsAsFactors = FALSE, ...
)
readPVAR <- function(path) {
  connection <- file(path, open = "rt")
  on.exit(close(connection), add = TRUE)
  metadataLINES <- 0L
  repeat {
    line <- readLines(connection, n = 1L, warn = FALSE)
    if (length(line) == 0L) stop("PVAR has no #CHROM header.", call. = FALSE)
    if (startsWith(line, "#CHROM\t")) break
    if (!startsWith(line, "##")) {
      stop("PVAR contains content before its #CHROM header that is not VCF metadata.", call. = FALSE)
    }
    metadataLINES <- metadataLINES + 1L
  }
  readTSV(path, skip = metadataLINES)
}
writeTSV <- function(value, path, header = TRUE) utils::write.table(
  value, path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = header, na = ""
)
complement <- function(value) chartr("ACGT", "TGCA", value)
validALLELE <- function(value) grepl("^[ACGT]$", value)
samePAIR <- function(first, second) identical(sort(first), sort(second))

pvar <- readPVAR(option[["pvar"]])
names(pvar)[names(pvar) == "#CHROM"] <- "CHROM"
requiredPVAR <- c("CHROM", "POS", "ID", "REF", "ALT")
if (!all(requiredPVAR %in% names(pvar))) stop("PVAR is missing CHROM, POS, ID, REF, or ALT.", call. = FALSE)
if (anyDuplicated(pvar$ID)) {
  stop("Raw marker identifiers must be unique before dbSNP resolution; duplicate source IDs cannot be audited safely.", call. = FALSE)
}

initial <- readTSV(option[["initial-decisions"]])
if (nrow(initial) != nrow(pvar)) stop("Initial marker decisions and PVAR row counts differ.", call. = FALSE)
if (!identical(as.character(initial$final_id), as.character(pvar$ID))) {
  stop("Initial marker decisions are not in the same order as the annotated PVAR.", call. = FALSE)
}

calls <- readTSV(option[["calls"]])
names(calls)[names(calls) == "#CHROM"] <- "CHROM"
if (!all(c("SNP", "COUNTED", "ALT") %in% names(calls)) || nrow(calls) != nrow(pvar)) {
  stop("PLINK A-transpose calls do not correspond one-to-one with the annotated PVAR.", call. = FALSE)
}
if (!identical(as.character(calls$SNP), as.character(pvar$ID))) {
  stop("PLINK A-transpose marker order differs from the annotated PVAR.", call. = FALSE)
}
sampleCOLUMN <- setdiff(names(calls), c("CHROM", "SNP", "(C)M", "CM", "POS", "COUNTED", "ALT"))
if (length(sampleCOLUMN) == 0L) stop("PLINK A-transpose output contains no participant dosages.", call. = FALSE)
dosage <- as.matrix(data.frame(lapply(calls[sampleCOLUMN], as.numeric), check.names = FALSE))

chromosomeMAP <- utils::read.delim(
  option[["chromosome-map"]], sep = "\t", header = FALSE, quote = "", comment.char = "",
  stringsAsFactors = FALSE, col.names = c("chromosome", "accession")
)
dbsnp <- utils::read.delim(
  option[["dbsnp-records"]], sep = "\t", header = FALSE, quote = "", comment.char = "",
  stringsAsFactors = FALSE, col.names = c("accession", "position", "candidate", "reference", "alternate")
)
dbsnp$chromosome <- chromosomeMAP$chromosome[match(dbsnp$accession, chromosomeMAP$accession)]
dbsnp <- dbsnp[
  !is.na(dbsnp$chromosome) & grepl("^rs[0-9]+$", dbsnp$candidate) &
    validALLELE(toupper(dbsnp$reference)),
  , drop = FALSE
]

# Expand multiallelic dbSNP records so each compatible biallelic pair is assessed
# independently while retaining the candidate rsID.
expanded <- list()
for (row in seq_len(nrow(dbsnp))) {
  alternate <- strsplit(toupper(dbsnp$alternate[[row]]), ",", fixed = TRUE)[[1L]]
  alternate <- alternate[validALLELE(alternate)]
  for (allele in alternate) {
    expanded[[length(expanded) + 1L]] <- data.frame(
      chromosome = as.character(dbsnp$chromosome[[row]]),
      position = as.integer(dbsnp$position[[row]]),
      candidate = dbsnp$candidate[[row]],
      candidate_ref = toupper(dbsnp$reference[[row]]),
      candidate_alt = allele,
      stringsAsFactors = FALSE
    )
  }
}
candidate <- if (length(expanded) > 0L) unique(do.call(rbind, expanded)) else data.frame(
  chromosome = character(), position = integer(), candidate = character(),
  candidate_ref = character(), candidate_alt = character(), stringsAsFactors = FALSE
)

result <- data.frame(
  source_id = initial$source_id,
  final_id = "",
  source_chr = initial$source_chr,
  source_pos = initial$source_pos,
  final_chr = pvar$CHROM,
  final_pos = pvar$POS,
  source_ref = initial$source_ref,
  source_alt = initial$source_alt,
  final_ref = pvar$REF,
  final_alt = pvar$ALT,
  coordinate_candidates = integer(nrow(pvar)),
  allele_compatible_candidates = integer(nrow(pvar)),
  call_count = integer(nrow(pvar)),
  call_rate = numeric(nrow(pvar)),
  probe_count = integer(nrow(pvar)),
  assay_count = integer(nrow(pvar)),
  overlap_count = integer(nrow(pvar)),
  concordance = rep(NA_real_, nrow(pvar)),
  decision = "EXCLUDED_UNRESOLVED_MARKER",
  reason = "No allele-compatible dbSNP match",
  stringsAsFactors = FALSE
)
matchTYPE <- rep("", nrow(pvar))
canonicalDOSAGE <- matrix(NA_real_, nrow = nrow(pvar), ncol = ncol(dosage))

for (row in seq_len(nrow(pvar))) {
  chromosome <- sub("^chr", "", as.character(pvar$CHROM[[row]]), ignore.case = TRUE)
  position <- suppressWarnings(as.integer(pvar$POS[[row]]))
  sourcePAIR <- toupper(c(pvar$REF[[row]], pvar$ALT[[row]]))
  result$call_count[[row]] <- sum(is.finite(dosage[row, ]))
  result$call_rate[[row]] <- result$call_count[[row]] / ncol(dosage)
  if (!all(validALLELE(sourcePAIR))) {
    result$reason[[row]] <- if (all(sourcePAIR %in% c("0", "."))) "All genotypes or assay alleles missing" else "Incompatible assay alleles"
    next
  }
  atCOORDINATE <- candidate[candidate$chromosome == chromosome & candidate$position == position, , drop = FALSE]
  result$coordinate_candidates[[row]] <- length(unique(atCOORDINATE$candidate))
  if (nrow(atCOORDINATE) == 0L) {
    result$reason[[row]] <- "No dbSNP coordinate match"
    next
  }
  compatible <- vapply(seq_len(nrow(atCOORDINATE)), function(index) {
    pair <- c(atCOORDINATE$candidate_ref[[index]], atCOORDINATE$candidate_alt[[index]])
    samePAIR(sourcePAIR, pair) || samePAIR(complement(sourcePAIR), pair)
  }, logical(1L))
  atCOORDINATE <- atCOORDINATE[compatible, , drop = FALSE]
  compatibleID <- unique(atCOORDINATE$candidate)
  result$allele_compatible_candidates[[row]] <- length(compatibleID)
  if (length(compatibleID) != 1L) {
    result$reason[[row]] <- if (length(compatibleID) == 0L) "No allele-compatible dbSNP match" else "Multiple allele-compatible dbSNP matches"
    next
  }
  selected <- atCOORDINATE[atCOORDINATE$candidate == compatibleID[[1L]], , drop = FALSE]
  selected$direct <- vapply(seq_len(nrow(selected)), function(index) {
    samePAIR(sourcePAIR, c(selected$candidate_ref[[index]], selected$candidate_alt[[index]]))
  }, logical(1L))
  selected <- selected[order(!selected$direct, selected$candidate_ref, selected$candidate_alt), , drop = FALSE][1L, , drop = FALSE]
  result$final_id[[row]] <- selected$candidate[[1L]]
  result$final_chr[[row]] <- chromosome
  result$final_pos[[row]] <- position
  result$final_ref[[row]] <- selected$candidate_ref[[1L]]
  result$final_alt[[row]] <- selected$candidate_alt[[1L]]
  matchTYPE[[row]] <- if (selected$direct[[1L]]) "DIRECT" else "COMPLEMENT"

  counted <- toupper(as.character(calls$COUNTED[[row]]))
  if (matchTYPE[[row]] == "COMPLEMENT") counted <- complement(counted)
  anchor <- min(result$final_ref[[row]], result$final_alt[[row]])
  other <- max(result$final_ref[[row]], result$final_alt[[row]])
  if (counted == anchor) canonicalDOSAGE[row, ] <- dosage[row, ]
  if (counted == other) canonicalDOSAGE[row, ] <- 2 - dosage[row, ]
  if (!counted %in% c(anchor, other)) {
    result$final_id[[row]] <- ""
    result$reason[[row]] <- "PLINK counted allele is incompatible with selected dbSNP alleles"
  } else if (result$call_count[[row]] == 0L) {
    result$final_id[[row]] <- ""
    result$reason[[row]] <- "All genotypes missing"
  } else {
    result$reason[[row]] <- "Unique allele-compatible dbSNP match"
  }
}

eligible <- which(nzchar(result$final_id))
for (rsid in unique(result$final_id[eligible])) {
  index <- eligible[result$final_id[eligible] == rsid]
  assayKEY <- paste(pmin(result$final_ref[index], result$final_alt[index]), pmax(result$final_ref[index], result$final_alt[index]), sep = "/")
  result$probe_count[index] <- length(index)
  result$assay_count[index] <- length(unique(assayKEY))
  if (length(index) == 1L) {
    result$decision[index] <- "RETAINED_UNIQUE"
    next
  }
  if (length(unique(assayKEY)) != 1L) {
    result$decision[index] <- "EXCLUDED_DIFFERENT_ASSAY_DUPLICATE_GROUP"
    result$reason[index] <- "Repeated rsID has more than one reference-oriented assay pair"
    next
  }
  exactID <- result$source_id[index] == rsid
  representative <- index[order(-result$call_rate[index], -as.integer(exactID), index)][[1L]]
  groupCONCORDANT <- TRUE
  for (row in index) {
    observed <- is.finite(canonicalDOSAGE[representative, ]) & is.finite(canonicalDOSAGE[row, ])
    result$overlap_count[[row]] <- sum(observed)
    result$concordance[[row]] <- if (any(observed)) mean(abs(canonicalDOSAGE[representative, observed] - canonicalDOSAGE[row, observed]) < 1e-8) else NA_real_
    if (!is.finite(result$concordance[[row]]) || result$concordance[[row]] != 1) groupCONCORDANT <- FALSE
  }
  if (!groupCONCORDANT) {
    result$decision[index] <- "EXCLUDED_DISCORDANT_DUPLICATE_GROUP"
    result$reason[index] <- "Same-assay duplicate probes are not completely concordant over observed calls"
  } else {
    result$decision[index] <- "EXCLUDED_REDUNDANT_DUPLICATE_PROBE"
    result$reason[index] <- "Completely concordant duplicate represented by higher-call-rate probe"
    result$decision[representative] <- "RETAINED_DUPLICATE_REPRESENTATIVE"
    result$reason[representative] <- "Highest call rate, then exact rsID, then original row order"
  }
}

retained <- result$decision %in% c("RETAINED_UNIQUE", "RETAINED_DUPLICATE_REPRESENTATIVE")
if (!any(retained)) stop("Raw marker resolution retained no variants.", call. = FALSE)
if (anyDuplicated(result$final_id[retained])) stop("Raw marker resolution did not produce unique final rsIDs.", call. = FALSE)

writeTSV(result, option[["output-decisions"]])
writeTSV(data.frame(id = result$final_id[retained]), option[["keep"]], header = FALSE)
writeTSV(data.frame(source_id = pvar$ID[retained], final_id = result$final_id[retained]), option[["rename"]], header = FALSE)
