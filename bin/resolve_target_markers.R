#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))
argument <- commandArgs(trailingOnly = TRUE)
if (length(argument) %% 2L != 0L) stop("Arguments must be --name value pairs.", call. = FALSE)
option <- setNames(as.list(argument[seq.int(2L, length(argument), 2L)]), sub("^--", "", argument[seq.int(1L, length(argument), 2L)]))
setDTthreads(as.integer(option[["threads"]]))
readTSV <- function(path, ...) fread(path, sep = "\t", quote = "", na.strings = c("", "NA"), ...)
writeTSV <- function(value, path, header = TRUE) fwrite(value, path, sep = "\t", quote = FALSE, col.names = header, na = "")
complement <- function(value) chartr("ACGT", "TGCA", value)
pairKEY <- function(ref, alt) paste(pmin(ref, alt), pmax(ref, alt), sep = "/")
validPAIR <- function(ref, alt) !is.na(ref) & !is.na(alt) & grepl("^[ACGT]$", ref) & grepl("^[ACGT]$", alt) & ref != alt
readCALLS <- function(path, ids) {
  calls <- readTSV(path)
  metadata <- c("#CHROM", "CHROM", "CHR", "SNP", "(C)M", "CM", "POS", "COUNTED", "ALT")
  if (!all(c("SNP", "COUNTED", "ALT") %in% names(calls)) || anyDuplicated(calls$SNP) ||
      !setequal(as.character(calls$SNP), ids)) {
    stop("PLINK genotype export must contain exactly the requested source marker IDs.", call. = FALSE)
  }
  sample <- readTSV(option[["psam"]], colClasses = "character")
  setnames(sample, sub("^#", "", names(sample)))
  if (!"IID" %in% names(sample)) stop("Imported PSAM requires IID.", call. = FALSE)
  if (!"FID" %in% names(sample)) sample[, FID := "0"]
  expected <- paste(sample$FID, sample$IID, sep = "_")
  sampleCOLUMN <- setdiff(names(calls), metadata)
  if (anyDuplicated(expected) || !identical(sampleCOLUMN, expected)) {
    stop("PLINK genotype export participant IDs or order disagree with the imported PSAM.", call. = FALSE)
  }
  calls <- calls[match(ids, SNP)]
  dosage <- as.matrix(calls[, lapply(.SD, as.numeric), .SDcols = sampleCOLUMN])
  if (any(!is.na(dosage) & (!is.finite(dosage) | dosage < 0 | dosage > 2))) {
    stop("PLINK genotype export requires missing values or dosages from zero to two.", call. = FALSE)
  }
  list(calls = calls, dosage = dosage)
}

if (option[["action"]] == "match") {
  pvar <- readTSV(option[["pvar"]], skip = "#CHROM", colClasses = "character")
  setnames(pvar, "#CHROM", "CHROM")
  if (!all(c("CHROM", "POS", "ID", "REF", "ALT") %in% names(pvar))) {
    stop("PVAR requires CHROM, POS, ID, REF and ALT.", call. = FALSE)
  }
  if (anyNA(pvar$ID) || anyDuplicated(pvar$ID)) {
    stop("Raw marker identifiers must be present and unique before dbSNP resolution.", call. = FALSE)
  }
  initial <- readTSV(option[["initial-decisions"]], colClasses = "character")
  if (!identical(initial$final_id, pvar$ID)) {
    stop("Initial marker decisions must follow the annotated PVAR source order.", call. = FALSE)
  }
  missing <- readTSV(option[["missingness"]])
  setnames(missing, "#ID", "ID", skip_absent = TRUE)
  setnames(missing, "MISSING_DOSAGE_CT", "MISSING_CT", skip_absent = TRUE)
  if (!all(c("ID", "MISSING_CT", "OBS_CT") %in% names(missing)) || !identical(as.character(missing$ID), pvar$ID)) {
    stop("PLINK variant missingness must contain ID, MISSING_CT and OBS_CT in PVAR order.", call. = FALSE)
  }
  sampleCOUNT <- as.integer(option[["sample-count"]])
  if (sampleCOUNT < 1L || anyNA(missing$OBS_CT) || anyNA(missing$MISSING_CT) ||
      any(missing$MISSING_CT < 0 | missing$OBS_CT < missing$MISSING_CT | missing$OBS_CT > sampleCOUNT)) {
    stop("PLINK missingness counts are invalid for the imported participant count.", call. = FALSE)
  }
  pvar[, source_row := .I]
  pvar[, c("CHROM", "REF", "ALT") := .(sub("^chr", "", CHROM, ignore.case = TRUE), toupper(REF), toupper(ALT))]
  pvar[, POS := as.integer(POS)]
  pvar[, direct_pair := pairKEY(REF, ALT)]
  pvar[, allele_key := pmin(direct_pair, pairKEY(complement(REF), complement(ALT)))]
  callCOUNT <- missing$OBS_CT - missing$MISSING_CT
  yROW <- which(pvar$CHROM %in% c("Y", "24"))
  if (length(yROW)) {
    yCALLS <- readCALLS(option[["y-calls"]], pvar$ID[yROW])
    if (ncol(yCALLS$dosage) != sampleCOUNT) stop("Y export and imported sample counts disagree.", call. = FALSE)
    callCOUNT[yROW] <- rowSums(is.finite(yCALLS$dosage))
  }
  for (column in c("manifest_a", "manifest_b", "top_a", "top_b", "ilmn_strand", "ref_strand", "assay_ref_a", "assay_ref_b")) {
    if (!column %in% names(initial)) initial[, (column) := ""]
  }
  if (!"assay_status" %in% names(initial)) initial[, assay_status := "NOT_SUPPLIED"]
  assayELIGIBLE <- initial$assay_status %in% c("NOT_SUPPLIED", "COMPATIBLE")

  chromosomeMAP <- readTSV(option[["chromosome-map"]], header = FALSE, col.names = c("chromosome", "accession"), colClasses = "character")
  candidate <- data.table(chromosome = character(), position = integer(), candidate = character(), candidate_ref = character(), candidate_alt = character())
  if (file.info(option[["dbsnp-records"]])$size > 0) {
    dbsnp <- readTSV(option[["dbsnp-records"]], header = FALSE,
      col.names = c("accession", "position", "candidate", "reference", "alternate"), colClasses = "character")
    dbsnp[, chromosome := chromosomeMAP$chromosome[match(accession, chromosomeMAP$accession)]]
    dbsnp <- dbsnp[!is.na(chromosome) & grepl("^rs[0-9]+$", candidate)]
    if (nrow(dbsnp)) {
      alt <- strsplit(toupper(dbsnp$alternate), ",", fixed = TRUE)
      candidate <- dbsnp[rep(seq_len(.N), lengths(alt)), .(chromosome, position = as.integer(position), candidate, candidate_ref = toupper(reference))]
      candidate[, candidate_alt := unlist(alt, use.names = FALSE)]
      candidate <- unique(candidate[validPAIR(candidate_ref, candidate_alt)])
    }
  }
  candidate[, direct_pair := pairKEY(candidate_ref, candidate_alt)]
  candidate[, allele_key := pmin(direct_pair, pairKEY(complement(candidate_ref), complement(candidate_alt)))]
  coordinateCOUNT <- candidate[, .(coordinate_candidates = uniqueN(candidate)), by = .(chromosome, position)]
  target <- pvar[validPAIR(REF, ALT) & assayELIGIBLE]
  compatible <- merge(target[, .(source_row, chromosome = CHROM, position = POS, allele_key, source_pair = direct_pair)],
    candidate, by = c("chromosome", "position", "allele_key"), allow.cartesian = TRUE, sort = FALSE)
  compatible[, allele_compatible_candidates := uniqueN(candidate), by = source_row]
  compatible[, direct := source_pair == direct_pair]
  setorderv(compatible, c("source_row", "direct", "candidate_ref", "candidate_alt"), c(1L, -1L, 1L, 1L))
  selected <- compatible[allele_compatible_candidates == 1L][!duplicated(source_row)]

  result <- data.table(
    source_id = initial$source_id, final_id = "", source_chr = initial$source_chr,
    source_pos = initial$source_pos, final_chr = pvar$CHROM, final_pos = pvar$POS,
    source_ref = initial$source_ref, source_alt = initial$source_alt,
    final_ref = pvar$REF, final_alt = pvar$ALT,
    coordinate_candidates = 0L, allele_compatible_candidates = 0L,
    call_count = callCOUNT,
    call_rate = callCOUNT / sampleCOUNT,
    probe_count = 0L, assay_count = 0L, overlap_count = 0L, concordance = NA_real_,
    decision = "EXCLUDED_UNRESOLVED_MARKER", reason = "No dbSNP coordinate match",
    source_row = pvar$source_row, import_id = pvar$ID, match_type = ""
  )
  for (column in c("manifest_a", "manifest_b", "top_a", "top_b", "ilmn_strand", "ref_strand", "assay_ref_a", "assay_ref_b", "assay_status")) {
    result[, (column) := initial[[column]]]
  }
  counts <- merge(target[, .(source_row, chromosome = CHROM, position = POS)],
    coordinateCOUNT, by = c("chromosome", "position"), sort = FALSE)
  result[counts$source_row, coordinate_candidates := counts$coordinate_candidates]
  counts <- unique(compatible[, .(source_row, allele_compatible_candidates)])
  result[counts$source_row, allele_compatible_candidates := counts$allele_compatible_candidates]
  result[coordinate_candidates > 0L, reason := "No allele-compatible dbSNP match"]
  result[allele_compatible_candidates > 1L, reason := "Multiple allele-compatible dbSNP matches"]
  result[!validPAIR(pvar$REF, pvar$ALT), reason := "Incompatible assay alleles"]
  result[!assayELIGIBLE, reason := "Source alleles are incompatible with or absent from the assay manifest"]
  result[pvar$REF %in% c("0", ".") & pvar$ALT %in% c("0", "."), reason := "All genotypes or assay alleles missing"]
  result[selected$source_row, c("final_id", "final_chr", "final_pos", "final_ref", "final_alt", "match_type", "reason") :=
    .(selected$candidate, selected$chromosome, selected$position, selected$candidate_ref, selected$candidate_alt,
      ifelse(selected$direct, "DIRECT", "COMPLEMENT"), "Unique allele-compatible dbSNP match")]
  result[nzchar(final_id) & call_count == 0L, c("final_id", "reason") := .("", "All genotypes missing")]
  result[assay_status == "NOT_SUPPLIED", c("assay_ref_a", "assay_ref_b") := .(final_ref, final_alt)]
  result[, assay_pair := pairKEY(assay_ref_a, assay_ref_b)]
  result[, assay_complement := ifelse(assay_status == "NOT_SUPPLIED", match_type == "COMPLEMENT", ref_strand == "-")]
  result[nzchar(final_id), c("probe_count", "assay_count") := .(.N, uniqueN(assay_pair)), by = final_id]
  result[probe_count == 1L, decision := "RETAINED_UNIQUE"]
  result[probe_count > 1L & assay_count > 1L, c("decision", "reason") :=
    .("EXCLUDED_DIFFERENT_ASSAY_DUPLICATE_GROUP", "Repeated rsID has more than one reference-oriented assay pair")]
  saveRDS(result, option[["candidates"]])
  writeTSV(result[probe_count > 1L & assay_count == 1L, .(import_id)], option[["duplicate-markers"]], FALSE)
} else if (option[["action"]] == "finalise") {
  result <- readRDS(option[["candidates"]])
  duplicates <- result[probe_count > 1L & assay_count == 1L]
  if (nrow(duplicates)) {
    exported <- readCALLS(option[["calls"]], duplicates$import_id)
    calls <- exported$calls
    dosage <- exported$dosage
    counted <- toupper(calls$COUNTED)
    counted[duplicates$assay_complement] <- complement(counted[duplicates$assay_complement])
    anchor <- pmin(duplicates$assay_ref_a, duplicates$assay_ref_b)
    other <- pmax(duplicates$assay_ref_a, duplicates$assay_ref_b)
    if (any(!counted %in% c("A", "C", "G", "T")) || any(counted != anchor & counted != other)) {
      stop("PLINK counted alleles do not match the reference-oriented assay pairs for duplicate probes.", call. = FALSE)
    }
    reverse <- which(counted == other)
    dosage[reverse, ] <- 2 - dosage[reverse, , drop = FALSE]
    observedCOUNT <- rowSums(is.finite(dosage))
    mismatch <- which(observedCOUNT != duplicates$call_count)
    if (length(mismatch)) {
      index <- mismatch[[1L]]
      stop(sprintf("Probe '%s' on chromosome %s has %s exported calls, but marker eligibility recorded %s. Review stored-call counting.",
        duplicates$source_id[index], duplicates$final_chr[index], observedCOUNT[index], duplicates$call_count[index]), call. = FALSE)
    }
    # Group once; only repeated probes need participant-level comparisons.
    groups <- split(seq_len(nrow(duplicates)), duplicates$final_id)
    for (index in groups) {
      sourceROW <- duplicates$source_row[index]
      exactID <- duplicates$source_id[index] == duplicates$final_id[index]
      representative <- index[order(-duplicates$call_rate[index], -as.integer(exactID), sourceROW)][[1L]]
      observed <- is.finite(dosage[index, , drop = FALSE]) &
        matrix(is.finite(dosage[representative, ]), nrow = length(index), ncol = ncol(dosage), byrow = TRUE)
      difference <- abs(sweep(dosage[index, , drop = FALSE], 2L, dosage[representative, ], "-"))
      overlap <- rowSums(observed)
      probeCONCORDANCE <- rowSums(observed & difference < 1e-8, na.rm = TRUE) / overlap
      result[sourceROW, c("overlap_count", "concordance") := .(overlap, probeCONCORDANCE)]
      if (any(!is.finite(probeCONCORDANCE) | probeCONCORDANCE != 1)) {
        result[sourceROW, c("decision", "reason") :=
          .("EXCLUDED_DISCORDANT_DUPLICATE_GROUP", "Same-assay duplicate probes are not completely concordant over observed calls")]
      } else {
        result[sourceROW, c("decision", "reason") :=
          .("EXCLUDED_REDUNDANT_DUPLICATE_PROBE", "Completely concordant duplicate represented by higher-call-rate probe")]
        result[duplicates$source_row[representative], c("decision", "reason") :=
          .("RETAINED_DUPLICATE_REPRESENTATIVE", "Highest call rate, then exact rsID, then original row order")]
      }
    }
  }
  retained <- result[startsWith(decision, "RETAINED_")]
  if (nrow(retained) && (anyNA(retained$final_id) || any(!nzchar(retained$final_id)) ||
      anyDuplicated(retained$final_id) || anyDuplicated(retained$import_id) ||
      any(!validPAIR(retained$final_ref, retained$final_alt)))) {
    stop("Retained markers require unique source and final IDs and distinct A/C/G/T REF and ALT alleles.", call. = FALSE)
  }
  writeTSV(result[, !c("source_row", "import_id", "match_type"), with = FALSE], option[["output-decisions"]])
  if (!nrow(retained)) stop("Raw marker resolution retained no variants; inspect the marker decision reasons.", call. = FALSE)
  writeTSV(retained[, .(import_id)], option[["keep"]], FALSE)
  writeTSV(retained[, .(import_id, final_id)], option[["rename"]], FALSE)
} else {
  stop("Marker resolution action must be match or finalise.", call. = FALSE)
}
