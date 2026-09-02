#!/usr/bin/env Rscript

argument <- commandArgs(trailingOnly = TRUE)
if (length(argument) %% 2L != 0L) stop("Arguments must be --name value pairs.", call. = FALSE)
option <- stats::setNames(
  as.list(argument[seq.int(2L, length(argument), 2L)]),
  sub("^--", "", argument[seq.int(1L, length(argument), 2L)])
)
if (!requireNamespace("data.table", quietly = TRUE)) stop("The data.table package is required.", call. = FALSE)

complement <- function(value) chartr("ACGT", "TGCA", toupper(value))
gwas <- data.table::fread(option[["cojo"]])
reference <- data.table::fread(option[["reference-pvar"]], skip = "#CHROM")
data.table::setnames(reference, "#CHROM", "CHR", skip_absent = TRUE)

requiredGWAS <- c("SNP", "CHR", "BP", "A1", "A2", "freq", "b", "se", "p", "N")
requiredREFERENCE <- c("CHR", "POS", "ID", "REF", "ALT")
if (!all(requiredGWAS %in% names(gwas))) stop("Harmonised GWAS is missing required COJO columns.", call. = FALSE)
if (!all(requiredREFERENCE %in% names(reference))) stop("PLINK reference PVAR is missing CHR, POS, ID, REF, or ALT.", call. = FALSE)

gwas[, `:=`(
  CHR = sub("^chr", "", as.character(CHR), ignore.case = TRUE),
  BP = as.integer(BP), A1 = toupper(A1), A2 = toupper(A2), source_row = .I
)]
reference <- reference[
  ID != "" & ID != "." & grepl("^[ACGT]$", REF) & grepl("^[ACGT]$", ALT),
  .(
    CHR = sub("^chr", "", as.character(CHR), ignore.case = TRUE),
    BP = as.integer(POS), reference_id = as.character(ID),
    reference_ref = toupper(REF), reference_alt = toupper(ALT)
  )
]
if (nrow(reference) == 0L) stop("PLINK reference PVAR contains no named biallelic SNVs.", call. = FALSE)
if (anyDuplicated(reference$reference_id)) stop("PLINK reference PVAR variant IDs must be unique.", call. = FALSE)

candidate <- merge(
  gwas[, .(source_row, source_id = SNP, CHR, BP, source_a1 = A1, source_a2 = A2)],
  reference,
  by = c("CHR", "BP"), all = FALSE, allow.cartesian = TRUE, sort = FALSE
)
candidate[, direct :=
  (source_a1 == reference_ref & source_a2 == reference_alt) |
  (source_a1 == reference_alt & source_a2 == reference_ref)
]
candidate[, strand :=
  (complement(source_a1) == reference_ref & complement(source_a2) == reference_alt) |
  (complement(source_a1) == reference_alt & complement(source_a2) == reference_ref)
]
candidate <- candidate[direct | strand]
candidate[, rank := data.table::fcase(
  source_id == reference_id, 0L,
  direct, 1L,
  default = 2L
)]
candidate <- candidate[, .SD[rank == min(rank)], by = source_row]

compatibleROWS <- unique(candidate$source_row)
ambiguousROWS <- candidate[, .(references = data.table::uniqueN(reference_id)), by = source_row][references != 1L, source_row]
selected <- candidate[!source_row %in% ambiguousROWS]
selected <- selected[, .SD[1L], by = source_row]
data.table::setorder(selected, source_row)

aligned <- gwas[match(selected$source_row, source_row)]
aligned[, SNP := selected$reference_id]
complemented <- !selected$direct
aligned[complemented, `:=`(A1 = complement(A1), A2 = complement(A2))]
duplicate <- duplicated(aligned$SNP) | duplicated(aligned$SNP, fromLast = TRUE)
duplicateCOUNT <- sum(duplicate)
aligned <- aligned[!duplicate]
aligned[, source_row := NULL]
if (nrow(aligned) == 0L) stop("No GWAS variants aligned uniquely to the PLINK LD reference.", call. = FALSE)

data.table::setorder(aligned, CHR, BP)
data.table::fwrite(aligned, option[["output-cojo"]], sep = "\t", quote = FALSE, na = "NA")
data.table::fwrite(
  aligned[, .(ID = SNP, CHR, POS = BP, A1, P = p)],
  option[["output-clump"]], sep = "\t", quote = FALSE
)
data.table::fwrite(
  data.table::data.table(
    trait_id = option[["trait-id"]],
    prs_name = option[["prs-name"]],
    input_variants = nrow(gwas),
    reference_aligned_variants = nrow(aligned),
    filtered_reference_missing = nrow(gwas) - length(compatibleROWS),
    filtered_reference_ambiguous = length(ambiguousROWS),
    filtered_reference_duplicate = duplicateCOUNT,
    complemented_alleles = sum(complemented[!duplicate]),
    status = "PASS"
  ),
  option[["output-qc"]], sep = "\t", quote = FALSE
)
