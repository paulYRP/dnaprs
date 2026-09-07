#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))
resolver <- normalizePath("bin/resolve_target_markers.R")
root <- tempfile("target-markers-")
dir.create(root)
writeTSV <- function(value, name, header = TRUE) fwrite(value, file.path(root, name), sep = "\t", quote = FALSE, col.names = header, na = "NA")
run <- function(action, expected = 0L) {
  paths <- c(pvar = "target.pvar", `initial-decisions` = "initial.tsv", missingness = "missing.vmiss",
    `dbsnp-records` = "dbsnp.tsv", `chromosome-map` = "map.tsv", candidates = "candidates.rds",
    `duplicate-markers` = "duplicates.txt", calls = "calls.traw", `output-decisions` = "decisions.tsv",
    keep = "keep.txt", rename = "rename.tsv")
  args <- c(resolver, "--action", action, "--threads", "1", "--sample-count", "4",
    unlist(Map(function(key, value) c(paste0("--", key), file.path(root, value)), names(paths), paths)))
  status <- system2(file.path(R.home("bin"), "Rscript"), shQuote(args), stdout = file.path(root, "run.log"), stderr = file.path(root, "run.err"))
  if (status != expected) stop(paste(readLines(file.path(root, "run.err")), collapse = "\n"))
}

pvar <- data.table(`#CHROM` = "1", POS = c(10, 20, 30, 40, 50, 60, 70, 70, 70, 80, 80, 90, 90, 110, 110, 120, 120, 140),
  ID = c("direct", "complement", "unresolved", "multiallelic", "ambiguous", "invalid", "rs107", "b107", "c107", "b108", "rs108", "d1", "d2", "assay1", "assay2", "no_overlap1", "no_overlap2", "missing"),
  REF = c("A", "T", "A", "A", "A", "A", rep("A", 12)),
  ALT = c("G", "C", "G", "C", "G", "A", rep("G", 8), "C", "G", "G", "G"))
writeTSV(pvar, "target.pvar")
initial <- pvar[, .(source_id = ID, final_id = ID, source_chr = `#CHROM`, source_pos = POS, source_ref = REF, source_alt = ALT)]
writeTSV(initial, "initial.tsv")
missing <- data.table(`#ID` = pvar$ID, MISSING_CT = c(rep(0L, 6), 1L, rep(0L, 8), 3L, 3L, 4L), OBS_CT = 4L)
writeTSV(missing, "missing.vmiss")
writeTSV(data.table(chromosome = "1", accession = "NC_1"), "map.tsv", FALSE)
dbsnp <- data.table(accession = "NC_1", position = c(10, 20, 40, 50, 50, 60, 70, 80, 90, 110, 120, 140),
  candidate = c("rs101", "rs102", "rs104", "rs105", "rs205", "rs106", "rs107", "rs108", "rs109", "rs110", "rs120", "rs140"),
  reference = "A", alternate = c("A,G", "G", "A,C,G", "G", "G", "A", "G", "G", "G", "G,C", "G", "G"))
writeTSV(dbsnp, "dbsnp.tsv", FALSE)
run("match")
candidates <- readRDS(file.path(root, "candidates.rds"))
duplicate <- readLines(file.path(root, "duplicates.txt"))
stopifnot(identical(duplicate, c("rs107", "b107", "c107", "b108", "rs108", "d1", "d2", "no_overlap1", "no_overlap2")))
calls <- data.table(CHR = "1", SNP = duplicate, CM = 0, POS = pvar$POS[match(duplicate, pvar$ID)],
  COUNTED = c("G", "A", rep("G", 7)), ALT = c("A", "G", rep("A", 7)),
  S1 = c(0, 2, 0, 0, 0, 0, 2, 0, NA),
  S2 = c(1, 1, 1, 1, 1, 1, 1, NA, NA),
  S3 = c(2, 0, 2, 2, 2, 2, 2, NA, 2),
  S4 = c(NA, 2, 0, 0, 0, 0, 0, NA, NA))
writeTSV(calls, "calls.traw")
run("finalise")
decision <- fread(file.path(root, "decisions.tsv"))
stopifnot(
  identical(decision$source_id, pvar$ID),
  identical(readLines(file.path(root, "keep.txt")), c("direct", "complement", "multiallelic", "b107", "rs108")),
  decision[source_id == "direct", final_ref] == "A",
  decision[source_id == "direct", final_alt] == "G",
  candidates[source_id == "complement", match_type] == "COMPLEMENT",
  decision[source_id == "unresolved", reason] == "No dbSNP coordinate match",
  decision[source_id == "ambiguous", allele_compatible_candidates] == 2L,
  decision[source_id == "invalid", reason] == "Incompatible assay alleles",
  all(decision[source_id %in% c("d1", "d2", "no_overlap1", "no_overlap2"), decision] == "EXCLUDED_DISCORDANT_DUPLICATE_GROUP"),
  all(decision[source_id %in% c("assay1", "assay2"), decision] == "EXCLUDED_DIFFERENT_ASSAY_DUPLICATE_GROUP"),
  decision[source_id == "missing", reason] == "All genotypes missing",
  decision[source_id == "b107", overlap_count] == 4L,
  all(decision[startsWith(decision, "RETAINED_"), final_ref != final_alt])
)

# Empty candidate queries retain the complete exclusion audit before failing.
writeLines(character(), file.path(root, "dbsnp.tsv"))
run("match")
run("finalise", 1L)
stopifnot(nrow(fread(file.path(root, "decisions.tsv"))) == nrow(pvar))
writeTSV(data.table(accession = "NC_1", position = 10, candidate = ".", reference = "A", alternate = "G"), "dbsnp.tsv", FALSE)
run("match")
run("finalise", 1L)
stopifnot(nrow(fread(file.path(root, "decisions.tsv"))) == nrow(pvar))

# Unique markers do not require a genotype export, including larger synthetic inputs.
count <- 10000L
pvar <- data.table(`#CHROM` = "1", POS = seq_len(count), ID = paste0("probe", seq_len(count)), REF = "A", ALT = "G")
writeTSV(pvar, "target.pvar")
writeTSV(pvar[, .(source_id = ID, final_id = ID, source_chr = `#CHROM`, source_pos = POS, source_ref = REF, source_alt = ALT)], "initial.tsv")
writeTSV(data.table(`#ID` = pvar$ID, MISSING_CT = 0L, OBS_CT = 4L), "missing.vmiss")
writeTSV(data.table(accession = "NC_1", position = pvar$POS, candidate = paste0("rs", seq_len(count)), reference = "A", alternate = "G"), "dbsnp.tsv", FALSE)
run("match")
stopifnot(file.info(file.path(root, "duplicates.txt"))$size == 0)
run("finalise")
stopifnot(length(readLines(file.path(root, "keep.txt"))) == count)
message("Marker joins, invalid pairs, restricted exports and duplicate selection passed.")
