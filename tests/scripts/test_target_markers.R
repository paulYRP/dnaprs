#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))
resolver <- normalizePath("bin/resolve_target_markers.R")
root <- tempfile("target-markers-")
dir.create(root)
writeTSV <- function(value, name, header = TRUE) fwrite(value, file.path(root, name), sep = "\t", quote = FALSE, col.names = header, na = "NA")
run <- function(action, expected = 0L, threads = 1L) {
  paths <- c(pvar = "target.pvar", `initial-decisions` = "initial.tsv", missingness = "missing.vmiss",
    psam = "target.psam", `y-calls` = "y.traw",
    `dbsnp-records` = "dbsnp.tsv", `chromosome-map` = "map.tsv", candidates = "candidates.rds",
    `duplicate-markers` = "duplicates.txt", calls = "calls.traw", `output-decisions` = "decisions.tsv",
    keep = "keep.txt", rename = "rename.tsv")
  args <- c(resolver, "--action", action, "--threads", as.character(threads), "--sample-count", "4",
    unlist(Map(function(key, value) c(paste0("--", key), file.path(root, value)), names(paths), paths)))
  status <- system2(file.path(R.home("bin"), "Rscript"), shQuote(args), stdout = file.path(root, "run.log"), stderr = file.path(root, "run.err"))
  if (status != expected) stop(paste(readLines(file.path(root, "run.err")), collapse = "\n"))
}
writeTSV(data.table(`#FID` = "0", IID = paste0("S", 1:4), SEX = c(1, 2, 2, 0)), "target.psam")

# Candidate compatibility uses the complete record, not individual REF/ALT pairs.
writeTSV(data.table(chromosome = "1", accession = "NC_1"), "map.tsv", FALSE)
edge <- data.table(`#CHROM` = "1", POS = c(10L, 20L, 30L, 40L, 50L, 60L),
  ID = c("alt_pair", "mixed_record", "new_ambiguity", "repeated_rsid", "manifest_only", "absent"),
  REF = c("C", "A", "A", "A", "A", "A"), ALT = c("G", "G", "G", "G", "C", "G"))
edgeINITIAL <- edge[, .(source_id = ID, final_id = ID, source_chr = `#CHROM`, source_pos = POS,
  source_ref = REF, source_alt = ALT, top_a = REF, top_b = ALT,
  assay_ref_a = REF, assay_ref_b = ALT, assay_status = "COMPATIBLE", ref_strand = "+")]
edgeINITIAL[source_id == "manifest_only", `:=`(top_b = "G", assay_ref_b = "G", assay_status = "INCOMPATIBLE")]
edge <- rbind(edge, data.table(`#CHROM` = "1", POS = 70L, ID = "both_strands", REF = "T", ALT = "G"))
edgeINITIAL <- rbind(edgeINITIAL, data.table(source_id = "both_strands", final_id = "both_strands",
  source_chr = "1", source_pos = 70L, source_ref = "T", source_alt = "G", top_a = "", top_b = "",
  assay_ref_a = "", assay_ref_b = "", assay_status = "NOT_SUPPLIED", ref_strand = ""))
writeTSV(edge, "target.pvar")
writeTSV(edgeINITIAL, "initial.tsv")
writeTSV(data.table(`#ID` = edge$ID, MISSING_CT = 0L, OBS_CT = 4L), "missing.vmiss")
edgeDB <- data.table(accession = "NC_1", position = c(10L, 20L, 30L, 30L, 40L, 40L, 40L, 50L),
  candidate = c("rs10", "rs20", "rs30", "rs31", "rs40", "rs40", "rs40", "rs50"),
  reference = c("A", "A", "A", "C", "A", "A", "A", "A"),
  alternate = c("C,G", "G,AT", "G", "A,G", "G", "G", "C,G,T", "G"))
edgeDB <- rbind(edgeDB, data.table(accession = "NC_1", position = 70L, candidate = "rs70", reference = "A", alternate = "C,G,T"))
writeTSV(edgeDB, "dbsnp.tsv", FALSE)
run("match")
edgeRESULT <- readRDS(file.path(root, "candidates.rds"))
stopifnot(
  edgeRESULT[source_id == "alt_pair", allele_compatible_candidates] == 1L,
  edgeRESULT[source_id == "mixed_record", coordinate_candidates] == 0L,
  edgeRESULT[source_id == "new_ambiguity", allele_compatible_candidates] == 2L,
  edgeRESULT[source_id == "repeated_rsid", allele_compatible_candidates] == 1L,
  edgeRESULT[source_id == "manifest_only", allele_compatible_candidates] == 1L,
  !startsWith(edgeRESULT[source_id == "manifest_only", decision], "RETAINED_"),
  edgeRESULT[source_id == "absent", coordinate_candidates] == 0L,
  edgeRESULT[source_id == "both_strands", match_type] == "COMPLEMENT",
  edgeRESULT[source_id == "both_strands", assay_ref_a] == "A",
  edgeRESULT[source_id == "both_strands", assay_ref_b] == "C"
)
writeTSV(edgeDB[rev(seq_len(.N))], "dbsnp.tsv", FALSE)
run("match", threads = 2L)
orderCHECK <- all.equal(as.data.frame(edgeRESULT), as.data.frame(readRDS(file.path(root, "candidates.rds"))))
if (!isTRUE(orderCHECK)) stop(paste(orderCHECK, collapse = "; "))

# Exhaustive base-set membership supplies independent expected mapping decisions.
bases <- c("A", "C", "G", "T")
pairs <- combn(bases, 2L, simplify = FALSE)
sets <- unlist(lapply(2:4, function(size) combn(bases, size, simplify = FALSE)), recursive = FALSE)
grid <- CJ(pair = seq_along(pairs), record = seq_along(sets))
for (reverseSTRAND in c(FALSE, TRUE)) {
  query <- lapply(grid$pair, function(index) pairs[[index]])
  records <- lapply(grid$record, function(index) sets[[index]])
  if (reverseSTRAND) query <- lapply(query, function(value) chartr("ACGT", "TGCA", value))
  expected <- mapply(function(pair, record) all(pair %in% record) || all(chartr("ACGT", "TGCA", pair) %in% record), query, records)
  pvar <- data.table(`#CHROM` = "1", POS = seq_len(nrow(grid)), ID = paste0("pair", seq_len(nrow(grid))),
    REF = vapply(query, `[`, character(1L), 1L), ALT = vapply(query, `[`, character(1L), 2L))
  writeTSV(pvar, "target.pvar")
  writeTSV(pvar[, .(source_id = ID, final_id = ID, source_chr = `#CHROM`, source_pos = POS, source_ref = REF, source_alt = ALT)], "initial.tsv")
  writeTSV(data.table(`#ID` = pvar$ID, MISSING_CT = 0L, OBS_CT = 4L), "missing.vmiss")
  writeTSV(data.table(accession = "NC_1", position = pvar$POS, candidate = paste0("rs", pvar$POS),
    reference = vapply(records, `[`, character(1L), 1L),
    alternate = vapply(records, function(value) paste(value[-1L], collapse = ","), character(1L))), "dbsnp.tsv", FALSE)
  run("match")
  result <- readRDS(file.path(root, "candidates.rds"))
  stopifnot(identical(result$allele_compatible_candidates, as.integer(expected)), identical(result$source_id, pvar$ID))
}

# Single-zero recovery requires a compatible observed TOP allele.
writeTSV(data.table(`#CHROM` = "1", POS = 10:13, ID = c("recover", "incompatible", "reverse", "unlisted"),
  REF = c("A", "C", "A", "A"), ALT = c(".", ".", "G", "G")), "assays.pvar")
writeLines(c("[Assay]", "Name,IlmnStrand,SNP,Chr,MapInfo,RefStrand",
  "recover,BOT,[T/C],1,10,-", "incompatible,TOP,[A/G],1,11,+", "reverse,TOP,[A/G],1,12,+"), file.path(root, "assays.csv"))
stopifnot(system2("perl", shQuote(c(normalizePath("bin/target_adapter.pl"), "annotate-pvar",
  file.path(root, "assays.pvar"), file.path(root, "assays.csv"), "",
  file.path(root, "annotated.pvar"), file.path(root, "assay_decisions.tsv")))) == 0L)
assay <- fread(file.path(root, "assay_decisions.tsv"))
stopifnot(assay[source_id == "recover", final_alt] == "G",
  assay[source_id == "recover", assay_ref_a] == "T", assay[source_id == "recover", assay_ref_b] == "C",
  assay[source_id == "incompatible", assay_status] == "INCOMPATIBLE",
  assay[source_id == "unlisted", assay_status] == "MISSING_MANIFEST_RECORD")

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
setnames(calls, paste0("S", 1:4), paste0("0_S", 1:4))
writeTSV(calls, "calls.traw")
run("finalise")
decision <- fread(file.path(root, "decisions.tsv"))
stopifnot(
  identical(decision$source_id, pvar$ID),
  identical(readLines(file.path(root, "keep.txt")), c("direct", "complement", "multiallelic", "b107", "rs108")),
  candidates[source_id == "direct", candidate_ref] == "A",
  candidates[source_id == "direct", candidate_alt] == "A,G",
  candidates[source_id == "complement", match_type] == "COMPLEMENT",
  decision[source_id == "unresolved", reason] == "No dbSNP coordinate match",
  decision[source_id == "ambiguous", allele_compatible_candidates] == 2L,
  decision[source_id == "invalid", reason] == "Incompatible assay alleles",
  all(decision[source_id %in% c("d1", "d2", "no_overlap1", "no_overlap2"), decision] == "EXCLUDED_DISCORDANT_DUPLICATE_GROUP"),
  all(decision[source_id %in% c("assay1", "assay2"), decision] == "EXCLUDED_DIFFERENT_ASSAY_DUPLICATE_GROUP"),
  decision[source_id == "missing", reason] == "All genotypes missing",
  decision[source_id == "b107", overlap_count] == 4L,
  all(is.na(decision[startsWith(decision, "RETAINED_"), final_ref]))
)

# Count all stored Y calls before eligibility and use manifest assay identity for groups.
pvar <- data.table(`#CHROM` = c("1", rep("Y", 5), "1", "1"), POS = c(10, 70, 70, 140, 150, 160, 80, 80),
  ID = c("direct", "rs107", "better_y", "unique_y", "missing_y", "unknown_y", "assay_plus", "assay_minus"), REF = "A", ALT = "G")
writeTSV(pvar, "target.pvar")
initial <- pvar[, .(source_id = ID, final_id = ID, source_chr = `#CHROM`, source_pos = POS, source_ref = REF, source_alt = ALT,
  assay_status = "COMPATIBLE", top_a = "A", top_b = "G", assay_ref_a = "A", assay_ref_b = "G", ref_strand = "+")]
initial[source_id == "assay_minus", c("assay_ref_a", "assay_ref_b", "ref_strand") := .("T", "C", "-")]
writeTSV(initial, "initial.tsv")
writeTSV(data.table(`#ID` = pvar$ID, MISSING_CT = c(0, 0, 0, 1, 1, 1, 0, 0), OBS_CT = c(4, rep(1, 5), 4, 4)), "missing.vmiss")
writeTSV(data.table(chromosome = c("1", "Y"), accession = c("NC_1", "NC_Y")), "map.tsv", FALSE)
writeTSV(data.table(accession = c("NC_1", rep("NC_Y", 4), "NC_1"), position = c(10, 70, 140, 150, 160, 80),
  candidate = c("rs101", "rs107", "rs140", "rs150", "rs160", "rs108"), reference = "A", alternate = "G"), "dbsnp.tsv", FALSE)
y <- data.table(CHR = "Y", SNP = pvar$ID[2:6], CM = 0, POS = pvar$POS[2:6], COUNTED = "A", ALT = "G",
  `0_S1` = c(2, 2, NA, NA, NA), `0_S2` = c(NA, 2, 2, NA, NA), `0_S3` = NA_real_, `0_S4` = c(NA, NA, NA, NA, 0))
writeTSV(y, "y.traw")
run("match")
stopifnot(identical(readLines(file.path(root, "duplicates.txt")), c("rs107", "better_y")))
writeTSV(y[1:2], "calls.traw")
run("finalise")
decision <- fread(file.path(root, "decisions.tsv"))
stopifnot(
  decision[source_id == "better_y", call_count] == 2,
  decision[source_id == "better_y", call_rate] == 0.5,
  decision[source_id == "better_y", decision] == "RETAINED_DUPLICATE_REPRESENTATIVE",
  all(decision[source_id %in% c("unique_y", "unknown_y"), decision] == "RETAINED_UNIQUE"),
  decision[source_id == "missing_y", reason] == "All genotypes missing",
  all(decision[source_id %in% c("assay_plus", "assay_minus"), decision] == "EXCLUDED_DIFFERENT_ASSAY_DUPLICATE_GROUP")
)
bad <- copy(y[1:2]); setcolorder(bad, c(names(bad)[1:6], "0_S2", "0_S1", "0_S3", "0_S4"))
writeTSV(bad, "calls.traw")
run("finalise", 1L)
stopifnot(any(grepl("participant IDs or order", readLines(file.path(root, "run.err")))))
writeTSV(y[1:2], "calls.traw")

# Complemented probes can share one manifest assay while their TOP alleles differ.
pvar[`#CHROM` == "1" & ID == "assay_minus", c("REF", "ALT") := .("T", "C")]
writeTSV(pvar, "target.pvar")
initial[source_id == "assay_minus", c("assay_ref_a", "assay_ref_b", "top_a", "top_b") := .("A", "G", "T", "C")]
writeTSV(initial, "initial.tsv")
run("match")
assayCALLS <- data.table(CHR = "1", SNP = c("assay_plus", "assay_minus"), CM = 0, POS = 80,
  COUNTED = c("A", "C"), ALT = c("G", "T"), `0_S1` = c(2, 0), `0_S2` = c(1, 1), `0_S3` = c(0, 2), `0_S4` = c(2, 0))
writeTSV(rbind(y[1:2], assayCALLS), "calls.traw")
run("finalise")
decision <- fread(file.path(root, "decisions.tsv"))
stopifnot(decision[source_id == "assay_plus", decision] == "RETAINED_DUPLICATE_REPRESENTATIVE",
  all(decision[source_id %in% c("assay_plus", "assay_minus"), concordance] == 1))

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
writeTSV(data.table(chromosome = "1", accession = "NC_1"), "map.tsv", FALSE)
writeTSV(pvar[, .(source_id = ID, final_id = ID, source_chr = `#CHROM`, source_pos = POS, source_ref = REF, source_alt = ALT)], "initial.tsv")
writeTSV(data.table(`#ID` = pvar$ID, MISSING_CT = 0L, OBS_CT = 4L), "missing.vmiss")
writeTSV(data.table(accession = "NC_1", position = pvar$POS, candidate = paste0("rs", seq_len(count)), reference = "A", alternate = "G"), "dbsnp.tsv", FALSE)
run("match")
stopifnot(file.info(file.path(root, "duplicates.txt"))$size == 0)
run("finalise")
stopifnot(length(readLines(file.path(root, "keep.txt"))) == count)
message("Marker joins, invalid pairs, restricted exports and duplicate selection passed.")
