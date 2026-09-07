#!/usr/bin/env Rscript

# Usage: Rscript tests/scripts/benchmark_target_markers.R [markers] [biallelic|mixed] [resolver]
suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE)
count <- if (length(args)) as.integer(args[[1L]]) else 700000L
mode <- if (length(args) >= 2L) args[[2L]] else "biallelic"
resolver <- normalizePath(if (length(args) >= 3L) args[[3L]] else "bin/resolve_target_markers.R")
stopifnot(is.finite(count), count > 0L, mode %in% c("biallelic", "mixed"))
setDTthreads(1L)
root <- tempfile("marker-benchmark-")
dir.create(root)
writeTSV <- function(value, name, header = TRUE) fwrite(value, file.path(root, name), sep = "\t", quote = FALSE, col.names = header)
pvar <- data.table(`#CHROM` = "1", POS = seq_len(count), ID = paste0("probe", seq_len(count)), REF = "A", ALT = "G")
if (mode == "mixed") pvar[POS %% 2L == 0L, REF := "C"]
writeTSV(pvar, "target.pvar")
writeTSV(pvar[, .(source_id = ID, final_id = ID, source_chr = `#CHROM`, source_pos = POS, source_ref = REF, source_alt = ALT)], "initial.tsv")
writeTSV(data.table(`#ID` = pvar$ID, MISSING_CT = 0L, OBS_CT = 4L), "missing.vmiss")
writeTSV(data.table(chromosome = "1", accession = "NC_1"), "map.tsv", FALSE)
records <- data.table(accession = "NC_1", position = pvar$POS, candidate = paste0("rs", pvar$POS), reference = "A", alternate = "G")
if (mode == "mixed") records[position %% 2L == 0L, alternate := "C,G"]
writeTSV(records, "dbsnp.tsv", FALSE)
paths <- c(pvar = "target.pvar", `initial-decisions` = "initial.tsv", missingness = "missing.vmiss",
  `dbsnp-records` = "dbsnp.tsv", `chromosome-map` = "map.tsv", candidates = "candidates.rds", `duplicate-markers` = "duplicates.txt")
resolverARGS <- c("--action", "match", "--threads", "1", "--sample-count", "4",
  unlist(Map(function(key, value) c(paste0("--", key), file.path(root, value)), names(paths), paths)))
rm(pvar, records)
invisible(gc())
environment <- new.env(parent = globalenv())
environment$commandArgs <- function(trailingOnly = FALSE) resolverARGS
elapsed <- system.time(sys.source(resolver, envir = environment))[["elapsed"]]
result <- readRDS(file.path(root, "candidates.rds"))
stopifnot(nrow(result) == count, all(result$allele_compatible_candidates == 1L),
  all(result$decision == "RETAINED_UNIQUE"), file.info(file.path(root, "duplicates.txt"))$size == 0)
memory <- if (file.exists("/proc/self/status")) {
  sub("^VmHWM:[[:space:]]*", "", grep("^VmHWM:", readLines("/proc/self/status"), value = TRUE))
} else NA_character_
print(data.table(markers = count, mode, threads = 1L, elapsed_seconds = elapsed, peak_rss = memory,
  resolver = basename(resolver), retained = nrow(result), status = "PASS"))
