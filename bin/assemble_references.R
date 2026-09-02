#!/usr/bin/env Rscript

argument <- commandArgs(trailingOnly = TRUE)
name <- sub("^--", "", argument[seq.int(1L, length(argument), 2L)])
value <- argument[seq.int(2L, length(argument), 2L)]
option <- stats::setNames(as.list(value), name)
for (package in c("jsonlite", "openssl")) {
  if (!requireNamespace(package, quietly = TRUE)) stop(sprintf("The %s package is required.", package), call. = FALSE)
}

readTSV <- function(path) utils::read.delim(
  path, sep = "\t", header = TRUE, quote = "", comment.char = "",
  check.names = FALSE, stringsAsFactors = FALSE
)
writeTSV <- function(value, path) utils::write.table(
  value, path, sep = "\t", quote = FALSE, row.names = FALSE, na = ""
)
sha256FILE <- function(path) {
  connection <- file(path, "rb")
  on.exit(close(connection))
  tolower(as.character(openssl::sha256(connection)))
}
resourceFILES <- function(path, companion = "") {
  value <- if (dir.exists(path)) list.files(path, recursive = TRUE, full.names = TRUE, all.files = FALSE, no.. = TRUE) else path
  if (nzchar(companion)) value <- c(value, companion)
  unique(value[file.exists(value) & !dir.exists(value)])
}
resourceDIGEST <- function(files, base) {
  files <- sort(normalizePath(files, winslash = "/", mustWork = TRUE))
  if (length(files) == 1L) return(sha256FILE(files))
  base <- normalizePath(base, winslash = "/", mustWork = TRUE)
  relative <- ifelse(startsWith(files, paste0(base, "/")), substring(files, nchar(base) + 2L), basename(files))
  inventory <- paste(vapply(files, sha256FILE, character(1L)), relative, sep = "  ")
  tolower(as.character(openssl::sha256(charToRaw(paste0(paste(inventory, collapse = "\n"), "\n")))))
}

provided <- readTSV(option[["provided"]])
assetJSON <- rawToChar(jsonlite::base64_dec(option[["assets-json"]]))
asset <- jsonlite::fromJSON(assetJSON, simplifyDataFrame = TRUE)
if (!is.data.frame(asset)) asset <- as.data.frame(asset, stringsAsFactors = FALSE)
asset <- asset[asset$asset_id != "cache_complete", , drop = FALSE]
bundleROOT <- normalizePath("reference_bundle", winslash = "/", mustWork = FALSE)
dir.create(bundleROOT, recursive = TRUE, showWarnings = FALSE)

for (row in seq_len(nrow(asset))) {
  source <- file.path(option[["asset-dir"]], asset$asset_id[[row]])
  if (!file.exists(source)) stop(sprintf("Downloaded asset is missing: %s", asset$asset_id[[row]]), call. = FALSE)
  destination <- file.path(bundleROOT, asset$relative_path[[row]])
  dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
  if (!file.symlink(normalizePath(source, winslash = "/", mustWork = TRUE), destination)) {
    if (!file.copy(source, destination, overwrite = TRUE)) {
      stop(sprintf("Could not assemble reference asset: %s", asset$asset_id[[row]]), call. = FALSE)
    }
  }
}

fastaGZIP <- file.path(bundleROOT, "fasta", "human_g1k_v37.fasta.gz")
fasta <- file.path(bundleROOT, "fasta", "human_g1k_v37.fasta")
if (file.exists(fastaGZIP) && !file.exists(fasta)) {
  input <- gzfile(fastaGZIP, "rb")
  output <- file(fasta, "wb")
  repeat {
    block <- readBin(input, what = "raw", n = 8L * 1024L * 1024L)
    if (length(block) == 0L) break
    writeBin(block, output)
  }
  close(input)
  close(output)
}

downloadedROLE <- unique(asset$reference_type)
downloadedROLE <- downloadedROLE[nzchar(downloadedROLE)]
provided <- provided[!provided$reference_type %in% downloadedROLE, , drop = FALSE]
referenceROW <- lapply(downloadedROLE, function(role) {
  selected <- asset[asset$reference_type == role, , drop = FALSE]
  first <- selected[1L, , drop = FALSE]
  path <- file.path(bundleROOT, first$reference_path[[1L]])
  companion <- if (nzchar(first$companion[[1L]])) file.path(bundleROOT, first$companion[[1L]]) else ""
  files <- resourceFILES(path, companion)
  if (length(files) == 0L) stop(sprintf("Assembled reference role is empty: %s", role), call. = FALSE)
  digestBASE <- if (dir.exists(path)) path else dirname(path)
  data.frame(
    bundle_id = paste0("dnaprs_", option[["genome-build"]]),
    bundle_version = option[["bundle"]],
    reference_id = paste0(toupper(role), "_", option[["genome-build"]]),
    reference_type = role,
    path = normalizePath(path, winslash = "/", mustWork = TRUE),
    companion = if (nzchar(companion)) normalizePath(companion, winslash = "/", mustWork = TRUE) else "",
    build = option[["genome-build"]],
    ancestry = first$ancestry[[1L]],
    version = first$version[[1L]],
    checksum = resourceDIGEST(files, digestBASE),
    source_format = first$source_format[[1L]],
    reference_stage = "source",
    stringsAsFactors = FALSE
  )
})
downloaded <- if (length(referenceROW) > 0L) do.call(rbind, referenceROW) else provided[0L, , drop = FALSE]
reference <- rbind(provided, downloaded)
reference <- reference[order(reference$reference_type), , drop = FALSE]
if (anyDuplicated(reference$reference_type)) stop("The assembled reference bundle contains duplicate roles.", call. = FALSE)

dbsnp <- reference$path[reference$reference_type == "dbsnp"]
if (length(dbsnp) == 1L) {
  required <- c("assembly_report.txt", "GCF_000001405.25.gz", "GCF_000001405.25.gz.tbi")
  if (!all(file.exists(file.path(dbsnp, required)))) stop("The assembled dbSNP source is incomplete.", call. = FALSE)
}
panel <- reference$path[reference$reference_type == "imputation_panel"]
if (length(panel) == 1L && length(list.files(panel, pattern = "\\.bref3$")) != 22L) {
  stop("The assembled imputation panel must contain 22 BREF3 files.", call. = FALSE)
}
writeTSV(reference, "references.tsv")

receipt <- if (nrow(asset) > 0L) data.frame(
  asset_id = asset$asset_id,
  reference_type = asset$reference_type,
  source_url = asset$url,
  checksum_algorithm = asset$checksum_algorithm,
  expected_checksum = asset$checksum,
  observed_checksum = tolower(asset$checksum),
  expected_size = asset$size,
  cache_path = file.path(option[["cache-root"]], asset$relative_path),
  status = "AUTHENTICATED",
  stringsAsFactors = FALSE
) else data.frame(
  asset_id = character(), reference_type = character(), source_url = character(),
  checksum_algorithm = character(), expected_checksum = character(), observed_checksum = character(), expected_size = character(),
  cache_path = character(), status = character(), stringsAsFactors = FALSE
)
writeTSV(receipt, "reference_receipt.tsv")
