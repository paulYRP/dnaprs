#!/usr/bin/env Rscript

# Parse named command-line arguments.
argument <- commandArgs(trailingOnly = TRUE)
if (length(argument) %% 2L != 0L || any(!startsWith(argument[seq.int(1L, length(argument), 2L)], "--"))) {
  stop("Arguments must be supplied as --name value pairs.", call. = FALSE)
}
argumentNAME <- sub("^--", "", argument[seq.int(1L, length(argument), 2L)])
argumentVALUE <- argument[seq.int(2L, length(argument), 2L)]
option <- stats::setNames(as.list(argumentVALUE), argumentNAME)
if (is.null(option[["target-imputation"]])) option[["target-imputation"]] <- "false"

# Define compact helpers for tabular manifests and clear failures.
readTABLE <- function(path) {
  utils::read.delim(
    path,
    header = TRUE,
    sep = "\t",
    quote = "\"",
    comment.char = "",
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

writeTABLE <- function(value, path) {
  utils::write.table(
    value,
    path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    na = ""
  )
}

requireCOLUMN <- function(value, required, label) {
  missingCOLUMN <- setdiff(required, names(value))
  if (length(missingCOLUMN) > 0L) {
    stop(
      sprintf("%s is missing required column(s): %s", label, paste(missingCOLUMN, collapse = ", ")),
      call. = FALSE
    )
  }
}

requireVALUE <- function(value, column, label) {
  missingVALUE <- is.na(value[[column]]) | trimws(as.character(value[[column]])) == ""
  if (any(missingVALUE)) {
    stop(sprintf("%s has an empty value in '%s'.", label, column), call. = FALSE)
  }
}

resolvePATH <- function(path, base) {
  path <- as.character(path)
  if (length(path) == 0L) {
    return(character())
  }
  if (anyNA(path)) {
    stop("A required path is missing.", call. = FALSE)
  }
  path <- path.expand(path)
  if (.Platform$OS.type == "windows") {
    wslPATH <- grepl("^/mnt/[A-Za-z]/", path)
    path[wslPATH] <- paste0(
      toupper(substr(path[wslPATH], 6L, 6L)),
      ":/",
      substring(path[wslPATH], 8L)
    )
  }
  isABSOLUTE <- grepl("^(/|[A-Za-z]:[/\\\\])", path)
  resolved <- ifelse(isABSOLUTE, path, file.path(base, path))
  normalizePath(resolved, winslash = "/", mustWork = FALSE)
}

pipelinePATH <- function(path, base) {
  path <- gsub("\\\\", "/", path.expand(as.character(path)))
  base <- gsub("\\\\", "/", as.character(base))
  isABSOLUTE <- grepl("^(/|[A-Za-z]:/)", path)
  separator <- if (endsWith(base, "/")) "" else "/"
  ifelse(isABSOLUTE, path, paste0(base, separator, path))
}

checkPATH <- function(path, label) {
  missingPATH <- !file.exists(path) & !dir.exists(path)
  if (any(missingPATH)) {
    stop(
      sprintf("%s does not exist: %s", label, paste(path[missingPATH], collapse = ", ")),
      call. = FALSE
    )
  }
}

sha256FILE <- function(path) {
  if (nzchar(Sys.which("sha256sum"))) {
    value <- system2("sha256sum", shQuote(path), stdout = TRUE, stderr = TRUE)
    if (!is.null(attr(value, "status")) || length(value) != 1L) {
      stop(sprintf("Could not calculate SHA-256 for %s", path), call. = FALSE)
    }
    return(strsplit(value, "[[:space:]]+")[[1L]][1L])
  }
  if (requireNamespace("openssl", quietly = TRUE)) {
    connection <- file(path, "rb")
    on.exit(close(connection))
    return(as.character(openssl::sha256(connection)))
  }
  stop("SHA-256 requires the sha256sum command or the R package openssl.", call. = FALSE)
}

expandPATTERN <- function(path) {
  pattern <- gsub("\\{chromosome\\}|\\{chr\\}|\\{CHR\\}", "*", path)
  pattern <- gsub("#", "*", pattern, fixed = TRUE)
  if (grepl("[*?]", pattern)) Sys.glob(pattern) else pattern
}

referenceFILES <- function(path, companion = "") {
  expanded <- expandPATTERN(path)
  if (length(companion) > 0L && nzchar(companion)) {
    expanded <- c(expanded, expandPATTERN(companion))
  }
  unique(unlist(lapply(expanded, function(value) {
    if (dir.exists(value)) {
      list.files(value, recursive = TRUE, full.names = TRUE, all.files = FALSE, no.. = TRUE)
    } else {
      value
    }
  }), use.names = FALSE))
}

resourceDIGEST <- function(files, base) {
  files <- sort(normalizePath(files, winslash = "/", mustWork = TRUE))
  if (length(files) == 1L) return(tolower(sha256FILE(files)))
  base <- normalizePath(base, winslash = "/", mustWork = TRUE)
  relative <- ifelse(startsWith(files, paste0(base, "/")), substring(files, nchar(base) + 2L), basename(files))
  inventory <- paste(vapply(files, sha256FILE, character(1L)), relative, sep = "  ")
  if (!requireNamespace("openssl", quietly = TRUE)) {
    stop("Composite reference checksums require the R package openssl.", call. = FALSE)
  }
  as.character(openssl::sha256(charToRaw(paste0(paste(inventory, collapse = "\n"), "\n"))))
}

checkTARGET <- function(path, format, label) {
  hasPATTERN <- grepl("\\{chromosome\\}|\\{chr\\}|\\{CHR\\}|#|[*?]", path)
  pattern <- gsub("\\{chromosome\\}|\\{chr\\}|\\{CHR\\}", "*", path)
  pattern <- gsub("#", "*", pattern, fixed = TRUE)
  if (hasPATTERN && format == "pgen" && !grepl("\\.pgen$", pattern, ignore.case = TRUE)) pattern <- paste0(pattern, ".pgen")
  if (hasPATTERN && format == "bed" && !grepl("\\.bed$", pattern, ignore.case = TRUE)) pattern <- paste0(pattern, ".bed")
  expanded <- if (hasPATTERN) Sys.glob(pattern) else path
  if (length(expanded) == 0L) {
    stop(sprintf("%s pattern matched no files: %s", label, path), call. = FALSE)
  }
  if (format == "pgen") {
    prefix <- ifelse(grepl("\\.pgen$", expanded, ignore.case = TRUE), sub("\\.pgen$", "", expanded, ignore.case = TRUE), expanded)
    companion <- as.vector(outer(prefix, c(".pgen", ".pvar", ".psam"), paste0))
    checkPATH(companion, paste(label, "PGEN companion"))
  } else if (format == "bed") {
    prefix <- ifelse(grepl("\\.bed$", expanded, ignore.case = TRUE), sub("\\.bed$", "", expanded, ignore.case = TRUE), expanded)
    companion <- as.vector(outer(prefix, c(".bed", ".bim", ".fam"), paste0))
    checkPATH(companion, paste(label, "BED companion"))
  } else if (format == "ped") {
    if (hasPATTERN) stop("PED/MAP input cannot use a chromosome pattern.", call. = FALSE)
    prefix <- ifelse(grepl("\\.ped$", expanded, ignore.case = TRUE), sub("\\.ped$", "", expanded, ignore.case = TRUE), expanded)
    companion <- as.vector(outer(prefix, c(".ped", ".map"), paste0))
    checkPATH(companion, paste(label, "PED/MAP companion"))
  } else {
    checkPATH(expanded, label)
  }
}

readHEADER <- function(path, required = character()) {
  connection <- if (grepl("\\.gz$", path, ignore.case = TRUE)) gzfile(path, "rt") else file(path, "rt")
  on.exit(close(connection))
  candidate <- readLines(connection, n = 500L, warn = FALSE)
  candidate <- candidate[nzchar(trimws(candidate))]
  if (length(candidate) == 0L) {
    stop(sprintf("GWAS file has no readable header: %s", path), call. = FALSE)
  }
  splitHEADER <- function(value) {
    value <- sub("^#", "", value)
    if (grepl("\t", value, fixed = TRUE)) {
      strsplit(value, "\t", fixed = TRUE)[[1L]]
    } else {
      strsplit(trimws(value), "[[:space:]]+")[[1L]]
    }
  }
  header <- lapply(candidate, splitHEADER)
  match <- which(vapply(header, function(value) all(required %in% value), logical(1L)))
  if (length(match) == 0L) {
    stop(
      sprintf(
        "GWAS file has no header containing the declared columns (%s): %s",
        paste(required, collapse = ", "),
        path
      ),
      call. = FALSE
    )
  }
  header[[match[1L]]]
}

readPHENOTYPE <- function(path) {
  firstLINE <- readLines(path, n = 1L, warn = FALSE)
  separator <- if (length(firstLINE) == 1L && grepl("\t", firstLINE, fixed = TRUE)) "\t" else ","
  utils::read.table(
    path,
    header = TRUE,
    sep = separator,
    quote = "\"",
    comment.char = "",
    check.names = FALSE,
    stringsAsFactors = FALSE,
    na.strings = c("", "NA")
  )
}

firstMATCH <- function(path) {
  pattern <- gsub("\\{chromosome\\}|\\{chr\\}|\\{CHR\\}", "*", path)
  pattern <- gsub("#", "*", pattern, fixed = TRUE)
  matched <- if (grepl("[*?]", pattern)) Sys.glob(pattern) else path
  if (length(matched) == 0L) return(NA_character_)
  sort(matched)[[1L]]
}

targetSAMPLEIDS <- function(target, nativePATH) {
  output <- character()
  for (row in seq_len(nrow(target))) {
    format <- target$source_format[[row]]
    genotype <- firstMATCH(nativePATH[[row]])
    if (is.na(genotype)) next
    if (format == "pgen") {
      prefix <- sub("\\.pgen$", "", genotype, ignore.case = TRUE)
      samplePATH <- paste0(prefix, ".psam")
      if (file.exists(samplePATH)) {
        sample <- utils::read.table(samplePATH, header = TRUE, comment.char = "", check.names = FALSE, stringsAsFactors = FALSE)
        iidCOLUMN <- intersect(c("IID", "#IID"), names(sample))
        if (length(iidCOLUMN) == 1L) output <- c(output, as.character(sample[[iidCOLUMN]]))
      }
    } else if (format %in% c("bed", "ped")) {
      extension <- if (format == "bed") "\\.bed$" else "\\.ped$"
      prefix <- sub(extension, "", genotype, ignore.case = TRUE)
      samplePATH <- paste0(prefix, if (format == "bed") ".fam" else ".ped")
      if (file.exists(samplePATH)) {
        sample <- utils::read.table(samplePATH, header = FALSE, stringsAsFactors = FALSE)
        if (ncol(sample) >= 2L) output <- c(output, as.character(sample[[2L]]))
      }
    } else if (format == "bgen") {
      samplePATH <- if (nzchar(target$sample[[row]])) resolvePATH(target$sample[[row]], launchNATIVE) else NA_character_
      if (length(samplePATH) == 1L && !is.na(samplePATH) && file.exists(samplePATH)) {
        sample <- utils::read.table(samplePATH, header = TRUE, skip = 1L, stringsAsFactors = FALSE)
        iidCOLUMN <- intersect(c("ID_2", "IID", "ID"), names(sample))
        if (length(iidCOLUMN) > 0L) output <- c(output, as.character(sample[[iidCOLUMN[[1L]]]]))
      }
    } else if (format == "vcf") {
      connection <- if (grepl("\\.gz$", genotype, ignore.case = TRUE)) gzfile(genotype, "rt") else file(genotype, "rt")
      header <- character()
      repeat {
        value <- readLines(connection, n = 1L, warn = FALSE)
        if (length(value) == 0L || startsWith(value, "#CHROM")) {
          header <- value
          break
        }
      }
      close(connection)
      if (length(header) == 1L) {
        fields <- strsplit(header, "\t", fixed = TRUE)[[1L]]
        if (length(fields) > 9L) output <- c(output, fields[10L:length(fields)])
      }
    }
  }
  unique(output[!is.na(output) & nzchar(trimws(output))])
}

# Read the four manifests and the run settings.
target <- readTABLE(option[["target-manifest"]])
gwas <- readTABLE(option[["gwas-manifest"]])
reference <- readTABLE(option[["reference-manifest"]])
phenotypeMODEL <- readTABLE(option[["phenotype-models"]])
method <- strsplit(option[["methods"]], ",", fixed = TRUE)[[1L]]
launchDIR <- pipelinePATH(option[["launch-dir"]], getwd())
launchNATIVE <- resolvePATH(launchDIR, getwd())
if (!dir.exists(launchNATIVE)) {
  stop(sprintf("Launch directory does not exist: %s", launchDIR), call. = FALSE)
}
genomeBUILD <- option[["genome-build"]]
referenceBASEVALUE <- option[["reference-base"]]
if (is.null(referenceBASEVALUE) || length(referenceBASEVALUE) == 0L || !nzchar(referenceBASEVALUE)) {
  referenceBASEVALUE <- dirname(resolvePATH(option[["reference-manifest"]], getwd()))
}
referenceBASE <- pipelinePATH(referenceBASEVALUE, getwd())
referenceBASENATIVE <- resolvePATH(referenceBASE, getwd())
if (!dir.exists(referenceBASENATIVE)) {
  stop(sprintf("Reference base directory does not exist: %s", referenceBASE), call. = FALSE)
}

# Validate both raw and completed target checkpoints. Input format and scientific stage
# are independent so a BED trio can enter as raw, corrected, QC-completed, or imputed.
if (!"source_format" %in% names(target) && "format" %in% names(target)) {
  target$source_format <- target$format
}
if (!"format" %in% names(target) && "source_format" %in% names(target)) {
  target$format <- target$source_format
}
targetOPTIONAL <- c("sample", "keep", "dosage", "assay_manifest", "marker_map")
for (column in setdiff(targetOPTIONAL, names(target))) target[[column]] <- ""
targetREQUIRED <- c("cohort", "role", "source_format", "genotype", "build", "ancestry")
requireCOLUMN(target, targetREQUIRED, "Target manifest")
for (column in targetREQUIRED) {
  requireVALUE(target, column, "Target manifest")
}
target$source_format <- tolower(target$source_format)
target$format <- target$source_format
target$role <- tolower(target$role)
if (!"input_stage" %in% names(target)) target$input_stage <- "raw"
target$input_stage <- tolower(target$input_stage)
target$input_stage[target$input_stage == ""] <- "raw"
if (any(!target$source_format %in% c("pgen", "vcf", "bgen", "bed", "ped", "genomestudio"))) {
  stop("Target source_format must be pgen, vcf, bgen, bed, ped, or genomestudio.", call. = FALSE)
}
if (any(!target$role %in% c("target", "validation"))) {
  stop("Target role must be target or validation.", call. = FALSE)
}
if (any(!target$input_stage %in% c("raw", "corrected", "qc_completed", "imputed"))) {
  stop("Target input_stage must be raw, corrected, qc_completed, or imputed.", call. = FALSE)
}
if (any(target$build != genomeBUILD)) {
  stop("Every target row must use the selected genome build.", call. = FALSE)
}
target$genotype <- pipelinePATH(target$genotype, launchDIR)
if (any(grepl("[*?]", target$genotype))) {
  stop(
    "Target genotype paths must use {chr}, {CHR}, {chromosome}, or # instead of a raw * or ? glob.",
    call. = FALSE
  )
}
samplePRESENT <- !is.na(target$sample) & trimws(target$sample) != ""
keepPRESENT <- !is.na(target$keep) & trimws(target$keep) != ""
assayPRESENT <- !is.na(target$assay_manifest) & trimws(target$assay_manifest) != ""
markerMAPRESENT <- !is.na(target$marker_map) & trimws(target$marker_map) != ""
target$sample[!samplePRESENT] <- ""
target$keep[!keepPRESENT] <- ""
target$assay_manifest[!assayPRESENT] <- ""
target$marker_map[!markerMAPRESENT] <- ""
target$sample[samplePRESENT] <- pipelinePATH(target$sample[samplePRESENT], launchDIR)
target$keep[keepPRESENT] <- pipelinePATH(target$keep[keepPRESENT], launchDIR)
target$assay_manifest[assayPRESENT] <- pipelinePATH(target$assay_manifest[assayPRESENT], launchDIR)
target$marker_map[markerMAPRESENT] <- pipelinePATH(target$marker_map[markerMAPRESENT], launchDIR)
targetGENOTYPENATIVE <- resolvePATH(target$genotype, launchNATIVE)
targetSAMPLENATIVE <- resolvePATH(target$sample[target$sample != ""], launchNATIVE)
targetKEEPNATIVE <- resolvePATH(target$keep[target$keep != ""], launchNATIVE)
targetASSAYNATIVE <- resolvePATH(target$assay_manifest[target$assay_manifest != ""], launchNATIVE)
targetMARKERMAPNATIVE <- resolvePATH(target$marker_map[target$marker_map != ""], launchNATIVE)
for (row in seq_len(nrow(target))) {
  checkTARGET(targetGENOTYPENATIVE[row], target$source_format[row], sprintf("Target '%s' genotype", target$cohort[row]))
  if (target$source_format[row] == "genomestudio" && !nzchar(target$assay_manifest[row])) {
    stop(sprintf("Target '%s' requires assay_manifest for GenomeStudio input.", target$cohort[row]), call. = FALSE)
  }
}
checkPATH(targetSAMPLENATIVE, "Target sample input")
checkPATH(targetKEEPNATIVE, "Target keep input")
checkPATH(targetASSAYNATIVE, "Target assay manifest")
checkPATH(targetMARKERMAPNATIVE, "Target marker map")
if (anyDuplicated(target$cohort)) {
  stop("Target cohort values must be unique.", call. = FALSE)
}

# Validate the GWAS contract and each explicitly mapped source column.
gwasOPTIONAL <- c(
  "source_format", "freq_col", "case_freq_col", "control_freq_col",
  "case_n_col", "control_n_col", "n_col", "info_col", "info_min", "maf_min"
)
for (column in setdiff(gwasOPTIONAL, names(gwas))) gwas[[column]] <- ""
gwas$source_format[gwas$source_format == "" | is.na(gwas$source_format)] <- "auto"
gwas$info_min[is.na(gwas$info_min)] <- ""
gwas$maf_min[gwas$maf_min == "" | is.na(gwas$maf_min)] <- "0.01"
gwasREQUIRED <- c(
  "trait_id", "prs_name", "path", "build", "ancestry", "effect_type", "sample_size",
  "snp_col", "chr_col", "bp_col", "effect_allele_col", "other_allele_col", "beta_col",
  "se_col", "p_col", "freq_col", "n_col"
)
requireCOLUMN(gwas, gwasREQUIRED, "GWAS manifest")
for (column in setdiff(gwasREQUIRED, c("freq_col", "n_col"))) {
  requireVALUE(gwas, column, "GWAS manifest")
}
if (anyDuplicated(gwas$trait_id) || anyDuplicated(gwas$prs_name)) {
  stop("GWAS trait_id and prs_name values must each be unique.", call. = FALSE)
}
if (any(gwas$build != genomeBUILD)) {
  stop("Every GWAS row must use the selected genome build.", call. = FALSE)
}
gwas$path <- pipelinePATH(gwas$path, launchDIR)
gwasPATHNATIVE <- resolvePATH(gwas$path, launchNATIVE)
checkPATH(gwasPATHNATIVE, "GWAS input")
for (row in seq_len(nrow(gwas))) {
  declaredCOLUMN <- unlist(gwas[row, c(
    "snp_col", "chr_col", "bp_col", "effect_allele_col", "other_allele_col", "beta_col", "se_col", "p_col"
  )], use.names = FALSE)
  optionalCOLUMN <- unlist(gwas[row, c(
    "freq_col", "case_freq_col", "control_freq_col", "case_n_col", "control_n_col",
    "n_col", "info_col"
  )], use.names = FALSE)
  optionalCOLUMN <- optionalCOLUMN[!is.na(optionalCOLUMN) & trimws(optionalCOLUMN) != ""]
  missingCOLUMN <- setdiff(c(declaredCOLUMN, optionalCOLUMN), readHEADER(gwasPATHNATIVE[row], declaredCOLUMN))
  if (length(missingCOLUMN) > 0L) {
    stop(
      sprintf("GWAS '%s' is missing declared column(s): %s", gwas$trait_id[row], paste(missingCOLUMN, collapse = ", ")),
      call. = FALSE
    )
  }
}

# Validate immutable reference sources. Tool-ready PLINK and extracted SBayesRC
# resources are generated by the workflow and must not be supplied as run inputs.
referenceREQUIRED <- c(
  "bundle_id", "bundle_version", "reference_id", "reference_type", "path",
  "build", "ancestry", "version", "checksum", "source_format", "reference_stage"
)
requireCOLUMN(reference, referenceREQUIRED, "Reference manifest")
for (column in setdiff(referenceREQUIRED, "checksum")) {
  requireVALUE(reference, column, "Reference manifest")
}
if (anyDuplicated(reference$reference_id)) {
  stop("Reference identifiers must be unique.", call. = FALSE)
}
if (length(unique(reference$bundle_id)) != 1L || length(unique(reference$bundle_version)) != 1L) {
  stop("Every reference row must use one bundle_id and one bundle_version.", call. = FALSE)
}
if (any(reference$build != genomeBUILD)) {
  stop("Every reference row must use the selected genome build.", call. = FALSE)
}
reference$reference_type <- tolower(reference$reference_type)
reference$source_format <- tolower(reference$source_format)
reference$reference_stage <- tolower(reference$reference_stage)
if (any(reference$reference_stage != "source")) {
  stop(
    "Every reference must have reference_stage=source; prepared references are workflow outputs, not pipeline inputs.",
    call. = FALSE
  )
}
allowedREFERENCE <- c(
  "dbsnp", "reference_fasta", "genetic_map", "imputation_panel",
  "population_panel", "related_samples", "sbayesrc_ld_source", "annotation_source",
  "beagle_jar", "unbref3_jar"
)
unsupportedREFERENCE <- setdiff(unique(reference$reference_type), allowedREFERENCE)
if (length(unsupportedREFERENCE) > 0L) {
  stop(
    sprintf("Unsupported source reference type(s): %s", paste(unsupportedREFERENCE, collapse = ", ")),
    call. = FALSE
  )
}
if (anyDuplicated(reference$reference_type)) {
  stop("Reference source roles must be unique within a bundle.", call. = FALSE)
}
allowedSOURCEFORMAT <- list(
  dbsnp = c("directory", "vcf", "vcf.gz"),
  reference_fasta = "fasta",
  genetic_map = c("zip", "directory"),
  imputation_panel = c("bref3", "bref3_directory", "vcf", "vcf_directory"),
  population_panel = c("tsv", "panel"),
  related_samples = c("tsv", "txt"),
  sbayesrc_ld_source = "zip",
  annotation_source = "zip",
  beagle_jar = "jar",
  unbref3_jar = "jar"
)
invalidSOURCEFORMAT <- vapply(seq_len(nrow(reference)), function(row) {
  !reference$source_format[row] %in% allowedSOURCEFORMAT[[reference$reference_type[row]]]
}, logical(1L))
if (any(invalidSOURCEFORMAT)) {
  row <- which(invalidSOURCEFORMAT)[1L]
  stop(
    sprintf(
      "Reference '%s' has source_format=%s; reference_type=%s accepts: %s.",
      reference$reference_id[row], reference$source_format[row], reference$reference_type[row],
      paste(allowedSOURCEFORMAT[[reference$reference_type[row]]], collapse = ", ")
    ),
    call. = FALSE
  )
}
if (!"companion" %in% names(reference)) reference$companion <- ""
reference$companion[is.na(reference$companion)] <- ""
reference$path <- pipelinePATH(reference$path, referenceBASE)
companionPRESENT <- nzchar(trimws(reference$companion))
reference$companion[companionPRESENT] <- pipelinePATH(reference$companion[companionPRESENT], referenceBASE)
referencePATHNATIVE <- resolvePATH(reference$path, referenceBASENATIVE)
referenceCOMPANIONNATIVE <- rep("", nrow(reference))
referenceCOMPANIONNATIVE[companionPRESENT] <- resolvePATH(
  reference$companion[companionPRESENT],
  referenceBASENATIVE
)
referenceFILELIST <- vector("list", nrow(reference))
verifiedREFERENCE <- logical(nrow(reference))
for (row in seq_len(nrow(reference))) {
  checkPATH(expandPATTERN(referencePATHNATIVE[row]), sprintf("Reference '%s'", reference$reference_id[row]))
  if (companionPRESENT[row]) {
    checkPATH(expandPATTERN(referenceCOMPANIONNATIVE[row]), sprintf("Reference '%s' companion", reference$reference_id[row]))
  }
  referenceFILELIST[[row]] <- referenceFILES(referencePATHNATIVE[row], referenceCOMPANIONNATIVE[row])
  checkPATH(referenceFILELIST[[row]], sprintf("Reference '%s' content", reference$reference_id[row]))
  fileNAME <- basename(referenceFILELIST[[row]])
  referenceTYPE <- reference$reference_type[row]
  hasCONTENT <- switch(
    referenceTYPE,
    dbsnp = any(grepl("\\.vcf(\\.gz)?$|\\.gz$", fileNAME, ignore.case = TRUE)) &&
      any(grepl("assembly[_-]?report", fileNAME, ignore.case = TRUE)),
    reference_fasta = any(grepl("\\.(fa|fasta|fna)(\\.gz)?$", fileNAME, ignore.case = TRUE)) &&
      any(grepl("\\.fai$", fileNAME, ignore.case = TRUE)),
    genetic_map = reference$source_format[row] == "zip" || any(grepl("\\.map(\\.gz)?$", fileNAME, ignore.case = TRUE)),
    imputation_panel = any(grepl("\\.(bref3|vcf\\.gz|bcf)$", fileNAME, ignore.case = TRUE)),
    population_panel = length(referenceFILELIST[[row]]) >= 1L,
    related_samples = length(referenceFILELIST[[row]]) >= 1L,
    sbayesrc_ld_source = any(grepl("\\.zip$", fileNAME, ignore.case = TRUE)),
    annotation_source = any(grepl("\\.zip$", fileNAME, ignore.case = TRUE)),
    beagle_jar = any(grepl("\\.jar$", fileNAME, ignore.case = TRUE)),
    unbref3_jar = any(grepl("\\.jar$", fileNAME, ignore.case = TRUE)),
    FALSE
  )
  if (!hasCONTENT) {
    stop(
      sprintf("Reference '%s' does not contain the files expected for reference_type=%s.", reference$reference_id[row], referenceTYPE),
      call. = FALSE
    )
  }
  suppliedCHECKSUM <- !is.na(reference$checksum[row]) &&
    nzchar(trimws(reference$checksum[row])) &&
    toupper(trimws(reference$checksum[row])) != "NA"
  if (suppliedCHECKSUM) {
    declaredCHECKSUM <- tolower(trimws(reference$checksum[row]))
    if (!grepl("^[0-9a-f]{64}$", declaredCHECKSUM)) {
      stop(sprintf("Reference '%s' checksum must be a 64-character SHA-256 value.", reference$reference_id[row]), call. = FALSE)
    }
    digestBASE <- if (dir.exists(referencePATHNATIVE[row])) referencePATHNATIVE[row] else dirname(referencePATHNATIVE[row])
    observedCHECKSUM <- resourceDIGEST(referenceFILELIST[[row]], digestBASE)
    if (!identical(declaredCHECKSUM, observedCHECKSUM)) {
      stop(sprintf("Reference '%s' failed its declared SHA-256 check.", reference$reference_id[row]), call. = FALSE)
    }
    reference$checksum[row] <- declaredCHECKSUM
    verifiedREFERENCE[row] <- TRUE
  }
}
requiredREFERENCE <- c(
  "dbsnp", "reference_fasta", "imputation_panel", "population_panel",
  "related_samples", "unbref3_jar"
)
if ("sbayesrc" %in% method) {
  requiredREFERENCE <- c(requiredREFERENCE, "sbayesrc_ld_source", "annotation_source")
  directFREQUENCY <- !is.na(gwas$freq_col) & trimws(gwas$freq_col) != ""
  weightedFREQUENCY <- apply(
    gwas[, c("case_freq_col", "control_freq_col", "case_n_col", "control_n_col"), drop = FALSE],
    1L,
    function(value) all(!is.na(value) & trimws(value) != "")
  )
  missingFREQUENCY <- !(directFREQUENCY | weightedFREQUENCY)
  if (any(missingFREQUENCY)) {
    stop(
      "SBayesRC requires either freq_col or all case/control frequency and sample-count columns for every GWAS.",
      call. = FALSE
    )
  }
}
if (tolower(option[["target-imputation"]]) == "true") {
  requiredREFERENCE <- c(requiredREFERENCE, "imputation_panel", "genetic_map", "beagle_jar")
}
missingREFERENCE <- setdiff(unique(requiredREFERENCE), reference$reference_type)
if (length(missingREFERENCE) > 0L) {
  stop(sprintf("Selected methods require reference type(s): %s", paste(missingREFERENCE, collapse = ", ")), call. = FALSE)
}

# Validate phenotype model declarations only when phenotype analysis is requested.
modelREQUIRED <- c(
  "model_id", "outcome", "prs_name", "family", "covariates", "participant_id",
  "group_id", "expected_direction", "primary"
)
if (nrow(phenotypeMODEL) > 0L) {
  for (column in setdiff(c("control_value", "case_value", "subset"), names(phenotypeMODEL))) {
    phenotypeMODEL[[column]] <- ""
  }
  requireCOLUMN(phenotypeMODEL, modelREQUIRED, "Phenotype-model manifest")
  for (column in c("model_id", "outcome", "prs_name", "family", "primary")) {
    requireVALUE(phenotypeMODEL, column, "Phenotype-model manifest")
  }
  if (anyDuplicated(phenotypeMODEL$model_id)) {
    stop("Phenotype model identifiers must be unique.", call. = FALSE)
  }
  phenotypePATH <- resolvePATH(option[["phenotype-file"]], getwd())
  if (!file.exists(phenotypePATH)) {
    phenotypePATH <- resolvePATH(option[["phenotype-file"]], launchNATIVE)
  }
  checkPATH(phenotypePATH, "Phenotype input")
  phenotypeDATA <- readPHENOTYPE(phenotypePATH)
  phenotypeHEADER <- names(phenotypeDATA)

  declaredPARTICIPANT <- unique(trimws(phenotypeMODEL$participant_id[!is.na(phenotypeMODEL$participant_id)]))
  declaredPARTICIPANT <- declaredPARTICIPANT[nzchar(declaredPARTICIPANT)]
  if (length(declaredPARTICIPANT) > 1L) {
    stop("All phenotype models in one run must use one participant-ID column.", call. = FALSE)
  }
  if (length(declaredPARTICIPANT) == 0L) {
    targetID <- targetSAMPLEIDS(target, targetGENOTYPENATIVE)
    if (length(targetID) == 0L) {
      stop("Could not read target sample IDs; provide --participant_id explicitly.", call. = FALSE)
    }
    minimumMATCH <- min(length(targetID), max(2L, min(10L, ceiling(length(targetID) * 0.8))))
    candidate <- phenotypeHEADER[vapply(phenotypeDATA, function(value) {
      value <- as.character(value)
      value <- value[!is.na(value) & nzchar(trimws(value))]
      overlap <- length(intersect(unique(value), targetID))
      overlap >= minimumMATCH && overlap / length(targetID) >= 0.8 && !anyDuplicated(value)
    }, logical(1L))]
    if (length(candidate) != 1L) {
      detail <- if (length(candidate) == 0L) "none matched" else paste(candidate, collapse = ", ")
      stop(
        sprintf("Participant-ID inference requires exactly one phenotype column matching target IDs; %s. Provide --participant_id.", detail),
        call. = FALSE
      )
    }
    phenotypeMODEL$participant_id <- candidate[[1L]]
    declaredPARTICIPANT <- candidate[[1L]]
  } else {
    phenotypeMODEL$participant_id[is.na(phenotypeMODEL$participant_id) | trimws(phenotypeMODEL$participant_id) == ""] <- declaredPARTICIPANT[[1L]]
  }

  modelCOLUMN <- unique(c(
    phenotypeMODEL$outcome,
    phenotypeMODEL$participant_id,
    phenotypeMODEL$group_id[!is.na(phenotypeMODEL$group_id) & phenotypeMODEL$group_id != ""],
    trimws(unlist(strsplit(phenotypeMODEL$covariates, ",", fixed = TRUE)))
  ))
  modelCOLUMN <- modelCOLUMN[!is.na(modelCOLUMN) & modelCOLUMN != ""]
  missingMODEL <- setdiff(modelCOLUMN, phenotypeHEADER)
  if (length(missingMODEL) > 0L) {
    stop(sprintf("Phenotype file is missing model column(s): %s", paste(missingMODEL, collapse = ", ")), call. = FALSE)
  }
  missingPRS <- setdiff(unique(phenotypeMODEL$prs_name), gwas$prs_name)
  if (length(missingPRS) > 0L) {
    stop(sprintf("Phenotype models refer to unknown PRS name(s): %s", paste(missingPRS, collapse = ", ")), call. = FALSE)
  }

  participantVALUE <- phenotypeDATA[[declaredPARTICIPANT[[1L]]]]
  if (anyNA(participantVALUE) || any(!nzchar(trimws(as.character(participantVALUE))))) {
    stop(sprintf("Phenotype participant-ID column '%s' contains missing values.", declaredPARTICIPANT[[1L]]), call. = FALSE)
  }
  if (anyDuplicated(as.character(participantVALUE))) {
    stop(sprintf("Phenotype participant-ID column '%s' contains duplicates.", declaredPARTICIPANT[[1L]]), call. = FALSE)
  }

  allowedFAMILY <- c("gaussian", "binomial")
  phenotypeMODEL$family <- tolower(trimws(phenotypeMODEL$family))
  if (any(!phenotypeMODEL$family %in% allowedFAMILY)) {
    stop("Phenotype model family must be gaussian or binomial after input resolution.", call. = FALSE)
  }
  for (row in seq_len(nrow(phenotypeMODEL))) {
    modelID <- phenotypeMODEL$model_id[[row]]
    outcome <- phenotypeMODEL$outcome[[row]]
    outcomeVALUE <- phenotypeDATA[[outcome]]
    observed <- outcomeVALUE[!is.na(outcomeVALUE)]
    if (length(unique(observed)) < 2L) {
      stop(sprintf("Outcome '%s' for model '%s' has fewer than two observed values.", outcome, modelID), call. = FALSE)
    }
    if (phenotypeMODEL$family[[row]] == "gaussian" && !is.numeric(outcomeVALUE)) {
      stop(sprintf("Gaussian outcome '%s' for model '%s' must be numeric.", outcome, modelID), call. = FALSE)
    }
    if (phenotypeMODEL$family[[row]] == "binomial") {
      control <- trimws(as.character(phenotypeMODEL$control_value[[row]]))
      case <- trimws(as.character(phenotypeMODEL$case_value[[row]]))
      if (is.na(control)) control <- ""
      if (is.na(case)) case <- ""
      if (xor(nzchar(control), nzchar(case))) {
        stop(sprintf("Binomial model '%s' must declare both control_value and case_value, or neither.", modelID), call. = FALSE)
      }
      observedTEXT <- unique(as.character(observed))
      if (!nzchar(control) && !nzchar(case)) {
        if (!setequal(observedTEXT, c("0", "1"))) {
          stop(sprintf("Binomial model '%s' requires 0/1 outcome coding or explicit control_value and case_value.", modelID), call. = FALSE)
        }
        phenotypeMODEL$control_value[[row]] <- "0"
        phenotypeMODEL$case_value[[row]] <- "1"
      } else if (identical(control, case) || !all(c(control, case) %in% observedTEXT) || length(observedTEXT) != 2L) {
        stop(sprintf("Binomial model '%s' must map exactly two observed outcome values to distinct control/case values.", modelID), call. = FALSE)
      }
    }
    covariate <- trimws(strsplit(phenotypeMODEL$covariates[[row]], ",", fixed = TRUE)[[1L]])
    covariate <- covariate[nzchar(covariate)]
    constantCOVARIATE <- covariate[vapply(covariate, function(column) {
      length(unique(phenotypeDATA[[column]][!is.na(phenotypeDATA[[column]])])) < 2L
    }, logical(1L))]
    if (length(constantCOVARIATE) > 0L) {
      stop(sprintf("Model '%s' has constant or empty covariate(s): %s", modelID, paste(constantCOVARIATE, collapse = ", ")), call. = FALSE)
    }
    group <- trimws(as.character(phenotypeMODEL$group_id[[row]]))
    if (is.na(group)) group <- ""
    if (nzchar(group) && length(unique(phenotypeDATA[[group]][!is.na(phenotypeDATA[[group]])])) < 2L) {
      stop(sprintf("Grouping column '%s' for model '%s' has fewer than two groups.", group, modelID), call. = FALSE)
    }
  }
}

# Save resolved manifests and a concise input-check record.
writeTABLE(target, "targets.tsv")
writeTABLE(gwas, "gwas.tsv")
writeTABLE(reference, "references.tsv")
writeTABLE(phenotypeMODEL, "models.tsv")

yamlVALUE <- c(
  sprintf("run_name: '%s'", option[["run-name"]]),
  sprintf("genome_build: '%s'", genomeBUILD),
  "methods:",
  paste0("  - '", method, "'"),
  sprintf("seed: %s", option[["seed"]]),
  sprintf("report_enabled: %s", tolower(option[["report-enabled"]])),
  sprintf("target_imputation: %s", tolower(option[["target-imputation"]])),
  sprintf("launch_directory: '%s'", launchDIR),
  sprintf("reference_bundle_id: '%s'", unique(reference$bundle_id)),
  sprintf("reference_bundle_version: '%s'", unique(reference$bundle_version)),
  sprintf("reference_base: '%s'", referenceBASE)
)
writeLines(yamlVALUE, "run_settings.yml")

checksumTABLE <- unique(data.frame(
  display = c(
    basename(option[["target-manifest"]]),
    basename(option[["gwas-manifest"]]),
    basename(option[["reference-manifest"]]),
    gwas$path,
    reference$path,
    reference$companion[companionPRESENT]
  ),
  native = c(
    resolvePATH(option[["target-manifest"]], getwd()),
    resolvePATH(option[["gwas-manifest"]], getwd()),
    resolvePATH(option[["reference-manifest"]], getwd()),
    gwasPATHNATIVE,
    referencePATHNATIVE,
    referenceCOMPANIONNATIVE[companionPRESENT]
  ),
  stringsAsFactors = FALSE
))
checksumTABLE <- checksumTABLE[file.exists(checksumTABLE$native) & !dir.exists(checksumTABLE$native), , drop = FALSE]
checksumDISPLAY <- checksumTABLE$display
checksumPATH <- checksumTABLE$native
checksumCOMMAND <- vapply(
  checksumPATH,
  sha256FILE,
  character(1L)
)
writeTABLE(
  data.frame(path = checksumDISPLAY, algorithm = "SHA-256", checksum = checksumCOMMAND, stringsAsFactors = FALSE),
  "input_checksums.tsv"
)

inputCHECK <- data.frame(
  check = c("target_cohorts", "gwas_traits", "reference_sources", "reference_stage", "reference_bundle", "reference_checksums", "phenotype_models", "genome_build", "methods"),
  value = c(
    nrow(target), nrow(gwas), nrow(reference),
    "source",
    paste0(unique(reference$bundle_id), ":", unique(reference$bundle_version)),
    sprintf("%d supplied and verified; %d not supplied", sum(verifiedREFERENCE), sum(!verifiedREFERENCE)),
    nrow(phenotypeMODEL), genomeBUILD, paste(method, collapse = ",")
  ),
  status = "PASS",
  stringsAsFactors = FALSE
)
writeTABLE(inputCHECK, "input_checks.tsv")
