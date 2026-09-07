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

# Define compact helpers for generated tabular records and clear failures.
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

# Nextflow may copy declared inputs into the task directory instead of bind-mounting
# their original parent directories. Keep the original record paths for provenance,
# but use the staged copies for validation inside Docker or Apptainer.
stagedMATCH <- new.env(parent = emptyenv())
stagedPATH <- function(path, root) {
  path <- as.character(path)
  if (length(path) == 0L || is.null(root) || length(root) == 0L || !nzchar(root) || !dir.exists(root)) {
    return(path)
  }
  vapply(path, function(value) {
    if (file.exists(value) || dir.exists(value)) return(value)
    patternVALUE <- grepl("\\{chromosome\\}|\\{chr\\}|\\{CHR\\}|#|[*?]", value)
    fileNAME <- basename(value)
    if (patternVALUE) {
      fileNAME <- gsub("\\{chromosome\\}|\\{chr\\}|\\{CHR\\}", "*", fileNAME)
      fileNAME <- gsub("#", "*", fileNAME, fixed = TRUE)
      return(file.path(root, "*", fileNAME))
    }
    candidate <- list.files(
      root,
      recursive = TRUE,
      full.names = TRUE,
      all.files = FALSE,
      include.dirs = TRUE,
      no.. = TRUE
    )
    candidate <- sort(candidate[basename(candidate) == fileNAME])
    if (length(candidate) == 0L) return(value)
    sourceKEY <- paste(root, value, sep = "\r")
    if (exists(sourceKEY, envir = stagedMATCH, inherits = FALSE)) {
      return(get(sourceKEY, envir = stagedMATCH, inherits = FALSE))
    }
    baseKEY <- paste(root, fileNAME, sep = "\r")
    occurrence <- if (exists(baseKEY, envir = stagedMATCH, inherits = FALSE)) {
      get(baseKEY, envir = stagedMATCH, inherits = FALSE) + 1L
    } else {
      1L
    }
    occurrence <- min(occurrence, length(candidate))
    assign(baseKEY, occurrence, envir = stagedMATCH)
    assign(sourceKEY, candidate[[occurrence]], envir = stagedMATCH)
    candidate[[occurrence]]
  }, character(1L), USE.NAMES = FALSE)
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

relativeFILENAMES <- function(files, root) {
  if (!dir.exists(root)) return(basename(files))
  normalROOT <- normalizePath(root, winslash = "/", mustWork = TRUE)
  normalFILE <- normalizePath(files, winslash = "/", mustWork = TRUE)
  prefix <- paste0(normalROOT, "/")
  ifelse(startsWith(normalFILE, prefix), substring(normalFILE, nchar(prefix) + 1L), basename(normalFILE))
}

archiveNAMES <- function(path, label) {
  value <- try(utils::unzip(path, list = TRUE), silent = TRUE)
  if (inherits(value, "try-error")) {
    stop(sprintf("%s is not a readable ZIP archive: %s", label, path), call. = FALSE)
  }
  value$Name
}

referenceCHROMOSOMES <- function(fileNAME, resource) {
  fileNAME <- sub("^\\./", "", gsub("\\\\", "/", fileNAME))
  if (resource == "genetic_map") {
    supported <- grepl(
      "^((maps/)?plink\\.chr[0-9]+\\.GRCh37\\.map|chr[0-9]+\\.map)$",
      fileNAME
    )
    chromosome <- basename(fileNAME[supported])
    chromosome <- sub("^plink\\.chr", "", chromosome)
    chromosome <- sub("^chr", "", chromosome)
    chromosome <- sub("\\.GRCh37\\.map$", "", chromosome)
    chromosome <- sub("\\.map$", "", chromosome)
  } else {
    supported <- grepl(
      "^(panel/)?chr[0-9]+((\\..+)?\\.bref3|\\.vcf(\\.gz)?|\\.bcf)$",
      fileNAME
    )
    chromosome <- sub("^chr([0-9]+).*$", "\\1", basename(fileNAME[supported]))
  }
  chromosome <- suppressWarnings(as.integer(chromosome))
  chromosome[!is.na(chromosome) & chromosome >= 1L & chromosome <= 22L]
}

checkCHROMOSOMES <- function(chromosome, referenceID, label, requireCOMPLETE = FALSE) {
  if (length(chromosome) == 0L) {
    stop(sprintf("Reference '%s' contains no supported autosomal %s files.", referenceID, label), call. = FALSE)
  }
  duplicateCHROMOSOME <- sort(unique(chromosome[duplicated(chromosome)]))
  if (length(duplicateCHROMOSOME) > 0L) {
    stop(
      sprintf(
        "Reference '%s' contains more than one %s file for chromosome(s): %s. Keep one supported file per chromosome.",
        referenceID,
        label,
        paste(duplicateCHROMOSOME, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  if (requireCOMPLETE) {
    missingCHROMOSOME <- setdiff(seq_len(22L), chromosome)
    if (length(missingCHROMOSOME) > 0L) {
      stop(
        sprintf(
          "Reference '%s' is missing %s file(s) for chromosome(s): %s. Provide chromosomes 1 to 22.",
          referenceID,
          label,
          paste(missingCHROMOSOME, collapse = ", ")
        ),
        call. = FALSE
      )
    }
  }
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

checkTARGET <- function(path, format, label, stagedROOT = "") {
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
    checkPATH(stagedPATH(companion, stagedROOT), paste(label, "PGEN companion"))
  } else if (format == "bed") {
    prefix <- ifelse(grepl("\\.bed$", expanded, ignore.case = TRUE), sub("\\.bed$", "", expanded, ignore.case = TRUE), expanded)
    companion <- as.vector(outer(prefix, c(".bed", ".bim", ".fam"), paste0))
    checkPATH(stagedPATH(companion, stagedROOT), paste(label, "BED companion"))
  } else if (format == "ped") {
    if (hasPATTERN) stop("PED/MAP input cannot use a chromosome pattern.", call. = FALSE)
    prefix <- ifelse(grepl("\\.ped$", expanded, ignore.case = TRUE), sub("\\.ped$", "", expanded, ignore.case = TRUE), expanded)
    companion <- as.vector(outer(prefix, c(".ped", ".map"), paste0))
    checkPATH(stagedPATH(companion, stagedROOT), paste(label, "PED/MAP companion"))
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

targetSAMPLEIDS <- function(target, nativePATH, stagedROOT = "") {
  output <- character()
  for (row in seq_len(nrow(target))) {
    format <- target$source_format[[row]]
    genotype <- firstMATCH(nativePATH[[row]])
    if (is.na(genotype)) next
    if (format == "pgen") {
      prefix <- sub("\\.pgen$", "", genotype, ignore.case = TRUE)
      samplePATH <- stagedPATH(paste0(prefix, ".psam"), stagedROOT)
      if (file.exists(samplePATH)) {
        sample <- utils::read.table(samplePATH, header = TRUE, comment.char = "", check.names = FALSE, stringsAsFactors = FALSE)
        iidCOLUMN <- intersect(c("IID", "#IID"), names(sample))
        if (length(iidCOLUMN) == 1L) output <- c(output, as.character(sample[[iidCOLUMN]]))
      }
    } else if (format %in% c("bed", "ped")) {
      extension <- if (format == "bed") "\\.bed$" else "\\.ped$"
      prefix <- sub(extension, "", genotype, ignore.case = TRUE)
      samplePATH <- stagedPATH(paste0(prefix, if (format == "bed") ".fam" else ".ped"), stagedROOT)
      if (file.exists(samplePATH)) {
        sample <- utils::read.table(samplePATH, header = FALSE, stringsAsFactors = FALSE)
        if (ncol(sample) >= 2L) output <- c(output, as.character(sample[[2L]]))
      }
    } else if (format == "bgen") {
      samplePATH <- if (nzchar(target$sample[[row]])) {
        stagedPATH(resolvePATH(target$sample[[row]], launchNATIVE), stagedROOT)
      } else {
        NA_character_
      }
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

# Read the four generated record tables and the run settings.
target <- readTABLE(option[["target-manifest"]])
gwas <- readTABLE(option[["gwas-manifest"]])
reference <- readTABLE(option[["reference-manifest"]])
phenotypeMODEL <- readTABLE(option[["phenotype-models"]])
runPLAN <- readTABLE(option[["run-plan"]])
if (!all(c("setting", "value") %in% names(runPLAN)) || anyDuplicated(runPLAN$setting)) {
  stop("run_plan.tsv must contain unique setting and value columns.", call. = FALSE)
}
planVALUE <- stats::setNames(as.list(runPLAN$value), runPLAN$setting)
runPRS <- identical(tolower(planVALUE[["run_prs"]]), "true")
runPHENOTYPE <- identical(tolower(planVALUE[["run_phenotype"]]), "true")
referenceONLY <- identical(tolower(planVALUE[["reference_only"]]), "true")
requiredREFERENCE <- unique(trimws(strsplit(planVALUE[["required_reference_roles"]], ",", fixed = TRUE)[[1L]]))
if (length(requiredREFERENCE) == 0L || any(!nzchar(requiredREFERENCE))) {
  stop("run_plan.tsv has no valid required_reference_roles.", call. = FALSE)
}
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
targetSTAGEROOT <- resolvePATH(option[["target-assets"]], getwd())
gwasSTAGEROOT <- resolvePATH(option[["gwas-assets"]], getwd())
referenceSTAGEROOT <- resolvePATH(option[["reference-assets"]], getwd())

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
requireCOLUMN(target, targetREQUIRED, "Generated target records")
for (column in targetREQUIRED) {
  requireVALUE(target, column, "Generated target records")
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
targetGENOTYPENATIVE <- stagedPATH(resolvePATH(target$genotype, launchNATIVE), targetSTAGEROOT)
targetSAMPLENATIVE <- stagedPATH(resolvePATH(target$sample[target$sample != ""], launchNATIVE), targetSTAGEROOT)
targetKEEPNATIVE <- stagedPATH(resolvePATH(target$keep[target$keep != ""], launchNATIVE), targetSTAGEROOT)
targetASSAYNATIVE <- stagedPATH(resolvePATH(target$assay_manifest[target$assay_manifest != ""], launchNATIVE), targetSTAGEROOT)
targetMARKERMAPNATIVE <- stagedPATH(resolvePATH(target$marker_map[target$marker_map != ""], launchNATIVE), targetSTAGEROOT)
for (row in seq_len(nrow(target))) {
  checkTARGET(
    targetGENOTYPENATIVE[row],
    target$source_format[row],
    sprintf("Target '%s' genotype", target$cohort[row]),
    targetSTAGEROOT
  )
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

# Validate the GWAS contract and each explicitly mapped source column when the selected
# workflow stage includes PRS generation.
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
requireCOLUMN(gwas, gwasREQUIRED, "Generated GWAS records")
gwasPATHNATIVE <- character()
if (runPRS) {
  if (nrow(gwas) == 0L) stop("The selected PRS stage requires at least one GWAS input.", call. = FALSE)
  for (column in setdiff(gwasREQUIRED, c("freq_col", "n_col"))) {
    requireVALUE(gwas, column, "Generated GWAS records")
  }
  if (anyDuplicated(gwas$trait_id) || anyDuplicated(gwas$prs_name)) {
    stop("GWAS trait_id and prs_name values must each be unique.", call. = FALSE)
  }
  if (any(gwas$build != genomeBUILD)) {
    stop("Every GWAS row must use the selected genome build.", call. = FALSE)
  }
  gwas$path <- pipelinePATH(gwas$path, launchDIR)
  gwasPATHNATIVE <- stagedPATH(resolvePATH(gwas$path, launchNATIVE), gwasSTAGEROOT)
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
}

# Validate immutable reference sources. Tool-ready PLINK and extracted SBayesRC
# resources are generated by the workflow and must not be supplied as run inputs.
referenceREQUIRED <- c(
  "bundle_id", "bundle_version", "reference_id", "reference_type", "path",
  "build", "ancestry", "version", "checksum", "source_format", "reference_stage"
)
requireCOLUMN(reference, referenceREQUIRED, "Generated reference records")
for (column in setdiff(referenceREQUIRED, "checksum")) {
  requireVALUE(reference, column, "Generated reference records")
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
referencePATHNATIVE <- stagedPATH(resolvePATH(reference$path, referenceBASENATIVE), referenceSTAGEROOT)
referenceCOMPANIONNATIVE <- rep("", nrow(reference))
referenceCOMPANIONNATIVE[companionPRESENT] <- resolvePATH(
  reference$companion[companionPRESENT],
  referenceBASENATIVE
)
referenceCOMPANIONNATIVE[companionPRESENT] <- stagedPATH(
  referenceCOMPANIONNATIVE[companionPRESENT],
  referenceSTAGEROOT
)
referenceFILELIST <- vector("list", nrow(reference))
verifiedREFERENCE <- logical(nrow(reference))
observedREFERENCE <- character(nrow(reference))
for (row in seq_len(nrow(reference))) {
  checkPATH(expandPATTERN(referencePATHNATIVE[row]), sprintf("Reference '%s'", reference$reference_id[row]))
  if (companionPRESENT[row]) {
    checkPATH(expandPATTERN(referenceCOMPANIONNATIVE[row]), sprintf("Reference '%s' companion", reference$reference_id[row]))
  }
  referenceFILELIST[[row]] <- referenceFILES(referencePATHNATIVE[row], referenceCOMPANIONNATIVE[row])
  checkPATH(referenceFILELIST[[row]], sprintf("Reference '%s' content", reference$reference_id[row]))
  fileNAME <- basename(referenceFILELIST[[row]])
  contentNAME <- relativeFILENAMES(referenceFILELIST[[row]], referencePATHNATIVE[row])
  referenceTYPE <- reference$reference_type[row]
  dbsnpVCF <- grepl(
    "(\\.vcf(\\.gz|\\.bgz)?|\\.bcf|^GCF_.*\\.gz)$",
    fileNAME,
    ignore.case = TRUE
  )
  assemblyREPORT <- grepl("assembly[_-]?report.*\\.txt$", fileNAME, ignore.case = TRUE)
  if (referenceTYPE == "dbsnp") {
    if (sum(dbsnpVCF) != 1L) {
      stop(
        sprintf(
          "Reference '%s' must contain exactly one dbSNP VCF; found %s.",
          reference$reference_id[row],
          sum(dbsnpVCF)
        ),
        call. = FALSE
      )
    }
    if (sum(assemblyREPORT) != 1L) {
      stop(
        sprintf(
          "Reference '%s' must contain exactly one dbSNP assembly report; found %s.",
          reference$reference_id[row],
          sum(assemblyREPORT)
        ),
        call. = FALSE
      )
    }
  }
  panelCHROMOSOME <- referenceCHROMOSOMES(contentNAME, "imputation_panel")
  if (referenceTYPE == "imputation_panel") {
    checkCHROMOSOMES(
      panelCHROMOSOME,
      reference$reference_id[row],
      "imputation-panel",
      requireCOMPLETE = referenceONLY
    )
  }
  mapNAME <- character()
  mapCHROMOSOME <- integer()
  if (referenceTYPE == "genetic_map") {
    if (reference$source_format[row] == "zip") {
      zipFILE <- referenceFILELIST[[row]][grepl("\\.zip$", referenceFILELIST[[row]], ignore.case = TRUE)]
      if (length(zipFILE) != 1L) {
        stop(
          sprintf("Reference '%s' must resolve to exactly one genetic-map ZIP archive.", reference$reference_id[row]),
          call. = FALSE
        )
      }
      mapNAME <- archiveNAMES(zipFILE[[1L]], sprintf("Reference '%s'", reference$reference_id[row]))
    } else {
      mapNAME <- contentNAME
    }
    mapCHROMOSOME <- referenceCHROMOSOMES(mapNAME, "genetic_map")
    checkCHROMOSOMES(
      mapCHROMOSOME,
      reference$reference_id[row],
      "genetic-map",
      requireCOMPLETE = referenceONLY
    )
  }
  hasCONTENT <- switch(
    referenceTYPE,
    dbsnp = sum(dbsnpVCF) == 1L && sum(assemblyREPORT) == 1L,
    reference_fasta = any(grepl("\\.(fa|fasta|fna)(\\.gz)?$", fileNAME, ignore.case = TRUE)) &&
      any(grepl("\\.fai$", fileNAME, ignore.case = TRUE)),
    genetic_map = length(mapCHROMOSOME) > 0L,
    imputation_panel = length(panelCHROMOSOME) > 0L,
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
  digestBASE <- if (dir.exists(referencePATHNATIVE[row])) referencePATHNATIVE[row] else dirname(referencePATHNATIVE[row])
  observedREFERENCE[row] <- resourceDIGEST(referenceFILELIST[[row]], digestBASE)
  suppliedCHECKSUM <- !is.na(reference$checksum[row]) &&
    nzchar(trimws(reference$checksum[row])) &&
    toupper(trimws(reference$checksum[row])) != "NA"
  if (suppliedCHECKSUM) {
    declaredCHECKSUM <- tolower(trimws(reference$checksum[row]))
    if (!grepl("^[0-9a-f]{64}$", declaredCHECKSUM)) {
      stop(sprintf("Reference '%s' checksum must be a 64-character SHA-256 value.", reference$reference_id[row]), call. = FALSE)
    }
    if (!identical(declaredCHECKSUM, observedREFERENCE[row])) {
      stop(sprintf("Reference '%s' failed its declared SHA-256 check.", reference$reference_id[row]), call. = FALSE)
    }
    reference$checksum[row] <- declaredCHECKSUM
    verifiedREFERENCE[row] <- TRUE
  }
}
writeTABLE(
  data.frame(
    reference_id = reference$reference_id,
    reference_type = reference$reference_type,
    declared_sha256 = ifelse(verifiedREFERENCE, reference$checksum, ""),
    observed_sha256 = observedREFERENCE,
    status = ifelse(verifiedREFERENCE, "AUTHENTICATED", "OBSERVED_NOT_AUTHENTICATED"),
    stringsAsFactors = FALSE
  ),
  "reference_integrity.tsv"
)
if (runPRS && "sbayesrc" %in% method) {
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
missingREFERENCE <- setdiff(unique(requiredREFERENCE), reference$reference_type)
if (length(missingREFERENCE) > 0L) {
  stop(sprintf("Selected methods require reference type(s): %s", paste(missingREFERENCE, collapse = ", ")), call. = FALSE)
}

# Validate phenotype model declarations when phenotype analysis is enabled.
modelREQUIRED <- c(
  "model_id", "outcome", "prs_name", "family", "model_type", "covariates",
  "participant_id", "timepoint_column", "timepoint_values", "group_id",
  "expected_direction", "primary"
)
if (nrow(phenotypeMODEL) > 0L) {
  if (!runPHENOTYPE) stop("Phenotype models were resolved for a run without phenotype analysis enabled.", call. = FALSE)
  for (column in setdiff(c("control_value", "case_value"), names(phenotypeMODEL))) {
    phenotypeMODEL[[column]] <- ""
  }
  requireCOLUMN(phenotypeMODEL, modelREQUIRED, "Generated phenotype-model records")
  for (column in c("model_id", "outcome", "prs_name", "family", "primary")) {
    requireVALUE(phenotypeMODEL, column, "Generated phenotype-model records")
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
    targetID <- targetSAMPLEIDS(target, targetGENOTYPENATIVE, targetSTAGEROOT)
    if (length(targetID) == 0L) {
      stop("Could not read target sample IDs; provide --participant_id explicitly.", call. = FALSE)
    }
    minimumMATCH <- min(length(targetID), max(2L, min(10L, ceiling(length(targetID) * 0.8))))
    candidate <- phenotypeHEADER[vapply(phenotypeDATA, function(value) {
      value <- as.character(value)
      value <- value[!is.na(value) & nzchar(trimws(value))]
      overlap <- length(intersect(unique(value), targetID))
      overlap >= minimumMATCH && overlap / length(targetID) >= 0.8
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
    phenotypeMODEL$timepoint_column[!is.na(phenotypeMODEL$timepoint_column) & phenotypeMODEL$timepoint_column != ""],
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
  allowedFAMILY <- c("gaussian", "binomial")
  phenotypeMODEL$family <- tolower(trimws(phenotypeMODEL$family))
  if (any(!phenotypeMODEL$family %in% allowedFAMILY)) {
    stop("Phenotype model family must be gaussian or binomial after input resolution.", call. = FALSE)
  }
  phenotypeMODEL$model_type <- tolower(trimws(phenotypeMODEL$model_type))
  if (any(!phenotypeMODEL$model_type %in% c("gaussian", "binomial", "mixed"))) {
    stop("Phenotype model_type must be gaussian, binomial, or mixed.", call. = FALSE)
  }
  for (row in seq_len(nrow(phenotypeMODEL))) {
    modelID <- phenotypeMODEL$model_id[[row]]
    outcome <- phenotypeMODEL$outcome[[row]]
    participantCOLUMN <- phenotypeMODEL$participant_id[[row]]
    group <- trimws(as.character(phenotypeMODEL$group_id[[row]]))
    if (is.na(group)) group <- ""
    modelTYPE <- tolower(trimws(as.character(phenotypeMODEL$model_type[[row]])))
    grouped <- identical(modelTYPE, "mixed") || nzchar(group)
    timepointCOLUMN <- trimws(as.character(phenotypeMODEL$timepoint_column[[row]]))
    if (is.na(timepointCOLUMN)) timepointCOLUMN <- ""
    timepointTEXT <- as.character(phenotypeMODEL$timepoint_values[[row]])
    if (is.na(timepointTEXT)) timepointTEXT <- ""
    timepointVALUE <- if (nzchar(trimws(timepointTEXT))) {
      trimws(strsplit(timepointTEXT, "|", fixed = TRUE)[[1L]])
    } else {
      character()
    }
    timepointVALUE <- timepointVALUE[!is.na(timepointVALUE) & nzchar(timepointVALUE)]
    modelDATA <- phenotypeDATA

    if (length(timepointVALUE) > 0L) {
      if (!nzchar(timepointCOLUMN)) {
        stop(sprintf("Model '%s' supplies timepoint_values but no timepoint_column.", modelID), call. = FALSE)
      }
      availableTIMEPOINT <- unique(as.character(phenotypeDATA[[timepointCOLUMN]]))
      availableTIMEPOINT <- availableTIMEPOINT[!is.na(availableTIMEPOINT) & nzchar(availableTIMEPOINT)]
      missingTIMEPOINT <- setdiff(timepointVALUE, availableTIMEPOINT)
      if (length(missingTIMEPOINT) > 0L) {
        stop(
          sprintf(
            "Model '%s' selects timepoint value(s) %s, but column '%s' contains %s. Select available values in timepoint_values.",
            modelID,
            paste(missingTIMEPOINT, collapse = ", "),
            timepointCOLUMN,
            if (length(availableTIMEPOINT) > 0L) paste(availableTIMEPOINT, collapse = ", ") else "no non-missing values"
          ),
          call. = FALSE
        )
      }
      if (!grouped && length(timepointVALUE) != 1L) {
        stop(sprintf("Fixed model '%s' requires exactly one timepoint_values entry.", modelID), call. = FALSE)
      }
      modelDATA <- phenotypeDATA[as.character(phenotypeDATA[[timepointCOLUMN]]) %in% timepointVALUE, , drop = FALSE]
      recordKEY <- paste(modelDATA[[participantCOLUMN]], modelDATA[[timepointCOLUMN]], sep = "\r")
      analysisCOLUMN <- unique(c(
        outcome,
        trimws(strsplit(phenotypeMODEL$covariates[[row]], ",", fixed = TRUE)[[1L]])
      ))
      analysisCOLUMN <- analysisCOLUMN[nzchar(analysisCOLUMN)]
      duplicatedKEY <- unique(recordKEY[duplicated(recordKEY) | duplicated(recordKEY, fromLast = TRUE)])
      for (key in duplicatedKEY) {
        index <- which(recordKEY == key)
        for (column in analysisCOLUMN) {
          observedVALUE <- unique(as.character(modelDATA[[column]][index]))
          observedVALUE <- observedVALUE[!is.na(observedVALUE)]
          if (length(observedVALUE) > 1L) {
            stop(
              sprintf(
                "Model '%s' has conflicting values for participant '%s', field '%s', timepoint '%s': %s. Correct the technical records or timepoint selection.",
                modelID, as.character(modelDATA[[participantCOLUMN]][index[[1L]]]), column,
                as.character(modelDATA[[timepointCOLUMN]][index[[1L]]]), paste(observedVALUE, collapse = ", ")
              ),
              call. = FALSE
            )
          }
        }
      }
      modelDATA <- modelDATA[!duplicated(recordKEY), , drop = FALSE]
    } else if (!grouped && anyDuplicated(as.character(modelDATA[[participantCOLUMN]]))) {
      stop(
        sprintf(
          "Fixed model '%s' has repeated participant IDs. Supply timepoint_column and one timepoint_values entry.",
          modelID
        ),
        call. = FALSE
      )
    }

    outcomeVALUE <- modelDATA[[outcome]]
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
      length(unique(modelDATA[[column]][!is.na(modelDATA[[column]])])) < 2L
    }, logical(1L))]
    if (length(constantCOVARIATE) > 0L) {
      stop(sprintf("Model '%s' has constant or empty covariate(s): %s", modelID, paste(constantCOVARIATE, collapse = ", ")), call. = FALSE)
    }
    if (nzchar(group)) {
      participantGROUP <- as.character(modelDATA[[participantCOLUMN]])
      groupVALUE <- as.character(modelDATA[[group]])
      missingGROUP <- is.na(groupVALUE) | !nzchar(trimws(groupVALUE))
      if (any(missingGROUP)) {
        stop(sprintf("Grouping column '%s' for model '%s' contains missing values. Use the participant grouping field.", group, modelID), call. = FALSE)
      }
      groupingMAP <- unique(data.frame(participant = participantGROUP, group = groupVALUE, stringsAsFactors = FALSE))
      if (anyDuplicated(groupingMAP$participant) || anyDuplicated(groupingMAP$group)) {
        stop(
          sprintf(
            "Grouping column '%s' for model '%s' does not identify the same participants as '%s'. Set group_column to the participant grouping field.",
            group, modelID, participantCOLUMN
          ),
          call. = FALSE
        )
      }
      if (length(unique(groupVALUE)) < 2L) {
        stop(sprintf("Grouping column '%s' for model '%s' has fewer than two groups.", group, modelID), call. = FALSE)
      }
    }
  }
}

# Save validated record tables and a concise input-check record.
writeTABLE(target, "targets.tsv")
writeTABLE(gwas, "gwas.tsv")
writeTABLE(reference, "references.tsv")
writeTABLE(phenotypeMODEL, "models.tsv")
writeTABLE(runPLAN, "run_plan.tsv")

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
  sprintf("reference_base: '%s'", referenceBASE),
  sprintf("stop_after: '%s'", planVALUE[["stop_after"]])
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
    sprintf("%d authenticated; %d observed but not authenticated", sum(verifiedREFERENCE), sum(!verifiedREFERENCE)),
    nrow(phenotypeMODEL), genomeBUILD, paste(method, collapse = ",")
  ),
  status = c("PASS", "PASS", "PASS", "PASS", "PASS", if (all(verifiedREFERENCE)) "PASS" else "REVIEW", "PASS", "PASS", "PASS"),
  stringsAsFactors = FALSE
)
writeTABLE(inputCHECK, "input_checks.tsv")
