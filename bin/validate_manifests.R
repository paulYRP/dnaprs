#!/usr/bin/env Rscript

# Parse named command-line arguments.
argument <- commandArgs(trailingOnly = TRUE)
if (length(argument) %% 2L != 0L || any(!startsWith(argument[seq.int(1L, length(argument), 2L)], "--"))) {
  stop("Arguments must be supplied as --name value pairs.", call. = FALSE)
}
argumentNAME <- sub("^--", "", argument[seq.int(1L, length(argument), 2L)])
argumentVALUE <- argument[seq.int(2L, length(argument), 2L)]
option <- stats::setNames(as.list(argumentVALUE), argumentNAME)

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

expandPATTERN <- function(path) {
  pattern <- gsub("\\{chromosome\\}|\\{chr\\}", "*", path)
  pattern <- gsub("#", "*", pattern, fixed = TRUE)
  if (grepl("[*?]", pattern)) Sys.glob(pattern) else pattern
}

checkTARGET <- function(path, format, label) {
  hasPATTERN <- grepl("\\{chromosome\\}|\\{chr\\}|#|[*?]", path)
  pattern <- gsub("\\{chromosome\\}|\\{chr\\}", "*", path)
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
  } else {
    checkPATH(expanded, label)
  }
}

readHEADER <- function(path) {
  connection <- if (grepl("\\.gz$", path, ignore.case = TRUE)) gzfile(path, "rt") else file(path, "rt")
  on.exit(close(connection))
  firstLINE <- readLines(connection, n = 1L, warn = FALSE)
  if (length(firstLINE) != 1L || firstLINE == "") {
    stop(sprintf("GWAS file has no readable header: %s", path), call. = FALSE)
  }
  strsplit(firstLINE, "\t", fixed = TRUE)[[1L]]
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

# Validate the target contract without repeating completed genotype QC.
targetREQUIRED <- c("cohort", "role", "format", "genotype", "sample", "keep", "build", "ancestry", "dosage")
requireCOLUMN(target, targetREQUIRED, "Target manifest")
for (column in c("cohort", "role", "format", "genotype", "build", "ancestry")) {
  requireVALUE(target, column, "Target manifest")
}
target$format <- tolower(target$format)
target$role <- tolower(target$role)
if (any(!target$format %in% c("pgen", "vcf", "bgen", "bed"))) {
  stop("Target format must be pgen, vcf, bgen, or bed.", call. = FALSE)
}
if (any(!target$role %in% c("target", "validation"))) {
  stop("Target role must be target or validation.", call. = FALSE)
}
if (any(target$build != genomeBUILD)) {
  stop("Every target row must use the selected genome build.", call. = FALSE)
}
target$genotype <- pipelinePATH(target$genotype, launchDIR)
samplePRESENT <- !is.na(target$sample) & trimws(target$sample) != ""
keepPRESENT <- !is.na(target$keep) & trimws(target$keep) != ""
target$sample[!samplePRESENT] <- ""
target$keep[!keepPRESENT] <- ""
target$sample[samplePRESENT] <- pipelinePATH(target$sample[samplePRESENT], launchDIR)
target$keep[keepPRESENT] <- pipelinePATH(target$keep[keepPRESENT], launchDIR)
targetGENOTYPENATIVE <- resolvePATH(target$genotype, launchNATIVE)
targetSAMPLENATIVE <- resolvePATH(target$sample[target$sample != ""], launchNATIVE)
targetKEEPNATIVE <- resolvePATH(target$keep[target$keep != ""], launchNATIVE)
for (row in seq_len(nrow(target))) {
  checkTARGET(targetGENOTYPENATIVE[row], target$format[row], sprintf("Target '%s' genotype", target$cohort[row]))
}
checkPATH(targetSAMPLENATIVE, "Target sample input")
checkPATH(targetKEEPNATIVE, "Target keep input")
if (anyDuplicated(target$cohort)) {
  stop("Target cohort values must be unique.", call. = FALSE)
}

# Validate the GWAS contract and each explicitly mapped source column.
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
  optionalCOLUMN <- unlist(gwas[row, c("freq_col", "n_col")], use.names = FALSE)
  optionalCOLUMN <- optionalCOLUMN[!is.na(optionalCOLUMN) & trimws(optionalCOLUMN) != ""]
  missingCOLUMN <- setdiff(c(declaredCOLUMN, optionalCOLUMN), readHEADER(gwasPATHNATIVE[row]))
  if (length(missingCOLUMN) > 0L) {
    stop(
      sprintf("GWAS '%s' is missing declared column(s): %s", gwas$trait_id[row], paste(missingCOLUMN, collapse = ", ")),
      call. = FALSE
    )
  }
}

# Validate versioned references required by the selected scoring methods.
referenceREQUIRED <- c("reference_id", "reference_type", "path", "build", "ancestry", "version", "checksum")
requireCOLUMN(reference, referenceREQUIRED, "Reference manifest")
for (column in setdiff(referenceREQUIRED, "checksum")) {
  requireVALUE(reference, column, "Reference manifest")
}
if (anyDuplicated(reference$reference_id)) {
  stop("Reference identifiers must be unique.", call. = FALSE)
}
if (any(reference$build != genomeBUILD)) {
  stop("Every reference row must use the selected genome build.", call. = FALSE)
}
reference$reference_type <- tolower(reference$reference_type)
reference$path <- pipelinePATH(reference$path, launchDIR)
referencePATHNATIVE <- resolvePATH(reference$path, launchNATIVE)
for (row in seq_len(nrow(reference))) {
  if (reference$reference_type[row] == "plink_ld") {
    checkTARGET(referencePATHNATIVE[row], "pgen", sprintf("Reference '%s'", reference$reference_id[row]))
  } else {
    checkPATH(expandPATTERN(referencePATHNATIVE[row]), sprintf("Reference '%s'", reference$reference_id[row]))
  }
}
requiredREFERENCE <- character()
if ("plink_ct" %in% method) {
  requiredREFERENCE <- c(requiredREFERENCE, "plink_ld")
}
if ("sbayesrc" %in% method) {
  requiredREFERENCE <- c(requiredREFERENCE, "sbayesrc_ld", "annotation")
  missingFREQUENCY <- is.na(gwas$freq_col) | trimws(gwas$freq_col) == ""
  if (any(missingFREQUENCY)) {
    stop("SBayesRC requires a declared effect-allele-frequency column for every GWAS.", call. = FALSE)
  }
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
  requireCOLUMN(phenotypeMODEL, modelREQUIRED, "Phenotype-model manifest")
  for (column in c("model_id", "outcome", "prs_name", "family", "participant_id", "primary")) {
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
  phenotypeHEADER <- names(utils::read.csv(phenotypePATH, nrows = 0L, check.names = FALSE))
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
}

# Save resolved manifests and a concise preflight record.
writeTABLE(target, "target_manifest.validated.tsv")
writeTABLE(gwas, "gwas_manifest.validated.tsv")
writeTABLE(reference, "reference_manifest.validated.tsv")
writeTABLE(phenotypeMODEL, "phenotype_models.validated.tsv")

yamlVALUE <- c(
  sprintf("run_name: '%s'", option[["run-name"]]),
  sprintf("genome_build: '%s'", genomeBUILD),
  "methods:",
  paste0("  - '", method, "'"),
  sprintf("seed: %s", option[["seed"]]),
  sprintf("report_enabled: %s", tolower(option[["report-enabled"]])),
  sprintf("launch_directory: '%s'", launchDIR)
)
writeLines(yamlVALUE, "resolved_params.yaml")

checksumTABLE <- unique(data.frame(
  display = c(
    basename(option[["target-manifest"]]),
    basename(option[["gwas-manifest"]]),
    basename(option[["reference-manifest"]]),
    target$genotype,
    target$sample[target$sample != ""],
    target$keep[target$keep != ""],
    gwas$path,
    reference$path
  ),
  native = c(
    resolvePATH(option[["target-manifest"]], getwd()),
    resolvePATH(option[["gwas-manifest"]], getwd()),
    resolvePATH(option[["reference-manifest"]], getwd()),
    targetGENOTYPENATIVE,
    targetSAMPLENATIVE,
    targetKEEPNATIVE,
    gwasPATHNATIVE,
    referencePATHNATIVE
  ),
  stringsAsFactors = FALSE
))
checksumTABLE <- checksumTABLE[file.exists(checksumTABLE$native) & !dir.exists(checksumTABLE$native), , drop = FALSE]
checksumDISPLAY <- checksumTABLE$display
checksumPATH <- checksumTABLE$native
checksumCOMMAND <- vapply(
  checksumPATH,
  function(path) {
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
  },
  character(1L)
)
writeTABLE(
  data.frame(path = checksumDISPLAY, algorithm = "SHA-256", checksum = checksumCOMMAND, stringsAsFactors = FALSE),
  "input_checksums.tsv"
)

preflight <- data.frame(
  check = c("target_cohorts", "gwas_traits", "reference_resources", "phenotype_models", "genome_build", "methods"),
  value = c(nrow(target), nrow(gwas), nrow(reference), nrow(phenotypeMODEL), genomeBUILD, paste(method, collapse = ",")),
  status = "PASS",
  stringsAsFactors = FALSE
)
writeTABLE(preflight, "preflight_qc.tsv")
