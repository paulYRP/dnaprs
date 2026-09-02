#!/usr/bin/env Rscript

# Convert conventional input directories or explicit YAML records into the canonical
# TSV records consumed by the scientific modules. These generated files provide
# internal hand-off and provenance; they are not launch parameters.
argument <- commandArgs(trailingOnly = TRUE)
if (length(argument) %% 2L != 0L || any(!startsWith(argument[seq.int(1L, length(argument), 2L)], "--"))) {
  stop("Arguments must be supplied as --name value pairs.", call. = FALSE)
}
name <- sub("^--", "", argument[seq.int(1L, length(argument), 2L)])
value <- argument[seq.int(2L, length(argument), 2L)]
option <- stats::setNames(as.list(value), name)
for (optional in c("beagle-jar", "unbref3-jar")) {
  if (is.null(option[[optional]])) option[[optional]] <- ""
  if (tolower(trimws(option[[optional]])) %in% c("none", "null", "na")) option[[optional]] <- ""
}

for (package in c("jsonlite", "openssl", "yaml")) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop(sprintf("The %s package is required to resolve inputs.", package), call. = FALSE)
  }
}

decodeJSON <- function(value) {
  text <- rawToChar(openssl::base64_decode(value))
  jsonlite::fromJSON(text, simplifyVector = FALSE)
}

inputSPEC <- decodeJSON(option[["input-spec"]])
gwasSPEC <- decodeJSON(option[["gwas-spec"]])
referenceSPEC <- decodeJSON(option[["reference-spec"]])
modelSPEC <- decodeJSON(option[["models-spec"]])
genomeBUILD <- option[["genome-build"]]
method <- strsplit(option[["methods"]], ",", fixed = TRUE)[[1L]]
method <- unique(trimws(method[nzchar(trimws(method))]))
allowedMETHOD <- c("plink_ct", "sbayesrc")
if (length(method) == 0L || any(!method %in% allowedMETHOD)) {
  stop(sprintf("Methods must contain one or more of: %s.", paste(allowedMETHOD, collapse = ", ")), call. = FALSE)
}
stopAFTER <- tolower(trimws(option[["stop-after"]]))
allowedSTAGE <- c("target_qc", "imputation", "prs", "phenotype", "report")
if (!stopAFTER %in% allowedSTAGE) {
  stop(sprintf("stop_after must be one of: %s.", paste(allowedSTAGE, collapse = ", ")), call. = FALSE)
}
referenceONLY <- tolower(option[["reference-only"]]) == "true"
targetIMPUTATION <- tolower(option[["target-imputation"]]) == "true"
runPRS <- !referenceONLY && stopAFTER %in% c("prs", "phenotype", "report")
runPHENOTYPE <- runPRS && stopAFTER %in% c("phenotype", "report") && nzchar(option[["phenotype"]])
runIMPUTATION <- !referenceONLY && targetIMPUTATION && stopAFTER != "target_qc"
if (stopAFTER == "imputation" && !targetIMPUTATION) {
  stop("stop_after=imputation requires target_imputation=true.", call. = FALSE)
}

normaliseSLASH <- function(value) gsub("\\\\", "/", as.character(value))

combinePATH <- function(root, relative) {
  root <- sub("/+$", "", normaliseSLASH(root))
  relative <- sub("^/+", "", normaliseSLASH(relative))
  if (!nzchar(relative) || relative == ".") root else paste(root, relative, sep = "/")
}

relativePATH <- function(path, root) {
  # normalizePath() may leave a relative path unchanged when the final path is a
  # PLINK prefix rather than an existing file. Anchor it first so prefixes such as
  # target_source/cohort are still checked against the staged directory.
  path <- normaliseSLASH(path)
  if (!grepl("^(/|[A-Za-z]:/)", path)) path <- combinePATH(getwd(), path)
  path <- normaliseSLASH(normalizePath(path, winslash = "/", mustWork = FALSE))
  root <- normaliseSLASH(normalizePath(root, winslash = "/", mustWork = TRUE))
  if (identical(path, root)) return("")
  prefix <- paste0(root, "/")
  if (!startsWith(path, prefix)) stop(sprintf("Resolved path is outside its staged root: %s", path), call. = FALSE)
  substring(path, nchar(prefix) + 1L)
}

safeID <- function(value, fallback = "dataset") {
  value <- iconv(as.character(value), to = "ASCII//TRANSLIT", sub = "")
  value <- gsub("[^A-Za-z0-9]+", "_", value)
  value <- gsub("^_+|_+$", "", value)
  if (!nzchar(value)) value <- fallback
  value
}

withoutEXTENSION <- function(value) {
  value <- basename(normaliseSLASH(value))
  sub("(?i)(\\.vcf|\\.bcf|\\.pgen|\\.bed|\\.ped|\\.bgen|\\.tsv|\\.txt|\\.csv)(\\.bgz|\\.gz)?$", "", value, perl = TRUE)
}

writeTABLE <- function(value, path) {
  utils::write.table(value, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
}

recordVALUE <- function(record, candidates, default = "") {
  for (candidate in candidates) {
    value <- record[[candidate]]
    if (!is.null(value) && length(value) > 0L && !is.na(value[[1L]]) && nzchar(as.character(value[[1L]]))) {
      return(value[[1L]])
    }
  }
  default
}

asTEXT <- function(value, collapse = ",") {
  if (is.null(value) || length(value) == 0L) return("")
  paste(unlist(value, use.names = FALSE), collapse = collapse)
}

listDEPTH <- function(root, maximum = 4L) {
  files <- list.files(root, recursive = TRUE, full.names = TRUE, all.files = FALSE, no.. = TRUE)
  relative <- vapply(files, relativePATH, character(1L), root = root)
  depth <- lengths(strsplit(relative, "/", fixed = TRUE))
  files[depth <= maximum]
}

# Discover one coherent raw target dataset. Companion sets are treated as a unit and
# ambiguity is a hard error rather than an alphabetical choice.
discoverTARGET <- function(stagedROOT, originalROOT) {
  if (!dir.exists(stagedROOT)) stop("--input must be a directory for automatic discovery.", call. = FALSE)
  files <- listDEPTH(stagedROOT)
  lower <- tolower(files)
  candidates <- list()

  addSETS <- function(extension, companions, format) {
    primary <- files[endsWith(lower, extension)]
    for (path in primary) {
      prefix <- substr(path, 1L, nchar(path) - nchar(extension))
      expected <- paste0(prefix, companions)
      if (all(file.exists(expected))) {
        candidates[[length(candidates) + 1L]] <<- list(path = prefix, format = format, sample = "", assay = "")
      }
    }
  }
  addSETS(".pgen", c(".pgen", ".pvar", ".psam"), "pgen")
  addSETS(".bed", c(".bed", ".bim", ".fam"), "bed")
  addSETS(".ped", c(".ped", ".map"), "ped")

  for (path in files[grepl("\\.bgen$", lower)]) {
    prefix <- sub("(?i)\\.bgen$", "", path, perl = TRUE)
    sample <- paste0(prefix, ".sample")
    candidates[[length(candidates) + 1L]] <- list(
      path = path, format = "bgen", sample = if (file.exists(sample)) sample else "", assay = ""
    )
  }
  for (path in files[grepl("(\\.vcf(\\.gz|\\.bgz)?|\\.bcf)$", lower)]) {
    candidates[[length(candidates) + 1L]] <- list(path = path, format = "vcf", sample = "", assay = "")
  }

  textFILES <- files[grepl("\\.(txt|csv)(\\.gz)?$", lower)]
  finalREPORT <- textFILES[grepl("final.?report", basename(textFILES), ignore.case = TRUE)]
  assay <- textFILES[grepl("manifest", basename(textFILES), ignore.case = TRUE)]
  if (length(finalREPORT) > 0L) {
    if (length(finalREPORT) != 1L || length(assay) != 1L) {
      stop("GenomeStudio discovery requires exactly one FinalReport and one assay manifest.", call. = FALSE)
    }
    candidates[[length(candidates) + 1L]] <- list(
      path = finalREPORT[[1L]], format = "genomestudio", sample = "", assay = assay[[1L]]
    )
  }

  if (length(candidates) == 0L) {
    stop("The target directory contains no supported coherent PGEN, BED, PED, BGEN, VCF/BCF, or GenomeStudio dataset.", call. = FALSE)
  }
  key <- vapply(candidates, function(item) paste(item$format, item$path, sep = "|"), character(1L))
  candidates <- candidates[!duplicated(key)]
  if (length(candidates) != 1L) {
    display <- vapply(candidates, function(item) sprintf("%s:%s", item$format, relativePATH(item$path, stagedROOT)), character(1L))
    stop(
      sprintf("The target directory resolves to %d datasets (%s). Declare an explicit YAML input list.", length(candidates), paste(display, collapse = ", ")),
      call. = FALSE
    )
  }

  item <- candidates[[1L]]
  relative <- relativePATH(item$path, stagedROOT)
  cohort <- safeID(withoutEXTENSION(item$path), "target")
  data.frame(
    cohort = cohort,
    original_id = withoutEXTENSION(item$path),
    role = "target",
    source_format = item$format,
    genotype = combinePATH(originalROOT, relative),
    sample = if (nzchar(item$sample)) combinePATH(originalROOT, relativePATH(item$sample, stagedROOT)) else "",
    keep = "",
    build = genomeBUILD,
    ancestry = "unspecified",
    dosage = if (item$format == "vcf") "DS" else "",
    input_stage = "raw",
    assay_manifest = if (nzchar(item$assay)) combinePATH(originalROOT, relativePATH(item$assay, stagedROOT)) else "",
    marker_map = "",
    stringsAsFactors = FALSE
  )
}

explicitTARGET <- function(specification) {
  if (!is.list(specification) || length(specification) == 0L) stop("Explicit input must be a non-empty YAML list.", call. = FALSE)
  value <- lapply(specification, function(record) {
    path <- recordVALUE(record, c("path", "genotype"))
    if (!nzchar(path)) stop("Every explicit target record requires path.", call. = FALSE)
    format <- tolower(recordVALUE(record, c("format", "source_format")))
    if (!nzchar(format)) stop("Every explicit target record requires format.", call. = FALSE)
    id <- safeID(recordVALUE(record, c("id", "cohort"), withoutEXTENSION(path)), "target")
    data.frame(
      cohort = id,
      original_id = recordVALUE(record, c("id", "cohort"), withoutEXTENSION(path)),
      role = tolower(recordVALUE(record, "role", "target")),
      source_format = format,
      genotype = path,
      sample = recordVALUE(record, "sample"),
      keep = recordVALUE(record, "keep"),
      build = recordVALUE(record, "build", genomeBUILD),
      ancestry = recordVALUE(record, "ancestry", "unspecified"),
      dosage = recordVALUE(record, "dosage", if (format == "vcf") "DS" else ""),
      input_stage = tolower(recordVALUE(record, c("stage", "input_stage"), "raw")),
      assay_manifest = recordVALUE(record, "assay_manifest"),
      marker_map = recordVALUE(record, "marker_map"),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, value)
}

# Read metadata-prefixed text releases without reading their data payload.
readHEAD <- function(path, maximum = 500L) {
  connection <- if (grepl("\\.(gz|bgz)$", path, ignore.case = TRUE)) gzfile(path, "rt") else file(path, "rt")
  on.exit(close(connection))
  readLines(connection, n = maximum, warn = FALSE)
}

splitHEADER <- function(value) {
  if (grepl("\t", value, fixed = TRUE)) {
    strsplit(value, "\t", fixed = TRUE)[[1L]]
  } else if (grepl(",", value, fixed = TRUE)) {
    strsplit(value, ",", fixed = TRUE)[[1L]]
  } else {
    strsplit(trimws(value), "[[:space:]]+")[[1L]]
  }
}

columnALIAS <- list(
  snp = c("SNP", "ID", "RSID", "MARKERNAME", "VARIANT_ID"),
  chromosome = c("CHROM", "CHR", "CHROMOSOME"),
  position = c("POS", "BP", "POSITION", "BASE_PAIR_LOCATION"),
  effect_allele = c("EA", "A1", "EFFECT_ALLELE", "ALT"),
  other_allele = c("NEA", "A2", "OTHER_ALLELE", "NON_EFFECT_ALLELE", "REF"),
  effect = c("BETA", "LOG_OR", "LOGODDS", "OR", "ODDS_RATIO"),
  standard_error = c("SE", "STDERR", "STANDARD_ERROR"),
  p_value = c("P", "PVAL", "P_VALUE", "PVALUE"),
  frequency = c("EAF", "FREQ", "AF", "ALT_FREQ", "A1FREQ"),
  case_frequency = c("FCAS", "CASE_FREQ", "FRQ_A"),
  control_frequency = c("FCON", "CONTROL_FREQ", "FRQ_U"),
  case_n = c("NCAS", "NCASE", "N_CASE"),
  control_n = c("NCON", "NCONTROL", "N_CONTROL"),
  n = c("N", "NTOT", "NEFF", "N_EFFECTIVE"),
  info = c("IMPINFO", "INFO", "INFO_SCORE", "RSQ")
)

matchCOLUMN <- function(header, role, required = FALSE) {
  normal <- toupper(sub("^#", "", header))
  index <- match(columnALIAS[[role]], normal, nomatch = 0L)
  index <- index[index > 0L]
  if (length(index) == 0L) {
    if (required) stop(sprintf("GWAS header has no recognised %s column.", role), call. = FALSE)
    return("")
  }
  sub("^#", "", header[index[[1L]]])
}

findGWASHEADER <- function(lines, label) {
  parsed <- lapply(lines[nzchar(trimws(lines))], splitHEADER)
  score <- vapply(parsed, function(header) {
    normal <- toupper(sub("^#", "", header))
    required <- c("snp", "chromosome", "position", "effect_allele", "other_allele", "effect", "standard_error", "p_value")
    sum(vapply(required, function(role) any(columnALIAS[[role]] %in% normal), logical(1L)))
  }, integer(1L))
  match <- which(score == 8L)
  if (length(match) != 1L) {
    stop(sprintf("GWAS '%s' has no unique recognised summary-statistic header; declare its columns in params.yml.", label), call. = FALSE)
  }
  sub("^#", "", parsed[[match[[1L]]]])
}

metadataNUMBER <- function(lines, keys) {
  candidate <- lines[grepl(paste(keys, collapse = "|"), lines, ignore.case = TRUE)]
  if (length(candidate) == 0L) return(NA_real_)
  number <- suppressWarnings(as.numeric(gsub(".*[=\"]([0-9.]+).*", "\\1", candidate[[1L]])))
  if (length(number) == 1L && is.finite(number)) number else NA_real_
}

discoverGWASRECORD <- function(path, stagedROOT, originalROOT) {
  lines <- readHEAD(path)
  header <- findGWASHEADER(lines, basename(path))
  metadata <- lines[startsWith(lines, "##")]
  short <- sub('^.*shortName=["]?([^"[:space:]]+).*$','\\1', metadata[grepl("shortName", metadata, ignore.case = TRUE)][1L])
  if (length(short) == 0L || is.na(short) || !nzchar(short) || identical(short, metadata[grepl("shortName", metadata, ignore.case = TRUE)][1L])) {
    short <- withoutEXTENSION(path)
  }
  trait <- safeID(short, "trait")
  betaCOLUMN <- matchCOLUMN(header, "effect", TRUE)
  metadataTEXT <- paste(metadata, collapse = " ")
  effectTYPE <- if (toupper(betaCOLUMN) %in% c("OR", "ODDS_RATIO")) {
    "or"
  } else if (grepl("log odds|ln\\(odds ratio\\)|ln\\(or\\)|log_or", metadataTEXT, ignore.case = TRUE)) {
    "log_or"
  } else {
    "beta"
  }
  nCOLUMN <- matchCOLUMN(header, "n")
  sampleSIZE <- metadataNUMBER(metadata, c("nEffective", "nTotal", "nCase", "nControl"))
  if (!nzchar(nCOLUMN) && !is.finite(sampleSIZE)) {
    stop(sprintf("GWAS '%s' has neither a recognised sample-size column nor declared sample-size metadata.", basename(path)), call. = FALSE)
  }
  if (!is.finite(sampleSIZE)) sampleSIZE <- 1

  data.frame(
    trait_id = trait,
    original_id = short,
    prs_name = paste0(toupper(trait), "_PRS"),
    path = combinePATH(originalROOT, relativePATH(path, stagedROOT)),
    build = genomeBUILD,
    ancestry = "unspecified",
    effect_type = effectTYPE,
    sample_size = format(sampleSIZE, scientific = FALSE, trim = TRUE),
    source_format = "auto",
    snp_col = matchCOLUMN(header, "snp", TRUE),
    chr_col = matchCOLUMN(header, "chromosome", TRUE),
    bp_col = matchCOLUMN(header, "position", TRUE),
    effect_allele_col = matchCOLUMN(header, "effect_allele", TRUE),
    other_allele_col = matchCOLUMN(header, "other_allele", TRUE),
    beta_col = betaCOLUMN,
    se_col = matchCOLUMN(header, "standard_error", TRUE),
    p_col = matchCOLUMN(header, "p_value", TRUE),
    freq_col = matchCOLUMN(header, "frequency"),
    case_freq_col = matchCOLUMN(header, "case_frequency"),
    control_freq_col = matchCOLUMN(header, "control_frequency"),
    case_n_col = matchCOLUMN(header, "case_n"),
    control_n_col = matchCOLUMN(header, "control_n"),
    n_col = nCOLUMN,
    info_col = matchCOLUMN(header, "info"),
    info_min = "",
    maf_min = "0.01",
    stringsAsFactors = FALSE
  )
}

discoverGWAS <- function(stagedROOT, originalROOT) {
  if (!dir.exists(stagedROOT)) stop("--gwas must be a directory for automatic discovery.", call. = FALSE)
  files <- listDEPTH(stagedROOT)
  files <- files[grepl("\\.(tsv|txt|csv|vcf)(\\.gz|\\.bgz)?$", files, ignore.case = TRUE)]
  files <- files[!grepl("(^|[/\\\\])(readme|md5|sha|manifest|receipt)", files, ignore.case = TRUE)]
  if (length(files) == 0L) stop("The GWAS directory contains no supported summary-statistic files.", call. = FALSE)
  records <- lapply(sort(files), discoverGWASRECORD, stagedROOT = stagedROOT, originalROOT = originalROOT)
  value <- do.call(rbind, records)
  if (anyDuplicated(value$trait_id)) {
    stop("Automatically derived GWAS trait identifiers are duplicated; provide explicit YAML IDs.", call. = FALSE)
  }
  value
}

explicitGWAS <- function(specification) {
  if (!is.list(specification) || length(specification) == 0L) stop("Explicit gwas must be a non-empty YAML list.", call. = FALSE)
  value <- lapply(specification, function(record) {
    columns <- record[["columns"]]
    if (is.null(columns)) columns <- list()
    path <- recordVALUE(record, "path")
    id <- safeID(recordVALUE(record, c("id", "trait_id"), withoutEXTENSION(path)), "trait")
    column <- function(role, aliases = character()) recordVALUE(columns, c(role, aliases), recordVALUE(record, c(paste0(role, "_col"), aliases)))
    data.frame(
      trait_id = id,
      original_id = recordVALUE(record, c("id", "trait_id"), id),
      prs_name = recordVALUE(record, "prs_name", paste0(toupper(id), "_PRS")),
      path = path,
      build = recordVALUE(record, "build", genomeBUILD),
      ancestry = recordVALUE(record, "ancestry", "unspecified"),
      effect_type = tolower(recordVALUE(record, "effect_type")),
      sample_size = asTEXT(recordVALUE(record, "sample_size")),
      source_format = recordVALUE(record, "source_format", "auto"),
      snp_col = column("variant", c("snp", "snp_col")),
      chr_col = column("chromosome", c("chr", "chr_col")),
      bp_col = column("position", c("bp", "bp_col")),
      effect_allele_col = column("effect_allele", "effect_allele_col"),
      other_allele_col = column("other_allele", "other_allele_col"),
      beta_col = column("effect", c("beta", "beta_col")),
      se_col = column("standard_error", c("se", "se_col")),
      p_col = column("p_value", c("p", "p_col")),
      freq_col = column("frequency", c("freq", "freq_col")),
      case_freq_col = column("case_frequency", "case_freq_col"),
      control_freq_col = column("control_frequency", "control_freq_col"),
      case_n_col = column("case_n", "case_n_col"),
      control_n_col = column("control_n", "control_n_col"),
      n_col = column("sample_size", c("n", "n_col")),
      info_col = column("info", "info_col"),
      info_min = asTEXT(recordVALUE(record, "info_min")),
      maf_min = asTEXT(recordVALUE(record, "maf_min", "0.01")),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, value)
}

archiveNAMES <- function(path) {
  value <- try(utils::unzip(path, list = TRUE), silent = TRUE)
  if (inherits(value, "try-error")) character() else value$Name
}

oneREFERENCE <- function(values, role) {
  values <- unique(values[file.exists(values) | dir.exists(values)])
  if (length(values) == 0L) return("")
  if (length(values) > 1L) {
    stop(sprintf("Reference discovery found several candidates for role '%s'; declare references.assets.%s explicitly.", role, role), call. = FALSE)
  }
  values[[1L]]
}

discoverREFERENCES <- function(stagedROOT, originalROOT, mode, required) {
  if (!dir.exists(stagedROOT)) {
    if (mode == "local") stop("reference_mode=local requires an existing --references directory.", call. = FALSE)
    return(data.frame())
  }
  files <- list.files(stagedROOT, recursive = TRUE, full.names = TRUE, all.files = FALSE, no.. = TRUE)
  directories <- unique(dirname(files))
  lower <- tolower(files)

  dbDIR <- directories[vapply(directories, function(path) {
    item <- list.files(path, full.names = TRUE)
    any(grepl("\\.(vcf\\.gz|gz)$", item, ignore.case = TRUE)) &&
      any(grepl("\\.tbi$", item, ignore.case = TRUE)) &&
      any(grepl("assembly[_-]?report", item, ignore.case = TRUE))
  }, logical(1L))]
  fasta <- files[grepl("\\.(fa|fasta|fna)$", lower) & file.exists(paste0(files, ".fai"))]
  mapDIR <- directories[vapply(directories, function(path) length(list.files(path, pattern = "\\.map$", ignore.case = TRUE)) >= 22L, logical(1L))]
  panelDIR <- directories[vapply(directories, function(path) length(list.files(path, pattern = "\\.bref3$", ignore.case = TRUE)) >= 22L, logical(1L))]
  population <- files[grepl("\\.(panel|tsv|txt)$", lower) & grepl("population|sample.*panel|integrated_call", lower)]
  related <- files[grepl("related", lower) & grepl("\\.(txt|tsv)$", lower)]
  zipFILES <- files[grepl("\\.zip$", lower)]
  zipCONTENT <- lapply(zipFILES, archiveNAMES)
  ldZIP <- zipFILES[vapply(zipCONTENT, function(value) any(grepl("ldm\\.info$|block[0-9]+\\.eigen\\.bin$", value, ignore.case = TRUE)), logical(1L))]
  annotationZIP <- zipFILES[vapply(zipCONTENT, function(value) any(grepl("annot.*baseline", value, ignore.case = TRUE)), logical(1L))]
  beagleJAR <- files[grepl("beagle[.].*[.]jar$", lower)]
  unbref3JAR <- files[grepl("unbref3[.].*[.]jar$", lower)]

  candidate <- list(
    dbsnp = oneREFERENCE(dbDIR, "dbsnp"),
    reference_fasta = oneREFERENCE(fasta, "reference_fasta"),
    genetic_map = oneREFERENCE(mapDIR, "genetic_map"),
    imputation_panel = oneREFERENCE(panelDIR, "imputation_panel"),
    population_panel = oneREFERENCE(population, "population_panel"),
    related_samples = oneREFERENCE(related, "related_samples"),
    sbayesrc_ld_source = oneREFERENCE(ldZIP, "sbayesrc_ld_source"),
    annotation_source = oneREFERENCE(annotationZIP, "annotation_source"),
    beagle_jar = oneREFERENCE(beagleJAR, "beagle_jar"),
    unbref3_jar = oneREFERENCE(unbref3JAR, "unbref3_jar")
  )
  if (nzchar(option[["beagle-jar"]])) candidate$beagle_jar <- option[["beagle-jar"]]
  if (nzchar(option[["unbref3-jar"]])) candidate$unbref3_jar <- option[["unbref3-jar"]]
  missing <- required[!nzchar(unlist(candidate[required], use.names = FALSE))]
  if (length(missing) > 0L) {
    if (mode == "local") {
      stop(sprintf("The local reference root is missing required role(s): %s", paste(missing, collapse = ", ")), call. = FALSE)
    }
    # Missing automatic assets are downloaded by the REFERENCES subworkflow. Empty rows
    # are not written here, preventing an invalid path from reaching scientific modules.
  }
  present <- names(candidate)[nzchar(unlist(candidate, use.names = FALSE))]
  if (length(present) == 0L) return(data.frame())
  format <- c(
    dbsnp = "directory", reference_fasta = "fasta", genetic_map = "directory",
    imputation_panel = "bref3_directory", population_panel = "tsv",
    related_samples = "tsv", sbayesrc_ld_source = "zip", annotation_source = "zip",
    beagle_jar = "jar", unbref3_jar = "jar"
  )
  data.frame(
    bundle_id = paste0("dnaprs_", genomeBUILD),
    bundle_version = option[["reference-bundle"]],
    reference_id = paste0(toupper(present), "_", genomeBUILD),
    reference_type = present,
    path = vapply(present, function(role) {
      if (role == "beagle_jar" && nzchar(option[["beagle-jar"]])) return(option[["beagle-jar"]])
      if (role == "unbref3_jar" && nzchar(option[["unbref3-jar"]])) return(option[["unbref3-jar"]])
      combinePATH(originalROOT, relativePATH(candidate[[role]], stagedROOT))
    }, character(1L)),
    companion = vapply(present, function(role) if (role == "reference_fasta") paste0(combinePATH(originalROOT, relativePATH(candidate[[role]], stagedROOT)), ".fai") else "", character(1L)),
    build = genomeBUILD,
    ancestry = ifelse(present %in% c("dbsnp", "reference_fasta", "genetic_map"), "All", ifelse(present == "imputation_panel", "Multiple", "European")),
    version = option[["reference-bundle"]],
    checksum = "",
    source_format = unname(format[present]),
    reference_stage = "source",
    stringsAsFactors = FALSE
  )
}

explicitREFERENCES <- function(specification) {
  root <- recordVALUE(specification, "root")
  assets <- specification[["assets"]]
  if (is.null(assets) || length(assets) == 0L) stop("references.assets must declare at least one role.", call. = FALSE)
  aliases <- c(
    dbsnp = "dbsnp", fasta = "reference_fasta", reference_fasta = "reference_fasta",
    genetic_maps = "genetic_map", genetic_map = "genetic_map",
    imputation_panel = "imputation_panel", population_metadata = "population_panel",
    population_panel = "population_panel", related_samples = "related_samples",
    sbayesrc_ld = "sbayesrc_ld_source", sbayesrc_ld_source = "sbayesrc_ld_source",
    sbayesrc_annotation = "annotation_source", annotation_source = "annotation_source",
    beagle_jar = "beagle_jar", unbref3_jar = "unbref3_jar"
  )
  namesASSET <- setdiff(names(assets), "fasta_index")
  role <- unname(aliases[namesASSET])
  if (any(is.na(role))) stop(sprintf("Unknown reference asset role(s): %s", paste(namesASSET[is.na(role)], collapse = ", ")), call. = FALSE)
  if (anyDuplicated(role)) stop("references.assets must declare each reference role once.", call. = FALSE)
  format <- c(
    dbsnp = "directory", reference_fasta = "fasta", genetic_map = "directory",
    imputation_panel = "bref3_directory", population_panel = "tsv",
    related_samples = "tsv", sbayesrc_ld_source = "zip", annotation_source = "zip",
    beagle_jar = "jar", unbref3_jar = "jar"
  )
  path <- vapply(namesASSET, function(name) combinePATH(root, asTEXT(assets[[name]])), character(1L))
  overrides <- c(beagle_jar = "beagle-jar", unbref3_jar = "unbref3-jar")
  for (overrideROLE in names(overrides)) {
    overridePATH <- option[[overrides[[overrideROLE]]]]
    if (!nzchar(overridePATH)) next
    index <- match(overrideROLE, role)
    if (is.na(index)) {
      role <- c(role, overrideROLE)
      path <- c(path, overridePATH)
    } else {
      path[[index]] <- overridePATH
    }
  }
  sourceFORMAT <- unname(format[role])
  sourceFORMAT[role == "genetic_map" & grepl("[.]zip$", path, ignore.case = TRUE)] <- "zip"
  data.frame(
    bundle_id = paste0("dnaprs_", genomeBUILD),
    bundle_version = option[["reference-bundle"]],
    reference_id = paste0(toupper(role), "_", genomeBUILD),
    reference_type = role,
    path = path,
    companion = vapply(seq_along(role), function(index) {
      if (role[[index]] == "reference_fasta" && !is.null(assets[["fasta_index"]])) combinePATH(root, asTEXT(assets[["fasta_index"]])) else ""
    }, character(1L)),
    build = genomeBUILD,
    ancestry = ifelse(role %in% c("dbsnp", "reference_fasta", "genetic_map"), "All", ifelse(role == "imputation_panel", "Multiple", "European")),
    version = option[["reference-bundle"]], checksum = "", source_format = sourceFORMAT,
    reference_stage = "source", stringsAsFactors = FALSE
  )
}

emptyMODELS <- function() data.frame(
  model_id = character(), outcome = character(), prs_name = character(), family = character(),
  covariates = character(), participant_id = character(), group_id = character(),
  expected_direction = character(), primary = logical(), control_value = character(),
  case_value = character(), stringsAsFactors = FALSE
)

buildMODELS <- function(specification, gwas) {
  phenotype <- option[["phenotype"]]
  scalarOUTCOME <- option[["outcome"]]
  hasSCALAR <- nzchar(phenotype) && nzchar(scalarOUTCOME)
  hasLIST <- is.list(specification) && length(specification) > 0L
  if (nzchar(phenotype) && !hasSCALAR && !hasLIST) {
    stop("--phenotype requires --outcome, --covariates, and --model_type, or a YAML models list.", call. = FALSE)
  }
  if (hasSCALAR && (!nzchar(option[["covariates"]]) || !nzchar(option[["model-type"]]))) {
    stop("The scalar phenotype interface requires --phenotype, --outcome, --covariates (or none), and --model_type.", call. = FALSE)
  }
  if (!hasSCALAR && !hasLIST) {
    return(emptyMODELS())
  }
  if (!nzchar(phenotype)) stop("Phenotype models require --phenotype.", call. = FALSE)
  if (hasSCALAR && hasLIST) stop("Use either scalar phenotype-model parameters or models, not both.", call. = FALSE)

  declarations <- if (hasSCALAR) {
    list(list(
      id = safeID(scalarOUTCOME, "model"), outcome = scalarOUTCOME,
      covariates = option[["covariates"]], model_type = option[["model-type"]],
      participant_id = option[["participant-id"]], group_column = option[["group-column"]],
      control_value = option[["control-value"]], case_value = option[["case-value"]], score_ids = "all",
      expected_direction = "", primary = TRUE
    ))
  } else specification

  rows <- list()
  for (record in declarations) {
    outcome <- recordVALUE(record, "outcome")
    family <- tolower(recordVALUE(record, c("model_type", "family"), "gaussian"))
    group <- recordVALUE(record, c("group_column", "group_id"))
    if (family == "mixed") {
      if (!nzchar(group)) stop(sprintf("Mixed model '%s' requires group_column.", recordVALUE(record, "id", outcome)), call. = FALSE)
      family <- "gaussian"
    }
    selected <- record[["score_ids"]]
    if (is.null(selected) || identical(selected, "all") || identical(unlist(selected), "all")) {
      index <- seq_len(nrow(gwas))
    } else {
      selected <- unlist(selected, use.names = FALSE)
      index <- which(gwas$trait_id %in% selected | gwas$prs_name %in% selected)
      missing <- setdiff(selected, c(gwas$trait_id, gwas$prs_name))
      if (length(missing) > 0L) stop(sprintf("Model refers to unknown score ID(s): %s", paste(missing, collapse = ", ")), call. = FALSE)
    }
    baseID <- safeID(recordVALUE(record, c("id", "model_id"), outcome), "model")
    covariateVALUE <- asTEXT(record[["covariates"]])
    if (tolower(trimws(covariateVALUE)) == "none") covariateVALUE <- ""
    for (item in index) {
      rows[[length(rows) + 1L]] <- data.frame(
        model_id = if (length(index) == 1L) baseID else paste(baseID, gwas$trait_id[[item]], sep = "__"),
        outcome = outcome,
        prs_name = gwas$prs_name[[item]],
        family = family,
        covariates = covariateVALUE,
        participant_id = recordVALUE(record, "participant_id", option[["participant-id"]]),
        group_id = group,
        expected_direction = recordVALUE(record, "expected_direction"),
        primary = as.logical(recordVALUE(record, "primary", TRUE)),
        control_value = asTEXT(recordVALUE(record, "control_value", option[["control-value"]])),
        case_value = asTEXT(recordVALUE(record, "case_value", option[["case-value"]])),
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}

emptyTARGET <- function() data.frame(
  cohort = character(), original_id = character(), role = character(), source_format = character(),
  genotype = character(), sample = character(), keep = character(), build = character(),
  ancestry = character(), dosage = character(), input_stage = character(),
  assay_manifest = character(), marker_map = character(), stringsAsFactors = FALSE
)

emptyGWAS <- function() data.frame(
  trait_id = character(), original_id = character(), prs_name = character(), path = character(),
  build = character(), ancestry = character(), effect_type = character(), sample_size = character(),
  snp_col = character(), chr_col = character(), bp_col = character(), effect_allele_col = character(),
  other_allele_col = character(), beta_col = character(), se_col = character(), p_col = character(),
  freq_col = character(), case_freq_col = character(), control_freq_col = character(),
  case_n_col = character(), control_n_col = character(), n_col = character(), info_col = character(),
  info_min = character(), source_format = character(), maf_min = character(), stringsAsFactors = FALSE
)

if (referenceONLY) {
  target <- emptyTARGET()
  gwas <- emptyGWAS()
} else {
  target <- if (is.character(inputSPEC) && length(inputSPEC) == 1L) {
    discoverTARGET(option[["target-staged"]], inputSPEC)
  } else {
    explicitTARGET(inputSPEC)
  }
  gwas <- if (!runPRS) {
    emptyGWAS()
  } else if (is.character(gwasSPEC) && length(gwasSPEC) == 1L) {
      discoverGWAS(option[["gwas-staged"]], gwasSPEC)
    } else {
      explicitGWAS(gwasSPEC)
    }
}

rawTARGET <- referenceONLY || any(target$input_stage == "raw")
requiredREFERENCE <- c("imputation_panel", "population_panel", "related_samples", "unbref3_jar")
if (rawTARGET) requiredREFERENCE <- c(requiredREFERENCE, "dbsnp", "reference_fasta")
if (runIMPUTATION) requiredREFERENCE <- c(requiredREFERENCE, "reference_fasta")
if (referenceONLY || runIMPUTATION) requiredREFERENCE <- c(requiredREFERENCE, "genetic_map", "beagle_jar")
if (referenceONLY || (runPRS && "sbayesrc" %in% method)) {
  requiredREFERENCE <- c(requiredREFERENCE, "sbayesrc_ld_source", "annotation_source")
}
requiredREFERENCE <- unique(requiredREFERENCE)

referenceMODE <- tolower(option[["reference-mode"]])
reference <- if (is.list(referenceSPEC) && !is.character(referenceSPEC) && length(referenceSPEC) > 0L) {
  explicitREFERENCES(referenceSPEC)
} else if (is.character(referenceSPEC) && length(referenceSPEC) == 1L && nzchar(referenceSPEC)) {
  discoverREFERENCES(option[["reference-staged"]], referenceSPEC, referenceMODE, requiredREFERENCE)
} else {
  data.frame()
}

models <- if (runPHENOTYPE) buildMODELS(modelSPEC, gwas) else emptyMODELS()

writeTABLE(target, "targets.tsv")
writeTABLE(gwas, "gwas.tsv")
if (nrow(reference) > 0L) writeTABLE(reference, "references.tsv") else {
  writeTABLE(data.frame(
    bundle_id = character(), bundle_version = character(), reference_id = character(),
    reference_type = character(), path = character(), companion = character(), build = character(),
    ancestry = character(), version = character(), checksum = character(), source_format = character(),
    reference_stage = character(), stringsAsFactors = FALSE
  ), "references.tsv")
}
writeTABLE(models, "models.tsv")
writeTABLE(
  data.frame(
    setting = c(
      "stop_after", "run_imputation", "run_prs", "run_phenotype", "reference_only",
      "required_reference_roles"
    ),
    value = c(
      stopAFTER, tolower(runIMPUTATION), tolower(runPRS), tolower(runPHENOTYPE),
      tolower(referenceONLY), paste(requiredREFERENCE, collapse = ",")
    ),
    stringsAsFactors = FALSE
  ),
  "run_plan.tsv"
)

resolutionPARTS <- list()
if (nrow(target) > 0L) {
  resolutionPARTS[[length(resolutionPARTS) + 1L]] <- data.frame(
    kind = "target", id = target$cohort, source = target$genotype,
    resolution = if (is.character(inputSPEC)) "directory" else "yaml",
    stringsAsFactors = FALSE
  )
}
if (nrow(gwas) > 0L) {
  resolutionPARTS[[length(resolutionPARTS) + 1L]] <- data.frame(
    kind = "gwas", id = gwas$trait_id, source = gwas$path,
    resolution = if (is.character(gwasSPEC)) "directory" else "yaml",
    stringsAsFactors = FALSE
  )
}
if (nrow(reference) > 0L) {
  resolutionPARTS[[length(resolutionPARTS) + 1L]] <- data.frame(
    kind = "reference", id = reference$reference_id, source = reference$path,
    resolution = if (is.character(referenceSPEC)) "directory" else "yaml",
    stringsAsFactors = FALSE
  )
}
resolution <- if (length(resolutionPARTS)) {
  do.call(rbind, resolutionPARTS)
} else {
  data.frame(kind = character(), id = character(), source = character(), resolution = character())
}
writeTABLE(resolution, "input_resolution.tsv")

settings <- list(
  genome = genomeBUILD,
  methods = method,
  reference_mode = referenceMODE,
  reference_bundle = option[["reference-bundle"]],
  reference_only = referenceONLY,
  stop_after = stopAFTER,
  run_imputation = runIMPUTATION,
  run_prs = runPRS,
  run_phenotype = runPHENOTYPE,
  required_reference_roles = requiredREFERENCE,
  target_count = nrow(target),
  gwas_count = nrow(gwas),
  reference_count = nrow(reference),
  model_count = nrow(models)
)
yaml::write_yaml(settings, "reference_settings.yml")
