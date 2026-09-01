# Purpose:
#   Read a completed report input without loading unavailable or empty files.
#
# Args:
#   fileNAME: Basename recorded in the pipeline output manifest.
#
# Returns:
#   A data.table containing the completed result, or an empty data.table.
readRESULT <- function(fileNAME) {
  path <- file.path(inputROOT, fileNAME)
  if (!file.exists(path) || file.info(path)$size == 0) {
    return(data.table())
  }
  fread(path, showProgress = FALSE)
}

# Read and combine completed cohort- or trait-specific result tables selected by name.
readRESULTS <- function(pattern) {
  fileVALUE <- manifest[grepl(pattern, file_name), file_name]
  if (!length(fileVALUE)) return(data.table())
  rbindlist(lapply(fileVALUE, function(fileNAME) {
    value <- readRESULT(fileNAME)
    if (!nrow(value)) return(NULL)
    if (!"cohort" %in% names(value)) value[, cohort := sub("\\..*$", "", fileNAME)]
    value[, source_file := fileNAME]
    value
  }), use.names = TRUE, fill = TRUE)
}

# Purpose:
#   Return the local report path for one completed pipeline file.
#
# Args:
#   fileNAME: Basename recorded in the output manifest.
#
# Returns:
#   A project-relative download path or an empty string.
resultPATH <- function(fileNAME) {
  path <- manifest[file_name == fileNAME, relative_path]
  if (!length(path)) "" else path[[1L]]
}

# Purpose:
#   Encode text for safe insertion into report HTML.
#
# Args:
#   value: Character value to encode.
#
# Returns:
#   HTML-safe character text.
htmlESCAPE <- function(value) {
  value <- gsub("&", "&amp;", as.character(value), fixed = TRUE)
  value <- gsub("<", "&lt;", value, fixed = TRUE)
  value <- gsub(">", "&gt;", value, fixed = TRUE)
  value <- gsub('"', "&quot;", value, fixed = TRUE)
  value
}

# Purpose:
#   Create one accessible file-download control.
#
# Args:
#   path: Project-relative file path.
#   label: Visible control text.
#
# Returns:
#   HTML for one download link or an empty string.
downloadBUTTON <- function(path, label) {
  if (!nzchar(path)) return("")
  sprintf(
    '<a class="dnaprs-download" href="%s" download>%s</a>',
    htmlESCAPE(path),
    htmlESCAPE(label)
  )
}

# Purpose:
#   Create a download control for one pipeline result.
#
# Args:
#   fileNAME: Basename recorded in the output manifest.
#   label: Visible control text.
#
# Returns:
#   HTML for the requested result download.
fileBUTTON <- function(fileNAME, label = paste("Download", fileNAME)) {
  downloadBUTTON(resultPATH(fileNAME), label)
}

# Purpose:
#   Display related download controls beside their report content.
#
# Args:
#   ...: Download-control HTML strings.
#
# Returns:
#   An as-is HTML toolbar.
buttonROW <- function(...) {
  value <- paste(c(...), collapse = "")
  knitr::asis_output(paste0('<div class="dnaprs-toolbar">', value, '</div>'))
}

# Purpose:
#   Add a readable PRS method label while preserving the method identifier.
#
# Args:
#   value: A data.table that may contain a method column.
#
# Returns:
#   The supplied table with method_label when method is available.
addMETHODLABEL <- function(value) {
  if (!"method" %in% names(value)) return(value)
  value[, method_label := unname(methodLABEL[method])]
  value[is.na(method_label), method_label := method]
  value
}

# Purpose:
#   Display a searchable and paginated table without loading large source datasets.
#
# Args:
#   value: Small result or QC table to display.
#   caption: Plain-language table caption.
#   escape: Whether cell HTML must be escaped.
#   page: Initial number of rows shown.
#
# Returns:
#   An interactive HTML table with client-side search and pagination.
tableVIEW <- function(value, caption, escape = TRUE, page = 25L) {
  if (!nrow(value)) {
    return(knitr::asis_output('<p class="dnaprs-note">No records were produced for this section.</p>'))
  }
  tableN <<- tableN + 1L
  tableID <- paste0("tbl-dnaprs-", tableN)
  tableHTML <- knitr::kable(
    value,
    format = "html",
    escape = escape,
    caption = caption,
    table.attr = sprintf('id="%s" class="table table-striped table-sm"', tableID)
  )
  controlHTML <- sprintf(
    paste0(
      '<div class="dnaprs-table-viewer" data-table-viewer>',
      '<div class="dnaprs-table-controls">',
      '<label>Rows <select data-table-size><option>25</option><option>50</option><option>100</option></select></label>',
      '<label class="dnaprs-table-search">Search <input type="search" data-table-search aria-label="Search %s"></label>',
      '<button type="button" data-table-previous aria-label="Previous table page">Previous</button>',
      '<span data-table-status aria-live="polite"></span>',
      '<button type="button" data-table-next aria-label="Next table page">Next</button>',
      '<button type="button" class="dnaprs-expand" data-table-expand aria-label="Expand table" title="Expand table">&#x26F6;</button>',
      '</div><div class="dnaprs-table-wrap">%s</div></div>'
    ),
    htmlESCAPE(caption),
    tableHTML
  )
  knitr::asis_output(controlHTML)
}

# Purpose:
#   Display the files from one or more pipeline result sections with local downloads.
#
# Args:
#   pattern: Regular expression applied to publish_path.
#   caption: Plain-language table caption.
#
# Returns:
#   An interactive file table with contextual downloads.
fileTABLE <- function(pattern, caption) {
  value <- copy(manifest[grepl(pattern, publish_path)])
  if (!nrow(value)) return(tableVIEW(value, caption))
  value[, download := sprintf(
    '<a class="dnaprs-download" href="%s" download>Download</a>',
    relative_path
  )]
  tableVIEW(
    value[, .(section = publish_path, file = file_name, description, download)],
    caption,
    escape = FALSE
  )
}

# Purpose:
#   Display a compact selector for the quantitative figures on one report page.
#
# Args:
#   page: Page value recorded in the figure manifest.
#
# Returns:
#   An interactive figure gallery with interpretation and four download formats.
figureGALLERY <- function(page) {
  pageVALUE <- page
  value <- figureMANIFEST[page == pageVALUE]
  if (!nrow(value)) {
    return(knitr::asis_output('<p class="dnaprs-note">No quantitative figures were produced for this page.</p>'))
  }
  galleryN <<- galleryN + 1L
  galleryID <- paste0("dnaprs-gallery-", galleryN)
  records <- jsonlite::toJSON(value, dataframe = "rows", auto_unbox = TRUE, na = "null")
  options <- paste(sprintf(
    '<option value="%s">%s</option>',
    seq_len(nrow(value)),
    htmlESCAPE(value$title)
  ), collapse = "")
  galleryHTML <- sprintf(
    paste0(
      '<section id="%s" class="dnaprs-gallery" data-figure-gallery>',
      '<div class="dnaprs-gallery-controls">',
      '<button type="button" data-figure-previous aria-label="Previous figure">Previous</button>',
      '<label>Figure <select data-figure-select>%s</select></label>',
      '<button type="button" data-figure-next aria-label="Next figure">Next</button>',
      '<span data-figure-count aria-live="polite"></span>',
      '</div>',
      '<div class="dnaprs-gallery-note"><p data-figure-description></p>',
      '<p><strong>What to inspect.</strong> <span data-figure-inspection></span></p></div>',
      '<div class="dnaprs-toolbar">',
      '<a class="dnaprs-download" data-figure-svg download>Download SVG</a>',
      '<a class="dnaprs-download" data-figure-tiff download>Download TIFF</a>',
      '<a class="dnaprs-download" data-figure-png download>Download PNG</a>',
      '<a class="dnaprs-download" data-figure-jpeg download>Download JPEG</a>',
      '<a class="dnaprs-download" data-figure-source download>Download source table</a>',
      '</div>',
      '<article class="dnaprs-figure-card">',
      '<div class="dnaprs-figure-heading"><h2 data-figure-title></h2>',
      '<div class="dnaprs-figure-actions">',
      '<button type="button" class="dnaprs-zoom-control" data-figure-zoom-out aria-label="Zoom out" title="Zoom out">&#x2212;</button>',
      '<button type="button" class="dnaprs-zoom-control" data-figure-zoom-reset aria-label="Fit figure to panel" title="Fit figure to panel">Fit</button>',
      '<span class="dnaprs-zoom-status" data-figure-zoom-status aria-live="polite">Fit</span>',
      '<button type="button" class="dnaprs-zoom-control" data-figure-zoom-in aria-label="Zoom in" title="Zoom in">+</button>',
      '<button type="button" class="dnaprs-expand" data-figure-expand aria-label="Expand selected figure" title="Expand selected figure">&#x26F6;</button>',
      '</div></div>',
      '<div class="dnaprs-figure-canvas" tabindex="0" role="region" aria-label="Interactive figure">',
      '<div class="dnaprs-figure-stage" data-figure-stage><img data-figure-image alt="" loading="lazy" draggable="false"></div>',
      '</div>',
      '</article>',
      '<script type="application/json" data-figure-data>%s</script>',
      '</section>'
    ),
    galleryID,
    options,
    records
  )
  knitr::asis_output(galleryHTML)
}

inputROOT <- Sys.getenv("DNAPRS_REPORT_INPUTS")
manifest <- fread(Sys.getenv("DNAPRS_OUTPUT_MANIFEST"))
manifest[, relative_path := gsub(
  "\\\\", "/",
  file.path("downloads", publish_path, file_name)
)]
outputFILE <- fread(file.path("provenance", "output_files.tsv"))
manifest <- outputFILE
figureMANIFEST <- fread(file.path("provenance", "figure_manifest.tsv"))

targetMANIFEST <- readRESULT("targets.tsv")
gwasMANIFEST <- readRESULT("gwas.tsv")
referenceMANIFEST <- readRESULT("references.tsv")
phenotypeMANIFEST <- readRESULT("models.tsv")
inputCHECK <- readRESULT("input_checks.tsv")
checksum <- readRESULT("input_checksums.tsv")
score <- readRESULT("prs_scores_long.tsv")
scoreQC <- readRESULT("score_qc.tsv")
concordance <- readRESULT("method_concordance.tsv")
variantFLOW <- readRESULT("variant_flow.tsv")
association <- readRESULT("phenotype_associations.tsv")
fittedMODEL <- readRESULT("phenotype_models_fitted.tsv")
phenotypePERMUTATION <- readRESULT("phenotype_permutations.tsv")
phenotypeINFLUENCE <- readRESULT("phenotype_influence.tsv")
genotypeEDA <- readRESULTS("\\.genotype_eda_summary\\.tsv$")
genotypeEDACHECK <- readRESULTS("\\.genotype_eda_checks\\.tsv$")
targetPREP <- readRESULTS("\\.target_prep_summary\\.tsv$")
targetQC <- readRESULTS("\\.target_qc\\.tsv$")
targetSAMPLEDECISION <- readRESULTS("\\.sample_decisions\\.tsv$")
targetVARIANTDECISION <- readRESULTS("\\.variant_decisions\\.tsv$")
participantDECISION <- readRESULTS("\\.participant_decisions\\.tsv$")
targetANCESTRY <- readRESULTS("\\.target_ancestry\\.tsv$")
referencePROJECTION <- readRESULTS("\\.reference_projection\\.tsv$")
ancestrySUMMARY <- readRESULTS("\\.ancestry_summary\\.tsv$")
imputationQC <- readRESULTS("\\.imputation_qc\\.tsv$")
imputationMANIFEST <- readRESULTS("\\.imputation_manifest\\.tsv$")
imputationDR2 <- readRESULTS("\\.imputation_dr2\\.tsv$")
plinkSENSITIVITY <- readRESULTS("\\.plink_ct\\.sensitivity\\.tsv$")
plinkSENSITIVITYQC <- readRESULTS("\\.plink_ct\\.sensitivity_qc\\.tsv$")
gwasSUMMARY <- readRESULTS("\\.harmonisation_qc\\.tsv$")

methodLABEL <- c(plink_ct = "PLINK C+T", sbayesrc = "SBayesRC")
score <- addMETHODLABEL(score)
scoreQC <- addMETHODLABEL(scoreQC)
variantFLOW <- addMETHODLABEL(variantFLOW)
tableN <- 0L
galleryN <- 0L
