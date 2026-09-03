#!/usr/bin/env Rscript

suppressWarnings(suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(openxlsx)
  library(scales)
}))

options(stringsAsFactors = FALSE, scipen = 999)

DNAPRS_COLOURS <- c(
  teal = "#0A8F78", blue = "#3B6FB6", orange = "#D97706",
  red = "#C74440", purple = "#7556A4", grey = "#68737D"
)

themeDNAPRS <- function(base_size = 12) {
  ggplot2::theme_minimal(base_size = base_size) +
    theme(
      plot.title.position = "plot",
      plot.title = element_text(face = "bold", colour = "#17212B"),
      axis.title = element_text(face = "bold", colour = "#26333F"),
      axis.text = element_text(colour = "#364553"),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = "#E2E8EC", linewidth = .35),
      strip.text = element_text(face = "bold", colour = "#26333F"),
      legend.position = "top"
    )
}

# Existing plot definitions call theme_minimal(); route them through the shared
# report theme so every page uses the same typography, grid, and contrast rules.
theme_minimal <- themeDNAPRS

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

# Read and combine report tables selected by their stable output filename.
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

# Save the exact records used to draw one report figure.
saveFIGUREDATA <- function(id, value) {
  path <- file.path("data", paste0(id, ".tsv"))
  fwrite(value, path, sep = "\t", na = "NA")
  gsub("\\\\", "/", path)
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
#   Save one quantitative plot in vector and publication-oriented raster formats.
#
# Args:
#   id: Stable figure identifier.
#   plot: ggplot object to save.
#   width: Figure width in inches.
#   height: Figure height in inches.
#   page: Report page containing the figure.
#   section: Page section containing the figure.
#   title: Human-readable figure title.
#   description: Technical description of the calculated values.
#   inspection: Plain-language guidance for visual review.
#
# Returns:
#   The four figure files and one row added to figureMANIFEST.
savePLOT <- function(id, plot, width, height, page, section, title,
                     description, inspection, source_table = "") {
  pathSVG <- file.path("figures", "svg", paste0(id, ".svg"))
  pathTIFF <- file.path("figures", "tiff", paste0(id, ".tiff"))
  pathPNG <- file.path("figures", "png", paste0(id, ".png"))
  pathJPEG <- file.path("figures", "jpeg", paste0(id, ".jpeg"))
  grDevices::svg(pathSVG, width = width, height = height, bg = "white", onefile = TRUE)
  tryCatch(print(plot), finally = grDevices::dev.off())
  ggsave(pathTIFF, plot, width = width, height = height, dpi = 360,
         compression = "lzw", bg = "white")
  ggsave(pathPNG, plot, width = width, height = height, dpi = 360,
         bg = "white")
  ggsave(pathJPEG, plot, width = width, height = height, dpi = 360,
         quality = 95, bg = "white")
  figureMANIFEST <<- rbind(
    figureMANIFEST,
    data.table(
      figure_id = id,
      page = page,
      section = section,
      title = title,
      description = description,
      inspection = inspection,
      width_in = width,
      height_in = height,
      dpi = 360L,
      svg = gsub("\\\\", "/", pathSVG),
      tiff = gsub("\\\\", "/", pathTIFF),
      png = gsub("\\\\", "/", pathPNG),
      jpeg = gsub("\\\\", "/", pathJPEG),
      source_table = source_table
    )
  )
}

# Purpose:
#   Describe the analytical role of a published result from its result section.
#
# Args:
#   path: Project-relative result section recorded by Nextflow.
#
# Returns:
#   One concise description for the file dictionary.
describeFILE <- function(path) {
  first <- strsplit(path, "/", fixed = TRUE)[[1L]][1L]
  description <- unname(c(
    run = "Validated run input or resolved setting.",
    genotype_eda = "Exploratory target-genotype result.",
    target_prep = "Target preparation result.",
    target_qc = "Target quality-control result.",
    checkpoints = "Reusable target-stage checkpoint or handoff manifest.",
    qc = "Quality-control result produced by the pipeline.",
    gwas = "Harmonised GWAS result.",
    reference = "Reference-panel result used by scoring.",
    plink_ct = "PLINK clumping-and-thresholding result.",
    sbayesrc = "SBayesRC summary or model result.",
    scores = "Participant polygenic score result.",
    phenotype = "Prespecified phenotype-association result.",
    logs = "Method execution log.",
    pipeline_info = "Pipeline software-version record."
  )[first])
  if (is.na(description)) "Pipeline result." else description
}

# Purpose:
#   Write the Quarto website definition for pages supported by completed results.
#
# Args:
#   pageVALUE: Named vector mapping page identifiers to navigation labels.
#
# Returns:
#   A project configuration that renders only available report pages.
writeCONFIG <- function(pageVALUE) {
  pageFILE <- c(
    overview = "index.qmd",
    genotype_eda = "genotype-eda.qmd",
    target_prep = "target-prep.qmd",
    target_qc = "target-qc.qmd",
    target_imputation = "target-imputation.qmd",
    gwas_qc = "gwas-qc.qmd",
    plink = "plink-prs.qmd",
    sbayesrc = "sbayesrc-prs.qmd",
    phenotype = "phenotype.qmd",
    logs = "logs.qmd"
  )
  renderLINE <- paste0("    - ", unname(pageFILE[names(pageVALUE)]))
  navigation <- unlist(lapply(names(pageVALUE), function(pageID) {
    c(
      paste0("      - href: ", pageFILE[[pageID]]),
      paste0("        text: \"", pageVALUE[[pageID]], "\"")
    )
  }), use.names = FALSE)
  writeLines(
    c(
      "project:",
      "  type: website",
      "  output-dir: _site",
      "  render:",
      renderLINE,
      "  resources:",
      "    - assets/",
      "    - downloads/",
      "    - data/",
      "    - figures/",
      "    - provenance/",
      "website:",
      "  title: \"dnaprs\"",
      "  favicon: assets/nf-core-dnaprs_logo_light.svg",
      "  search: false",
      "  page-navigation: false",
      "  bread-crumbs: false",
      "  reader-mode: false",
      "  navbar:",
      "    logo: assets/nf-core-dnaprs_logo_dark.svg",
      "    right:",
      navigation,
      "format:",
      "  html:",
      "    theme: cosmo",
      "    toc: false",
      "    code-fold: false",
      "    css: assets/dnaprs-report.css",
      "    include-after-body: report-footer.html",
      "execute:",
      "  echo: false",
      "  warning: false",
      "  message: false",
      "  error: false"
    ),
    "_quarto.yml",
    useBytes = TRUE
  )
}

inputROOT <- Sys.getenv("DNAPRS_REPORT_INPUTS")
manifest <- fread(Sys.getenv("DNAPRS_OUTPUT_MANIFEST"))
stopifnot(
  dir.exists(inputROOT),
  all(c("publish_path", "file_name") %in% names(manifest)),
  !anyDuplicated(manifest$file_name),
  !anyDuplicated(file.path(manifest$publish_path, manifest$file_name))
)

dir.create("downloads", recursive = TRUE, showWarnings = FALSE)
dir.create("data", recursive = TRUE, showWarnings = FALSE)
dir.create(file.path("figures", "svg"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path("figures", "tiff"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path("figures", "png"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path("figures", "jpeg"), recursive = TRUE, showWarnings = FALSE)
dir.create("provenance", recursive = TRUE, showWarnings = FALSE)
dir.create("assets", recursive = TRUE, showWarnings = FALSE)
reportASSET <- c(
  "nf-core-dnaprs_logo_light.svg",
  "nf-core-dnaprs_logo_dark.svg",
  "dnaprs-report.css",
  "dnaprs-report.js"
)
stopifnot(all(file.exists(reportASSET)))
stopifnot(all(file.copy(
  reportASSET,
  file.path("assets", reportASSET),
  overwrite = TRUE
)))

manifest[, relative_path := gsub(
  "\\\\", "/",
  file.path("downloads", publish_path, file_name)
)]
for (row in seq_len(nrow(manifest))) {
  sourceP <- file.path(inputROOT, manifest$file_name[row])
  destinationP <- manifest$relative_path[row]
  dir.create(dirname(destinationP), recursive = TRUE, showWarnings = FALSE)
  stopifnot(file.exists(sourceP), file.copy(sourceP, destinationP, overwrite = TRUE))
}

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
phenotypePLOT <- readRESULT("phenotype_plot_data.tsv")
phenotypePERMUTATION <- readRESULT("phenotype_permutations.tsv")
phenotypeINFLUENCE <- readRESULT("phenotype_influence.tsv")
phenotypeDATA <- readRESULT("phenoPRS.csv")
genotypeEDA <- readRESULTS("\\.genotype_eda_summary\\.tsv$")
chromosomeCOUNT <- readRESULTS("\\.chromosome_counts\\.tsv$")
markerDENSITY <- readRESULTS("\\.marker_density\\.tsv$")
identifierCLASS <- readRESULTS("\\.identifier_classes\\.tsv$")
alleleSTATE <- readRESULTS("\\.allele_states\\.tsv$")
sampleMISSING <- readRESULTS("\\.sample_missingness\\.tsv$")
variantMISSING <- readRESULTS("\\.variant_missingness\\.tsv$")
variantMISSINGBIN <- readRESULTS("\\.variant_missingness_bins\\.tsv$")
alleleFREQUENCY <- readRESULTS("\\.allele_frequency\\.tsv$")
alleleFREQUENCYBIN <- readRESULTS("\\.allele_frequency_bins\\.tsv$")
heterozygosity <- readRESULTS("\\.heterozygosity\\.tsv$")
sexCHECK <- readRESULTS("\\.sex_check\\.tsv$")
relatedness <- readRESULTS("\\.relatedness\\.tsv$")
relatednessBIN <- readRESULTS("\\.relatedness_bins\\.tsv$")
internalPCA <- readRESULTS("\\.internal_pca\\.tsv$")
pcaEIGENVALUE <- readRESULTS("\\.pca_eigenvalues\\.tsv$")
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
gwasDATA <- readRESULTS("\\.cojo\\.ma$")
sbayesrcWEIGHT <- readRESULTS("\\.sbayesrc\\.txt$")

methodLABEL <- c(plink_ct = "PLINK C+T", sbayesrc = "SBayesRC")
methodCOLOUR <- c("PLINK C+T" = DNAPRS_COLOURS[["blue"]], "SBayesRC" = DNAPRS_COLOURS[["teal"]])
score <- addMETHODLABEL(score)
scoreQC <- addMETHODLABEL(scoreQC)
variantFLOW <- addMETHODLABEL(variantFLOW)

hasPLINK <- nrow(score) > 0L && any(score$method == "plink_ct")
hasSBAYESRC <- nrow(score) > 0L && any(score$method == "sbayesrc")
hasPHENOTYPE <- nrow(association) > 0L
pageVALUE <- c(
  overview = "Overview"
)
if (nrow(genotypeEDA)) pageVALUE <- c(pageVALUE, genotype_eda = "Genotype EDA")
if (nrow(targetPREP)) pageVALUE <- c(pageVALUE, target_prep = "Target PREP")
if (nrow(targetQC)) pageVALUE <- c(pageVALUE, target_qc = "Target QC")
if (nrow(imputationQC)) pageVALUE <- c(pageVALUE, target_imputation = "Target Imputation")
if (nrow(gwasSUMMARY)) pageVALUE <- c(pageVALUE, gwas_qc = "GWAS QC")
if (hasPLINK) pageVALUE <- c(pageVALUE, plink = "PLINK PRS")
if (hasSBAYESRC) pageVALUE <- c(pageVALUE, sbayesrc = "SBayesRC PRS")
if (hasPHENOTYPE) pageVALUE <- c(pageVALUE, phenotype = "Phenotype")
pageVALUE <- c(pageVALUE, logs = "Logs")
writeCONFIG(pageVALUE)

figureMANIFEST <- data.table(
  figure_id = character(),
  page = character(),
  section = character(),
  title = character(),
  description = character(),
  inspection = character(),
  width_in = numeric(),
  height_in = numeric(),
  dpi = integer(),
  svg = character(),
  tiff = character(),
  png = character(),
  jpeg = character(),
  source_table = character()
)

if (nrow(chromosomeCOUNT)) {
  chromosomeCOUNT[, chromosome := factor(chromosome, levels = unique(chromosome[order(suppressWarnings(as.integer(chromosome)))]))]
  chromosomePLOT <- ggplot(chromosomeCOUNT, aes(x = chromosome, y = variants)) +
    geom_col(fill = "#111111", width = .72) +
    geom_text(aes(label = comma(variants)), vjust = -.45, size = 3) +
    facet_wrap(~cohort, scales = "free_y") +
    scale_y_continuous(labels = comma, expand = expansion(mult = c(0, .16))) +
    labs(x = "Chromosome", y = "Variants") +
    theme_minimal(base_size = 11) +
    theme(panel.grid.major.x = element_blank())
  savePLOT(
    "genotype_variants_by_chromosome", chromosomePLOT, 11, 6,
    "Genotype EDA", "Genome coverage", "Variants by chromosome",
    "Autosomal variants available in each supplied target cohort.",
    "Inspect empty chromosomes and unusually uneven coverage before score generation.",
    saveFIGUREDATA("genotype_variants_by_chromosome", chromosomeCOUNT)
  )
}

if (nrow(markerDENSITY)) {
  markerDENSITY[, `:=`(
    chromosome = factor(as.character(chromosome), levels = as.character(1:22)),
    bin_midpoint_mb = (as.numeric(bin_start) + as.numeric(bin_end)) / 2000000
  )]
  markerDensityPLOT <- ggplot(markerDENSITY, aes(x = bin_midpoint_mb, y = variants)) +
    geom_col(fill = "#111111", width = 1) +
    facet_grid(cohort ~ chromosome, scales = "free_x", space = "free_x") +
    scale_y_continuous(labels = comma, expand = expansion(mult = c(0, .08))) +
    labs(x = "Position within chromosome (Mb)", y = "Variants per 1 Mb bin") +
    theme_minimal(base_size = 9) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      axis.text.x = element_blank(),
      panel.spacing.x = grid::unit(.12, "lines")
    )
  savePLOT(
    "genotype_marker_density", markerDensityPLOT, 14,
    max(5.5, uniqueN(markerDENSITY$cohort) * 2.5),
    "Genotype EDA", "Genome coverage", "Genomic marker density",
    "Raw target markers counted in one-megabase autosomal bins.",
    "Inspect empty chromosomes, long regions without markers, and isolated density spikes.",
    saveFIGUREDATA("genotype_marker_density", markerDENSITY)
  )
}

if (nrow(identifierCLASS)) {
  identifierCLASS[, identifier_class := factor(
    identifier_class,
    levels = c("rsid", "array_probe", "coordinate", "missing", "other")
  )]
  identifierPLOT <- ggplot(identifierCLASS, aes(x = identifier_class, y = variants)) +
    geom_col(fill = "#111111", width = .68) +
    geom_text(aes(label = comma(variants)), vjust = -.4, size = 3) +
    facet_wrap(~cohort, scales = "free_y") +
    scale_y_continuous(labels = comma, expand = expansion(mult = c(0, .16))) +
    labs(x = "Identifier class", y = "Variants") +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1), panel.grid.major.x = element_blank())
  savePLOT(
    "genotype_identifier_classes", identifierPLOT, 10, 6,
    "Genotype EDA", "Variant identifiers and alleles", "Identifier classes",
    "Raw marker identifiers classified as rsID, array-probe-like, coordinate, missing, or other.",
    "Inspect whether assay-manifest-assisted identifier recovery may be required; non-rsID does not itself mean invalid.",
    saveFIGUREDATA("genotype_identifier_classes", identifierCLASS)
  )
}

if (nrow(alleleSTATE)) {
  alleleSTATE[, allele_state := factor(
    allele_state,
    levels = c(
      "acgt_snp", "one_allele_missing", "both_alleles_missing",
      "insertion_deletion", "multiallelic", "other"
    )
  )]
  alleleStatePLOT <- ggplot(alleleSTATE, aes(x = allele_state, y = variants)) +
    geom_col(fill = "#111111", width = .68) +
    geom_text(aes(label = comma(variants)), vjust = -.4, size = 3) +
    facet_wrap(~cohort, scales = "free_y") +
    scale_y_continuous(labels = comma, expand = expansion(mult = c(0, .16))) +
    labs(x = "Allele state", y = "Variants") +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1), panel.grid.major.x = element_blank())
  savePLOT(
    "genotype_allele_states", alleleStatePLOT, 11, 6,
    "Genotype EDA", "Variant identifiers and alleles", "Allele content",
    "Raw target variants classified by allele representation.",
    "Inspect missing alleles, indels, multiallelic records, and other representations before preparation.",
    saveFIGUREDATA("genotype_allele_states", alleleSTATE)
  )
}

if (nrow(sampleMISSING) && "F_MISS" %in% names(sampleMISSING)) {
  sampleMISSING[, participant := if ("IID" %in% names(sampleMISSING)) IID else get("#IID")]
  samplePLOT <- ggplot(sampleMISSING, aes(x = reorder(participant, F_MISS), y = F_MISS)) +
    geom_hline(yintercept = .02, linetype = 2, colour = "#777777") +
    geom_point(colour = "#111111", size = 2.2) +
    facet_wrap(~cohort, scales = "free_x") +
    scale_y_continuous(labels = percent) +
    labs(x = "Participant", y = "Missing genotype proportion") +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 55, hjust = 1))
  savePLOT(
    "genotype_sample_missingness", samplePLOT, 11, 6,
    "Genotype EDA", "Genotype availability", "Sample missingness",
    "Missing genotype proportion for every participant in the supplied target.",
    "Inspect participants near or above the declared missingness threshold.",
    saveFIGUREDATA("genotype_sample_missingness", sampleMISSING)
  )
}

if (nrow(variantMISSING) && "F_MISS" %in% names(variantMISSING)) {
  variantMissingPLOT <- ggplot(variantMISSING, aes(x = F_MISS)) +
    geom_histogram(binwidth = .01, boundary = 0, fill = "#111111", colour = "white") +
    geom_vline(xintercept = .01, linetype = 1, colour = "#777777") +
    geom_vline(xintercept = .10, linetype = 2, colour = "#777777") +
    facet_wrap(~cohort, scales = "free_y") +
    scale_x_continuous(labels = percent) +
    labs(x = "Missing genotype proportion", y = "Variants") +
    theme_minimal(base_size = 11)
  savePLOT(
    "genotype_variant_missingness", variantMissingPLOT, 10, 6,
    "Genotype EDA", "Genotype availability", "Variant missingness",
    "Distribution of missing genotype calls across target variants.",
    "Inspect the mass near zero and the number of variants beyond declared thresholds.",
    saveFIGUREDATA("genotype_variant_missingness", variantMISSING)
  )
}

if (nrow(alleleFREQUENCY) && "ALT_FREQS" %in% names(alleleFREQUENCY)) {
  alleleFREQUENCY[, minor_allele_frequency := pmin(as.numeric(ALT_FREQS), 1 - as.numeric(ALT_FREQS))]
  frequencyPLOT <- ggplot(alleleFREQUENCY[is.finite(minor_allele_frequency)], aes(x = minor_allele_frequency)) +
    geom_histogram(binwidth = .01, boundary = 0, fill = "#111111", colour = "white") +
    facet_wrap(~cohort, scales = "free_y") +
    scale_x_continuous(labels = percent) +
    labs(x = "Minor allele frequency", y = "Variants") +
    theme_minimal(base_size = 11)
  savePLOT(
    "genotype_allele_frequency", frequencyPLOT, 10, 6,
    "Genotype EDA", "Allele frequency", "Minor allele frequency",
    "Target minor allele frequencies calculated by PLINK.",
    "In small cohorts, inspect the discrete frequency pattern without using it alone for exclusion.",
    saveFIGUREDATA("genotype_allele_frequency", alleleFREQUENCY)
  )
}

if (nrow(heterozygosity) && all(c("heterozygosity_rate", "missingness") %in% names(heterozygosity))) {
  heterozygosity[, participant := IID]
  heterozygosityPLOT <- ggplot(
    heterozygosity[is.finite(heterozygosity_rate) & is.finite(missingness)],
    aes(x = missingness, y = heterozygosity_rate, label = participant, colour = status)
  ) +
    geom_point(size = 2.4) +
    geom_text(size = 3, vjust = -.7, check_overlap = TRUE, show.legend = FALSE) +
    facet_wrap(~cohort, scales = "free") +
    scale_colour_manual(values = c(PASS = "#111111", REVIEW = "#777777", FAIL = "#000000")) +
    scale_x_continuous(labels = percent) +
    scale_y_continuous(labels = percent) +
    labs(x = "Missing genotype proportion", y = "Observed heterozygosity", colour = "Status") +
    theme_minimal(base_size = 11) +
    theme(legend.position = "top")
  savePLOT(
    "genotype_heterozygosity_missingness", heterozygosityPLOT, 10, 6,
    "Genotype EDA", "Heterozygosity", "Heterozygosity and missingness",
    "Participant autosomal heterozygosity plotted against raw genotype missingness.",
    "Inspect participants outside the main cluster; ancestry, relatives, contamination, inbreeding, and technical effects can all contribute.",
    saveFIGUREDATA("genotype_heterozygosity_missingness", heterozygosity)
  )
}

if (nrow(sexCHECK) && "x_inbreeding_coefficient" %in% names(sexCHECK)) {
  sexCHECK[, x_inbreeding_coefficient := suppressWarnings(as.numeric(x_inbreeding_coefficient))]
  sexPLOT <- ggplot(
    sexCHECK[is.finite(x_inbreeding_coefficient)],
    aes(x = reorder(IID, x_inbreeding_coefficient), y = x_inbreeding_coefficient, colour = recorded_sex)
  ) +
    geom_hline(yintercept = c(.2, .8), linetype = 2, colour = "#777777") +
    geom_point(size = 2.4) +
    coord_flip() +
    facet_wrap(~cohort, scales = "free_y") +
    labs(x = "Participant", y = "X-chromosome inbreeding coefficient", colour = "Recorded sex") +
    theme_minimal(base_size = 11)
  savePLOT(
    "genotype_reported_sex", sexPLOT, 9, 6,
    "Genotype EDA", "Reported sex", "Recorded and genotype-inferred sex",
    "X-chromosome inbreeding coefficients shown with the recorded sex field.",
    "Inspect ambiguous coefficients and recorded/genetic discordance; the EDA stage does not recode participant sex.",
    saveFIGUREDATA("genotype_reported_sex", sexCHECK)
  )
}

if (nrow(relatedness) && "pi_hat" %in% names(relatedness)) {
  relatedness[, pi_hat_numeric := suppressWarnings(as.numeric(pi_hat))]
  relatednessPLOT <- ggplot(relatedness[is.finite(pi_hat_numeric)], aes(x = pi_hat_numeric)) +
    geom_histogram(bins = 40, fill = "#111111", colour = "white") +
    geom_vline(xintercept = .1875, linetype = 2, colour = "#777777") +
    facet_wrap(~cohort, scales = "free_y") +
    labs(x = "PLINK PI_HAT", y = "Participant pairs") +
    theme_minimal(base_size = 11)
  savePLOT(
    "genotype_relatedness_distribution", relatednessPLOT, 10, 6,
    "Genotype EDA", "Relatedness", "Pairwise relatedness distribution",
    "PLINK 1.9 PI_HAT estimates from the linkage-disequilibrium-pruned marker set.",
    "Pairs with PI_HAT at least 0.1875 are flagged; both participants are excluded from the primary analysis.",
    saveFIGUREDATA("genotype_relatedness_distribution", relatedness)
  )

  relationshipMATRIX <- rbindlist(list(
    relatedness[is.finite(pi_hat_numeric), .(cohort, participant_1 = IID1, participant_2 = IID2, pi_hat = pi_hat_numeric)],
    relatedness[is.finite(pi_hat_numeric), .(cohort, participant_1 = IID2, participant_2 = IID1, pi_hat = pi_hat_numeric)]
  ))
  diagonal <- unique(rbindlist(list(
    relatedness[, .(cohort, participant = IID1)],
    relatedness[, .(cohort, participant = IID2)]
  )))[, .(cohort, participant_1 = participant, participant_2 = participant, pi_hat = 1)]
  relationshipMATRIX <- rbind(relationshipMATRIX, diagonal, fill = TRUE)
  relatednessHeatmapPLOT <- ggplot(
    relationshipMATRIX,
    aes(x = participant_1, y = participant_2, fill = pi_hat)
  ) +
    geom_tile() +
    scale_fill_gradient(low = "white", high = "#111111", limits = c(0, 1), oob = squish) +
    coord_equal() +
    facet_wrap(~cohort, scales = "free") +
    labs(x = "Participant", y = "Participant", fill = "PI_HAT") +
    theme_minimal(base_size = 9) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5),
      panel.grid = element_blank()
    )
  savePLOT(
    "genotype_relatedness_heatmap", relatednessHeatmapPLOT, 9, 8,
    "Genotype EDA", "Relatedness", "Relatedness heatmap",
    "Symmetric participant-by-participant view of PLINK 1.9 PI_HAT.",
    "Inspect flagged pairs while recognising that estimates from small marker sets can be unstable.",
    saveFIGUREDATA("genotype_relatedness_heatmap", relationshipMATRIX)
  )
}

if (nrow(internalPCA) && all(c("PC1", "PC2") %in% names(internalPCA))) {
  internalPCA[, `:=`(PC1 = suppressWarnings(as.numeric(PC1)), PC2 = suppressWarnings(as.numeric(PC2)))]
  pcaPLOT <- ggplot(
    internalPCA[is.finite(PC1) & is.finite(PC2)],
    aes(x = PC1, y = PC2, label = IID)
  ) +
    geom_point(size = 2.4, colour = "#111111") +
    geom_text(size = 3, vjust = -.7, check_overlap = TRUE) +
    facet_wrap(~cohort, scales = "free") +
    labs(x = "Internal PC1", y = "Internal PC2") +
    theme_minimal(base_size = 11)
  savePLOT(
    "genotype_internal_pca", pcaPLOT, 10, 6,
    "Genotype EDA", "Internal population structure", "Internal target PCA",
    "Leading principal components calculated within each target from its exploratory LD-pruned marker set.",
    "Inspect major genotype structure, while recognising that close relatives can dominate small-cohort PCs and that these axes do not assign ancestry.",
    saveFIGUREDATA("genotype_internal_pca", internalPCA)
  )
}

if (nrow(targetPREP)) {
  targetPrepPLOT <- ggplot(targetPREP, aes(x = step, y = variants, fill = cohort)) +
    geom_col(position = position_dodge(width = .72), width = .65) +
    geom_text(aes(label = comma(variants)), position = position_dodge(width = .72), vjust = -.45, size = 3) +
    scale_y_continuous(labels = comma, expand = expansion(mult = c(0, .18))) +
    labs(x = NULL, y = "Variants", fill = "Cohort") +
    theme_minimal(base_size = 11) +
    theme(panel.grid.major.x = element_blank())
  savePLOT(
    "target_preparation_summary", targetPrepPLOT, 9, 5.5,
    "Target PREP", "Preparation outcome", "Target preparation summary",
    "Participants and variants in the normalised target PGEN files.",
    "Confirm that the prepared target contains the expected participants and non-zero autosomal variants.",
    saveFIGUREDATA("target_preparation_summary", targetPREP)
  )
}

if (nrow(targetQC)) {
  targetQCLONG <- melt(
    targetQC,
    id.vars = intersect(c("cohort", "input_stage", "status", "source_file"), names(targetQC)),
    measure.vars = intersect(
      c("source_participants", "retained_participants", "source_variants", "retained_variants", "chromosomes"),
      names(targetQC)
    ),
    variable.name = "measure",
    value.name = "value"
  )
  targetQCPLOT <- ggplot(targetQCLONG, aes(x = cohort, y = value)) +
    geom_col(fill = "#111111", width = .65) +
    geom_text(aes(label = comma(value)), vjust = -.45, size = 3) +
    facet_wrap(~measure, scales = "free_y") +
    scale_y_continuous(labels = comma, expand = expansion(mult = c(0, .18))) +
    labs(x = NULL, y = NULL) +
    theme_minimal(base_size = 11) +
    theme(panel.grid.major.x = element_blank())
  savePLOT(
    "target_qc_summary", targetQCPLOT, 10, 5.5,
    "Target QC", "Scoring readiness", "Target QC summary",
    "Participant, variant, and chromosome counts in the scoring-ready target.",
    "Confirm PASS status and expected final dimensions before reviewing PRS results.",
    saveFIGUREDATA("target_qc_summary", targetQCLONG)
  )
}

if (nrow(targetSAMPLEDECISION)) {
  sampleDECISIONCOUNT <- targetSAMPLEDECISION[, .N, by = .(cohort, decision, reason)]
  sampleDecisionPLOT <- ggplot(sampleDECISIONCOUNT, aes(x = decision, y = N, fill = decision)) +
    geom_col(width = .65) +
    geom_text(aes(label = comma(N)), vjust = -.4, size = 3) +
    facet_wrap(~cohort, scales = "free_y") +
    scale_y_continuous(labels = comma, expand = expansion(mult = c(0, .18))) +
    labs(x = NULL, y = "Participants", fill = "Decision") +
    theme_minimal(base_size = 11) +
    theme(panel.grid.major.x = element_blank(), legend.position = "top")
  savePLOT(
    "target_sample_decisions", sampleDecisionPLOT, 9, 5.5,
    "Target QC", "Participant decisions", "Target participant QC decisions",
    "Participant retention and filtering decisions under the declared target missingness threshold.",
    "Review every filtered participant and all EDA REVIEW records before accepting the checkpoint.",
    saveFIGUREDATA("target_sample_decisions", sampleDECISIONCOUNT)
  )
}

if (nrow(targetVARIANTDECISION)) {
  variantDECISIONCOUNT <- targetVARIANTDECISION[, .N, by = .(cohort, decision, reason)]
  variantDecisionPLOT <- ggplot(variantDECISIONCOUNT, aes(x = reorder(reason, N), y = N, fill = decision)) +
    geom_col(width = .68) +
    coord_flip() +
    facet_wrap(~cohort, scales = "free_y") +
    scale_y_continuous(labels = comma) +
    labs(x = NULL, y = "Variants", fill = "Decision") +
    theme_minimal(base_size = 10) +
    theme(legend.position = "top")
  savePLOT(
    "target_variant_decisions", variantDecisionPLOT, 11, 6.5,
    "Target QC", "Variant decisions", "Target variant QC decisions",
    "Variant retention and filtering reasons from missingness, allele, MAF, HWE, and chromosome checks.",
    "Inspect unexpectedly large losses and confirm that declared thresholds match the analysis plan.",
    saveFIGUREDATA("target_variant_decisions", variantDECISIONCOUNT)
  )
}

if (nrow(participantDECISION) && "primary_analysis" %in% names(participantDECISION)) {
  participantDECISION[, primary_analysis := toupper(as.character(primary_analysis)) == "TRUE"]
  participantDECISIONCOUNT <- participantDECISION[, .N, by = .(
    cohort,
    analysis_set = ifelse(primary_analysis, "Primary analysis", "Sensitivity only"),
    reason
  )]
  participantDecisionPLOT <- ggplot(
    participantDECISIONCOUNT,
    aes(x = analysis_set, y = N, fill = analysis_set)
  ) +
    geom_col(width = .65) +
    geom_text(aes(label = comma(N)), vjust = -.4, size = 3) +
    facet_wrap(~cohort, scales = "free_y") +
    scale_y_continuous(labels = comma, expand = expansion(mult = c(0, .18))) +
    scale_fill_manual(values = c("Primary analysis" = DNAPRS_COLOURS[["teal"]], "Sensitivity only" = DNAPRS_COLOURS[["orange"]])) +
    labs(x = NULL, y = "Participants", fill = "Analysis set") +
    theme_minimal(base_size = 11) +
    theme(panel.grid.major.x = element_blank(), legend.position = "none")
  savePLOT(
    "participant_analysis_eligibility", participantDecisionPLOT, 9, 5.5,
    "Target QC", "Integrated participant decisions", "Primary-analysis eligibility",
    "Participant eligibility after technical QC, relatedness review, and reference-ancestry classification.",
    "Review every sensitivity-only participant and its recorded reason before accepting the primary analysis set.",
    saveFIGUREDATA("participant_analysis_eligibility", participantDECISIONCOUNT)
  )
}

if (nrow(referencePROJECTION) && nrow(targetANCESTRY) && all(c("PC1", "PC2") %in% names(referencePROJECTION)) && all(c("PC1", "PC2") %in% names(targetANCESTRY))) {
  referencePLOT <- copy(referencePROJECTION)
  targetPLOT <- copy(targetANCESTRY)
  referencePLOT[, `:=`(
    PC1 = suppressWarnings(as.numeric(PC1)),
    PC2 = suppressWarnings(as.numeric(PC2)),
    group = if ("super_population" %in% names(referencePLOT)) super_population else "Reference",
    sample_type = "Reference"
  )]
  targetPLOT[, `:=`(
    PC1 = suppressWarnings(as.numeric(PC1)),
    PC2 = suppressWarnings(as.numeric(PC2)),
    group = if ("ancestry_flag" %in% names(targetPLOT)) ancestry_flag else "Target",
    sample_type = "Target"
  )]
  ancestryPLOTDATA <- rbindlist(list(
    referencePLOT[, .(cohort, IID, PC1, PC2, group, sample_type)],
    targetPLOT[, .(cohort, IID, PC1, PC2, group, sample_type)]
  ), use.names = TRUE, fill = TRUE)
  ancestryPLOT <- ggplot(
    ancestryPLOTDATA[is.finite(PC1) & is.finite(PC2)],
    aes(x = PC1, y = PC2, colour = group, shape = sample_type)
  ) +
    geom_point(alpha = .8, size = 2.4) +
    facet_wrap(~cohort, scales = "free") +
    scale_shape_manual(values = c(Reference = 16, Target = 17)) +
    labs(x = "Reference PC1 projection", y = "Reference PC2 projection", colour = "Population / flag", shape = NULL) +
    theme_minimal(base_size = 11)
  savePLOT(
    "reference_ancestry_projection", ancestryPLOT, 10.5, 6.5,
    "Target QC", "Reference ancestry", "Reference-anchored ancestry projection",
    "Target participants projected onto principal components calculated from the unrelated population reference.",
    "Inspect target placement relative to the European reference and any participant classified as an ancestry outlier.",
    saveFIGUREDATA("reference_ancestry_projection", ancestryPLOTDATA)
  )
}

if (nrow(targetANCESTRY) && all(c("ancestry_distance", "ancestry_threshold") %in% names(targetANCESTRY))) {
  ancestryDISTANCE <- copy(targetANCESTRY)
  ancestryDISTANCE[, `:=`(
    ancestry_distance = suppressWarnings(as.numeric(ancestry_distance)),
    ancestry_threshold = suppressWarnings(as.numeric(ancestry_threshold))
  )]
  ancestryDistancePLOT <- ggplot(
    ancestryDISTANCE[is.finite(ancestry_distance)],
    aes(x = reorder(IID, ancestry_distance), y = ancestry_distance, colour = ancestry_flag)
  ) +
    geom_point(size = 2.8) +
    geom_hline(aes(yintercept = ancestry_threshold), linetype = "dashed", colour = DNAPRS_COLOURS[["red"]]) +
    coord_flip() +
    facet_wrap(~cohort, scales = "free_y") +
    labs(x = "Participant", y = "European reference distance", colour = "Classification") +
    theme_minimal(base_size = 10)
  savePLOT(
    "reference_ancestry_distance", ancestryDistancePLOT, 9.5, 6.5,
    "Target QC", "Reference ancestry", "Participant ancestry distance",
    "Mahalanobis distance from the European reference centre; the dashed line is the empirical declared percentile threshold.",
    "Inspect participants beyond the threshold and confirm their sensitivity-only status in the integrated decision table.",
    saveFIGUREDATA("reference_ancestry_distance", ancestryDISTANCE)
  )
}

if (nrow(imputationQC) && "chromosome" %in% names(imputationQC)) {
  imputationLONG <- melt(
    imputationQC,
    id.vars = intersect(c("cohort", "chromosome", "status", "source_file"), names(imputationQC)),
    measure.vars = intersect(c("input_variants", "retained_variants", "imputed_variants"), names(imputationQC)),
    variable.name = "measure",
    value.name = "variants"
  )
  imputationLONG[, variants := suppressWarnings(as.numeric(variants))]
  imputationPLOT <- ggplot(imputationLONG[is.finite(variants)], aes(x = factor(chromosome), y = variants, fill = measure)) +
    geom_col(position = position_dodge(width = .75), width = .68) +
    facet_wrap(~cohort, scales = "free_y") +
    scale_y_continuous(labels = comma) +
    labs(x = "Chromosome", y = "Variants", fill = "Stage") +
    theme_minimal(base_size = 10) +
    theme(legend.position = "top", panel.grid.major.x = element_blank())
  savePLOT(
    "target_imputation_counts", imputationPLOT, 12, 6,
    "Target Imputation", "Chromosome results", "Target imputation variant counts",
    "Typed input, retained post-imputation, and newly imputed variant counts by chromosome.",
    "Inspect chromosome-specific failures and unexpected differences before scoring.",
    saveFIGUREDATA("target_imputation_counts", imputationLONG)
  )
}

if (nrow(imputationDR2) && all(c("chromosome", "dr2_bin", "variants") %in% names(imputationDR2))) {
  imputationDR2[, variants := suppressWarnings(as.numeric(variants))]
  imputationDR2[, dr2_bin := factor(dr2_bin, levels = c("[0,0.3)", "[0.3,0.5)", "[0.5,0.8)", "[0.8,0.9)", "[0.9,1]"))]
  imputationDR2PLOT <- ggplot(
    imputationDR2[is.finite(variants)],
    aes(x = factor(chromosome), y = variants, fill = dr2_bin)
  ) +
    geom_col(width = .72) +
    facet_wrap(~cohort, scales = "free_y") +
    scale_y_continuous(labels = comma) +
    scale_fill_brewer(palette = "YlGnBu", direction = 1, na.value = DNAPRS_COLOURS[["grey"]]) +
    labs(x = "Chromosome", y = "Imputed variants", fill = "DR2 bin") +
    theme_minimal(base_size = 10) +
    theme(panel.grid.major.x = element_blank())
  savePLOT(
    "target_imputation_dr2", imputationDR2PLOT, 12, 6,
    "Target Imputation", "Imputation quality", "Imputation DR2 distribution",
    "Imputed variant counts grouped by Beagle dosage R-squared bin for each chromosome.",
    "Confirm that retained variants are concentrated in the declared acceptable DR2 range and investigate chromosome-specific shifts.",
    saveFIGUREDATA("target_imputation_dr2", imputationDR2)
  )
}

if (nrow(plinkSENSITIVITY) && all(c("imputed_prs_z", "direct_prs_z") %in% names(plinkSENSITIVITY))) {
  plinkSENSITIVITY[, `:=`(
    imputed_prs_z = suppressWarnings(as.numeric(imputed_prs_z)),
    direct_prs_z = suppressWarnings(as.numeric(direct_prs_z))
  )]
  sensitivityPLOT <- ggplot(
    plinkSENSITIVITY[is.finite(imputed_prs_z) & is.finite(direct_prs_z)],
    aes(x = direct_prs_z, y = imputed_prs_z, label = IID)
  ) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = DNAPRS_COLOURS[["grey"]]) +
    geom_point(colour = DNAPRS_COLOURS[["blue"]], size = 2.5) +
    geom_text(size = 2.8, vjust = -.7, check_overlap = TRUE) +
    facet_wrap(~cohort + trait_id) +
    coord_equal() +
    labs(x = "Direct-genotype PRS (z)", y = "Imputed-target PRS (z)") +
    theme_minimal(base_size = 11)
  savePLOT(
    "plink_imputation_sensitivity", sensitivityPLOT, 10, 6.5,
    "PLINK PRS", "Imputation sensitivity", "Imputed versus direct-genotype PRS",
    "Participant PLINK C+T scores calculated from the imputed checkpoint and the stricter direct-genotype checkpoint.",
    "Inspect rank agreement, deviations from the identity line, and participants whose score changes materially after imputation.",
    saveFIGUREDATA("plink_imputation_sensitivity", plinkSENSITIVITY)
  )
}

if (nrow(plinkSENSITIVITYQC)) {
  coverageCOLUMN <- intersect(c("shared_variants", "imputed_only_variants", "direct_only_variants"), names(plinkSENSITIVITYQC))
  if (length(coverageCOLUMN)) {
    sensitivityCOVERAGE <- melt(
      plinkSENSITIVITYQC,
      id.vars = intersect(c("cohort", "trait_id", "prs_name", "source_file"), names(plinkSENSITIVITYQC)),
      measure.vars = coverageCOLUMN,
      variable.name = "coverage",
      value.name = "variants"
    )
    sensitivityCOVERAGE[, variants := suppressWarnings(as.numeric(variants))]
    coveragePLOT <- ggplot(sensitivityCOVERAGE, aes(x = trait_id, y = variants, fill = coverage)) +
      geom_col(width = .68) +
      facet_wrap(~cohort, scales = "free_y") +
      scale_y_continuous(labels = comma) +
      labs(x = "Trait", y = "Scoring variants", fill = "Coverage") +
      theme_minimal(base_size = 10) +
      theme(panel.grid.major.x = element_blank())
    savePLOT(
      "plink_typed_imputed_coverage", coveragePLOT, 10, 6,
      "PLINK PRS", "Imputation sensitivity", "Typed and imputed scoring coverage",
      "Scoring variants shared by the imputed and direct checkpoints or unique to either checkpoint.",
      "Confirm that added imputed coverage is plausible and direct-only variants are absent or understood.",
      saveFIGUREDATA("plink_typed_imputed_coverage", sensitivityCOVERAGE)
    )
  }
}

if (nrow(gwasSUMMARY)) {
  gwasRetention <- melt(
    gwasSUMMARY,
    id.vars = intersect(c("trait_id", "prs_name", "source_file"), names(gwasSUMMARY)),
    measure.vars = intersect(c("source_variants", "harmonised_variants"), names(gwasSUMMARY)),
    variable.name = "stage",
    value.name = "variants"
  )
  gwasRetention[, stage := factor(stage, levels = c("source_variants", "harmonised_variants"), labels = c("Source", "Post-QC"))]
  gwasRetentionPLOT <- ggplot(gwasRetention, aes(x = trait_id, y = variants, fill = stage)) +
    geom_col(position = position_dodge(width = .72), width = .65) +
    geom_text(aes(label = comma(variants)), position = position_dodge(width = .72), vjust = -.45, size = 3) +
    scale_y_continuous(labels = comma, expand = expansion(mult = c(0, .18))) +
    labs(x = "GWAS", y = "Variants", fill = "Stage") +
    theme_minimal(base_size = 11) +
    theme(panel.grid.major.x = element_blank(), legend.position = "top")
  savePLOT(
    "gwas_variant_retention", gwasRetentionPLOT, 10, 6,
    "GWAS QC", "Filtering results", "GWAS variant retention",
    "Source and analysis-ready variants for each GWAS.",
    "Inspect the number and percentage retained and review any unexpected large loss.",
    saveFIGUREDATA("gwas_variant_retention", gwasRetention)
  )
}

if (nrow(gwasDATA)) {
  gwasDATA[, trait_id := sub("\\.cojo\\.ma$", "", source_file)]
  gwasQQ <- gwasDATA[is.finite(p) & p > 0 & p <= 1, {
    observed <- -log10(sort(p))
    .(expected = -log10(stats::ppoints(.N)), observed)
  }, by = trait_id]
  qqPLOT <- ggplot(gwasQQ, aes(x = expected, y = observed)) +
    geom_abline(slope = 1, intercept = 0, colour = "#777777", linetype = 2) +
    geom_point(colour = "#111111", alpha = .55, size = 1.2) +
    facet_wrap(~trait_id, scales = "free") +
    labs(x = "Expected -log10(P)", y = "Observed -log10(P)") +
    theme_minimal(base_size = 11)
  savePLOT(
    "gwas_qq", qqPLOT, 10, 6,
    "GWAS QC", "Source GWAS diagnostics", "GWAS QQ plots",
    "Expected and observed GWAS association probabilities.",
    "Inspect broad departures from the diagonal while retaining the discovery-study context.",
    saveFIGUREDATA("gwas_qq", gwasQQ)
  )

  gwasEFFECTPLOT <- ggplot(gwasDATA[is.finite(b)], aes(x = b)) +
    geom_histogram(bins = 50, fill = "#111111", colour = "white") +
    geom_vline(xintercept = 0, colour = "#777777", linetype = 2) +
    facet_wrap(~trait_id, scales = "free_y") +
    labs(x = "GWAS effect estimate", y = "Variants") +
    theme_minimal(base_size = 11)
  savePLOT(
    "gwas_effect_distribution", gwasEFFECTPLOT, 10, 6,
    "GWAS QC", "Source GWAS diagnostics", "GWAS effect distributions",
    "Signed analysis-ready GWAS effect estimates.",
    "Inspect centring around zero, extreme effects, and marked trait differences.",
    saveFIGUREDATA("gwas_effect_distribution", gwasDATA[, .(trait_id, SNP, CHR, BP, b, se, p, freq, N)])
  )

  gwasMAF <- gwasDATA[is.finite(freq) & freq > 0 & freq < 1]
  gwasMAF[, minor_allele_frequency := pmin(freq, 1 - freq)]
  if (nrow(gwasMAF)) {
    gwasMAFPLOT <- ggplot(gwasMAF, aes(x = minor_allele_frequency)) +
      geom_histogram(bins = 50, fill = DNAPRS_COLOURS[["teal"]], colour = "white") +
      facet_wrap(~trait_id, scales = "free_y") +
      scale_x_continuous(labels = label_number(accuracy = .01)) +
      coord_cartesian(xlim = c(0, .5)) +
      labs(x = "Minor allele frequency", y = "Variants") +
      theme_minimal(base_size = 11)
    savePLOT(
      "gwas_minor_allele_frequency", gwasMAFPLOT, 10, 6,
      "GWAS QC", "Source GWAS diagnostics", "GWAS minor allele frequencies",
      "Minor allele frequencies derived from the analysis-ready effect-allele frequencies.",
      "Inspect the spectrum, rare-variant concentration, and marked differences among GWAS inputs.",
      saveFIGUREDATA("gwas_minor_allele_frequency", gwasMAF[, .(trait_id, SNP, CHR, BP, freq, minor_allele_frequency)])
    )
  }

  gwasMANHATTAN <- gwasDATA[
    is.finite(suppressWarnings(as.numeric(CHR))) & is.finite(BP) & is.finite(p) & p > 0 & p <= 1
  ]
  if (nrow(gwasMANHATTAN)) {
    gwasMANHATTAN[, chromosome := as.integer(CHR)]
    chromosomeLAYOUT <- gwasMANHATTAN[, .(chromosome_length = max(BP)), by = chromosome][order(chromosome)]
    chromosomeLAYOUT[, chromosome_offset := shift(cumsum(chromosome_length), fill = 0)]
    chromosomeLAYOUT[, chromosome_midpoint := chromosome_offset + chromosome_length / 2]
    gwasMANHATTAN <- merge(gwasMANHATTAN, chromosomeLAYOUT, by = "chromosome", all.x = TRUE)
    gwasMANHATTAN[, `:=`(
      genomic_position = chromosome_offset + BP,
      log10_p = -log10(p),
      chromosome_group = factor(chromosome %% 2L)
    )]
    manhattanPLOT <- ggplot(
      gwasMANHATTAN,
      aes(x = genomic_position, y = log10_p, colour = chromosome_group)
    ) +
      geom_hline(yintercept = -log10(5e-8), colour = DNAPRS_COLOURS[["red"]], linetype = 2, linewidth = .55) +
      geom_point(alpha = .72, size = 1.15) +
      facet_wrap(~trait_id, scales = "free_y") +
      scale_colour_manual(values = c(DNAPRS_COLOURS[["blue"]], DNAPRS_COLOURS[["teal"]]), guide = "none") +
      scale_x_continuous(
        breaks = chromosomeLAYOUT$chromosome_midpoint,
        labels = chromosomeLAYOUT$chromosome,
        expand = expansion(mult = c(.01, .01))
      ) +
      labs(x = "Chromosome", y = "-log10(P)") +
      theme_minimal(base_size = 11) +
      theme(panel.grid.major.x = element_blank())
    savePLOT(
      "gwas_manhattan", manhattanPLOT, 14, max(6, uniqueN(gwasMANHATTAN$trait_id) * 2.8),
      "GWAS QC", "Source GWAS diagnostics", "GWAS Manhattan plots",
      "Genome-wide association probabilities across autosomes; the dashed line marks P = 5e-8.",
      "Inspect the genomic distribution of association signals and isolated extreme values in each source GWAS.",
      saveFIGUREDATA("gwas_manhattan", gwasMANHATTAN[, .(trait_id, SNP, chromosome, BP, genomic_position, p, log10_p)])
    )
  }
}

stageLEVEL <- c(
  "Source GWAS", "Harmonised", "LD-aligned", "LD-clumped", "Summary imputed",
  "Model weights", "Non-zero effects", "Scored"
)
variantFLOW[, stage := factor(stage, levels = stageLEVEL)]
variantFlowPLOT <- ggplot(
  variantFLOW,
  aes(
    x = stage,
    y = variant_count,
    group = interaction(cohort, trait_id, method),
    colour = method_label
  )
) +
  geom_line(linewidth = .8) +
  geom_point(size = 2.5) +
  geom_text(
    aes(label = label_number(big.mark = ",")(variant_count)),
    vjust = -.75,
    size = 3,
    show.legend = FALSE
  ) +
  facet_grid(cohort + prs_name ~ method_label, scales = "free_x", space = "free_x") +
  scale_colour_manual(values = methodCOLOUR, guide = "none") +
  scale_y_log10(
    labels = label_number(big.mark = ","),
    expand = expansion(mult = c(.06, .18))
  ) +
  labs(x = NULL, y = "Variants (log scale)") +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 28, hjust = 1),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank()
  )
savePLOT(
  "variant_generation_flow",
  variantFlowPLOT,
  12,
  max(6, uniqueN(variantFLOW[, .(cohort, prs_name)]) * 2.2),
  "Overview",
  "Variant flow",
  "Variant generation flow",
  "Variant counts across GWAS harmonisation, method-specific modelling, and final scoring.",
  "Inspect large unexpected losses, differences between traits, and the final number of variants that contributed to scoring.",
  saveFIGUREDATA("variant_generation_flow", variantFLOW)
)

for (methodVALUE in intersect(c("plink_ct", "sbayesrc"), unique(score$method))) {
  methodPAGE <- if (methodVALUE == "plink_ct") "PLINK PRS" else "SBayesRC PRS"
  methodID <- if (methodVALUE == "plink_ct") "plink" else "sbayesrc"
  methodSCORE <- score[method == methodVALUE]
  methodQC <- scoreQC[method == methodVALUE]

  coveragePLOT <- ggplot(methodQC, aes(x = prs_name, y = used_variants)) +
    geom_col(width = .68, fill = "#111111", na.rm = TRUE) +
    geom_text(aes(label = comma(used_variants)), hjust = -.08, size = 3.2, na.rm = TRUE) +
    coord_flip(clip = "off") +
    facet_wrap(~cohort, scales = "free_x") +
    scale_y_continuous(labels = comma, expand = expansion(mult = c(0, .24))) +
    labs(x = NULL, y = "Variants contributing to the score") +
    theme_minimal(base_size = 11) +
    theme(panel.grid.major.y = element_blank())
  savePLOT(
    paste0(methodID, "_scoring_variant_counts"), coveragePLOT, 10,
    max(5, uniqueN(methodQC$cohort) * 2.6), methodPAGE, "Target scoring coverage",
    paste(methodLABEL[[methodVALUE]], "scoring variant counts"),
    "Variants contributing to each participant score.",
    "Compare coverage across traits and cohorts and review low counts with the variant-flow table.",
    saveFIGUREDATA(paste0(methodID, "_scoring_variant_counts"), methodQC)
  )

  distributionPLOT <- ggplot(methodSCORE, aes(x = prs_z, y = prs_name)) +
    geom_vline(xintercept = 0, colour = "grey75", linewidth = .4) +
    geom_boxplot(width = .32, outlier.shape = NA, alpha = .12, colour = "grey45") +
    geom_jitter(height = .12, width = 0, alpha = .75, size = 2, colour = "#111111") +
    facet_wrap(~cohort, scales = "free_y") +
    labs(x = "PRS (within-cohort standard deviations)", y = NULL) +
    theme_minimal(base_size = 11) +
    theme(panel.grid.major.y = element_blank())
  savePLOT(
    paste0(methodID, "_participant_prs_distributions"), distributionPLOT, 11,
    max(6, uniqueN(methodSCORE$cohort) * 3.2), methodPAGE, "Participant PRS review",
    paste(methodLABEL[[methodVALUE]], "participant PRS distributions"),
    "Participant-level standardised PRS values.",
    "Inspect every participant for missing values, isolated observations, or compressed score ranges.",
    saveFIGUREDATA(paste0(methodID, "_participant_prs_distributions"), methodSCORE)
  )

  methodSCORE[, score_label := prs_name]
  correlationLIST <- lapply(unique(methodSCORE$cohort), function(cohortVALUE) {
    wide <- dcast(methodSCORE[cohort == cohortVALUE], IID ~ score_label, value.var = "prs_z")
    scoreCOLUMN <- setdiff(names(wide), "IID")
    matrixVALUE <- stats::cor(wide[, ..scoreCOLUMN], use = "pairwise.complete.obs")
    result <- as.data.table(as.table(matrixVALUE))
    setnames(result, c("score_1", "score_2", "correlation"))
    result[, cohort := cohortVALUE]
    result
  })
  correlation <- rbindlist(correlationLIST)
  correlationPLOT <- ggplot(correlation, aes(x = score_1, y = score_2, fill = correlation)) +
    geom_tile(colour = "white", linewidth = .5) +
    geom_text(aes(label = number(correlation, accuracy = .01)), size = 3) +
    facet_wrap(~cohort) +
    scale_fill_gradient2(low = "#202020", mid = "white", high = "#8a8a8a", midpoint = 0, limits = c(-1, 1)) +
    coord_equal() +
    labs(x = NULL, y = NULL, fill = "Pearson r") +
    theme_minimal(base_size = 10) +
    theme(axis.text.x = element_text(angle = 38, hjust = 1), panel.grid = element_blank())
  savePLOT(
    paste0(methodID, "_cross_score_correlations"), correlationPLOT, 10, 8,
    methodPAGE, "Participant PRS review", paste(methodLABEL[[methodVALUE]], "cross-score correlations"),
    "Pearson correlations among standardised trait scores from the same method.",
    "Inspect shared information among scores; correlation describes similarity rather than validation.",
    saveFIGUREDATA(paste0(methodID, "_cross_score_correlations"), correlation)
  )
}

if (nrow(sbayesrcWEIGHT) && all(c("BETA", "PIP") %in% names(sbayesrcWEIGHT))) {
  sbayesrcWEIGHT[, trait_id := sub("\\.sbayesrc\\.txt$", "", source_file)]
  posteriorEFFECT <- sbayesrcWEIGHT[is.finite(BETA)]
  if (nrow(posteriorEFFECT)) {
    posteriorEffectPLOT <- ggplot(posteriorEFFECT, aes(x = BETA)) +
      geom_histogram(bins = 60, fill = DNAPRS_COLOURS[["teal"]], colour = "white") +
      geom_vline(xintercept = 0, colour = DNAPRS_COLOURS[["grey"]], linetype = 2) +
      facet_wrap(~trait_id, scales = "free_y") +
      labs(x = "Posterior joint effect", y = "Variants") +
      theme_minimal(base_size = 11)
    savePLOT(
      "sbayesrc_posterior_effect_distribution", posteriorEffectPLOT, 10, 6,
      "SBayesRC PRS", "Posterior model", "SBayesRC posterior effect distributions",
      "Signed posterior joint effects estimated by the completed SBayesRC model.",
      "Inspect centring, tails, zero effects, and large isolated weights for every trait.",
      saveFIGUREDATA("sbayesrc_posterior_effect_distribution", posteriorEFFECT)
    )
  }
  posteriorPIP <- sbayesrcWEIGHT[is.finite(PIP) & PIP >= 0 & PIP <= 1]
  if (nrow(posteriorPIP)) {
    posteriorPIPPLOT <- ggplot(posteriorPIP, aes(x = PIP)) +
      geom_histogram(bins = 50, fill = DNAPRS_COLOURS[["purple"]], colour = "white") +
      facet_wrap(~trait_id, scales = "free_y") +
      scale_x_continuous(labels = percent) +
      coord_cartesian(xlim = c(0, 1)) +
      labs(x = "Posterior inclusion probability", y = "Variants") +
      theme_minimal(base_size = 11)
    savePLOT(
      "sbayesrc_pip_distribution", posteriorPIPPLOT, 10, 6,
      "SBayesRC PRS", "Posterior model", "SBayesRC posterior inclusion probabilities",
      "Posterior inclusion probabilities for variants in the completed SBayesRC model.",
      "Inspect concentration near zero and the number and range of variants with appreciable posterior support.",
      saveFIGUREDATA("sbayesrc_pip_distribution", posteriorPIP)
    )
  }
}

if (hasPLINK && hasSBAYESRC && nrow(concordance) > 0L) {
  agreementLIST <- lapply(seq_len(nrow(concordance)), function(row) {
    record <- concordance[row]
    left <- score[
      cohort == record$cohort & trait_id == record$trait_id & method == record$method_1,
      .(IID, method_1_z = prs_z)
    ]
    right <- score[
      cohort == record$cohort & trait_id == record$trait_id & method == record$method_2,
      .(IID, method_2_z = prs_z)
    ]
    result <- merge(left, right, by = "IID")
    result[, `:=`(cohort = record$cohort, trait_id = record$trait_id)]
    result
  })
  agreement <- rbindlist(agreementLIST)
  agreementPLOT <- ggplot(agreement, aes(x = method_1_z, y = method_2_z)) +
    geom_hline(yintercept = 0, colour = "grey80") +
    geom_vline(xintercept = 0, colour = "grey80") +
    geom_point(colour = "#202020", alpha = .8, size = 2) +
    geom_smooth(
      method = "lm",
      formula = y ~ x,
      se = TRUE,
      colour = "#666666",
      linewidth = .7
    ) +
    facet_grid(cohort ~ trait_id, scales = "free") +
    labs(x = "PLINK C+T PRS", y = "SBayesRC PRS") +
    theme_minimal(base_size = 11)
  savePLOT(
    "method_agreement",
    agreementPLOT,
    10,
    max(5.5, uniqueN(agreement$cohort) * 3),
    "SBayesRC PRS",
    "Method agreement",
    "Method agreement",
    "PLINK C+T and SBayesRC participant scores.",
    "Inspect whether the two methods rank participants similarly and whether any participant departs strongly from the general relationship.",
    saveFIGUREDATA("method_agreement", agreement)
  )
}

if (hasPLINK && hasSBAYESRC) {
  score[, percentile_rank := (frank(prs_z, ties.method = "average") - .5) / .N,
        by = .(cohort, trait_id, method)]
  rankPLOT <- ggplot(
    score,
    aes(
      x = method_label,
      y = percentile_rank,
      group = interaction(cohort, trait_id, IID)
    )
  ) +
    geom_line(colour = "#777777", alpha = .25, linewidth = .45) +
    geom_point(aes(colour = method_label), alpha = .72, size = 1.8) +
    facet_grid(cohort ~ prs_name) +
    scale_colour_manual(values = methodCOLOUR, guide = "none") +
    scale_y_continuous(labels = percent, limits = c(0, 1)) +
    labs(x = NULL, y = "Within-cohort percentile rank") +
    theme_minimal(base_size = 11) +
    theme(panel.grid.major.x = element_blank())
  savePLOT(
    "participant_method_ranks",
    rankPLOT,
    10,
    max(5.5, uniqueN(score$cohort) * 3),
    "SBayesRC PRS",
    "Method agreement",
    "Participant method ranks",
    "Participant rank under each PRS method.",
    "Inspect large rank changes between methods. Rank agreement supports consistency but does not establish prediction.",
    saveFIGUREDATA("participant_method_ranks", score[, .(cohort, trait_id, prs_name, IID, method, method_label, prs_z, percentile_rank)])
  )
}

if (hasPHENOTYPE) {
  outcomeCOLUMN <- intersect(unique(phenotypeMANIFEST$outcome), names(phenotypeDATA))
  participantCOLUMN <- intersect(unique(phenotypeMANIFEST$participant_id), names(phenotypeDATA))
  if (length(outcomeCOLUMN) && length(participantCOLUMN) == 1L) {
    phenotypeOUTCOME <- copy(phenotypeDATA[, c(participantCOLUMN, outcomeCOLUMN), with = FALSE])
    phenotypeOUTCOME[, (outcomeCOLUMN) := lapply(
      .SD,
      function(value) suppressWarnings(as.numeric(value))
    ), .SDcols = outcomeCOLUMN]
    phenotypeOUTCOME <- melt(
      phenotypeOUTCOME,
      id.vars = participantCOLUMN,
      measure.vars = outcomeCOLUMN,
      variable.name = "outcome",
      value.name = "recorded_value",
      variable.factor = FALSE
    )
    setnames(phenotypeOUTCOME, participantCOLUMN, "IID")
    phenotypeOUTCOME <- phenotypeOUTCOME[is.finite(recorded_value)]
    if (nrow(phenotypeOUTCOME)) {
      phenotypeDistributionPLOT <- ggplot(
        phenotypeOUTCOME,
        aes(x = recorded_value, y = outcome)
      ) +
        geom_boxplot(width = .34, outlier.shape = NA, fill = "#E7E0F0", colour = DNAPRS_COLOURS[["purple"]]) +
        geom_jitter(height = .11, width = 0, alpha = .76, size = 1.9, colour = DNAPRS_COLOURS[["purple"]]) +
        labs(x = "Recorded phenotype value", y = NULL) +
        theme_minimal(base_size = 11) +
        theme(panel.grid.major.y = element_blank())
      savePLOT(
        "phenotype_outcome_distributions", phenotypeDistributionPLOT, 11,
        max(5.5, 2.5 + .55 * uniqueN(phenotypeOUTCOME$outcome)),
        "Phenotype", "Outcome review", "Phenotype outcome distributions",
        "Every finite value for each prespecified phenotype outcome, shown on its recorded scale.",
        "Inspect range, clustering, isolated observations, and participant coverage without comparing absolute scales across outcomes.",
        saveFIGUREDATA("phenotype_outcome_distributions", phenotypeOUTCOME)
      )

      outcomePAIR <- CJ(outcome_1 = outcomeCOLUMN, outcome_2 = outcomeCOLUMN)
      phenotypeCORRELATION <- rbindlist(lapply(seq_len(nrow(outcomePAIR)), function(row) {
        first <- suppressWarnings(as.numeric(phenotypeDATA[[outcomePAIR$outcome_1[row]]]))
        second <- suppressWarnings(as.numeric(phenotypeDATA[[outcomePAIR$outcome_2[row]]]))
        complete <- is.finite(first) & is.finite(second)
        data.table(
          outcome_1 = outcomePAIR$outcome_1[row],
          outcome_2 = outcomePAIR$outcome_2[row],
          participants = sum(complete),
          correlation = if (sum(complete) > 1L) stats::cor(first[complete], second[complete]) else NA_real_
        )
      }))
      phenotypeCorrelationPLOT <- ggplot(
        phenotypeCORRELATION,
        aes(x = outcome_1, y = outcome_2, fill = correlation)
      ) +
        geom_tile(colour = "white", linewidth = .6) +
        geom_text(aes(label = paste0(number(correlation, accuracy = .01), "\nn=", participants)), size = 3.2) +
        scale_fill_gradient2(
          low = DNAPRS_COLOURS[["blue"]], mid = "white", high = DNAPRS_COLOURS[["red"]],
          midpoint = 0, limits = c(-1, 1), na.value = "#E5E7E9"
        ) +
        coord_equal() +
        labs(x = NULL, y = NULL, fill = "Pearson r") +
        theme_minimal(base_size = 10) +
        theme(axis.text.x = element_text(angle = 35, hjust = 1), panel.grid = element_blank())
      savePLOT(
        "phenotype_outcome_correlations", phenotypeCorrelationPLOT, 9, 8,
        "Phenotype", "Outcome review", "Phenotype outcome correlations",
        "Pairwise Pearson correlations and complete-case participant counts among prespecified phenotype outcomes.",
        "Inspect shared outcome information and complete-case counts; correlation does not make the outcome scales interchangeable.",
        saveFIGUREDATA("phenotype_outcome_correlations", phenotypeCORRELATION)
      )
    }
  }

  estimated <- association[status == "ESTIMATED" & is.finite(beta)]
  if (nrow(estimated) > 0L) {
    estimated <- addMETHODLABEL(estimated)
    estimated[, model_label := paste(model_id, method_label, sep = " | ")]
    associationPLOT <- ggplot(
      estimated,
      aes(x = beta, y = model_label, colour = method_label)
    ) +
      geom_vline(xintercept = 0, linetype = 2, colour = "grey55") +
      geom_errorbar(
        aes(xmin = ci_low, xmax = ci_high),
        width = .18,
        orientation = "y"
      ) +
      geom_point(size = 2.4) +
      facet_wrap(~cohort, scales = "free_y") +
      scale_colour_manual(values = methodCOLOUR) +
      labs(
        x = "Phenotype coefficient per PRS standard deviation",
        y = NULL,
        colour = "Method"
      ) +
      theme_minimal(base_size = 11) +
      theme(legend.position = "top")
    savePLOT(
      "phenotype_prs_effects",
      associationPLOT,
      10,
      min(14, max(5.5, 2 + .45 * nrow(estimated))),
      "Phenotype",
      "Effect estimates",
      "Phenotype PRS effects",
      "Adjusted phenotype coefficients and 95 percent confidence intervals.",
      "Inspect effect direction, confidence-interval width, and the participant count. Small-sample uncertainty should remain visible.",
      saveFIGUREDATA("phenotype_prs_effects", estimated)
    )

    fitPLOT <- ggplot(
      estimated,
      aes(x = incremental_fit, y = model_label, fill = method_label)
    ) +
      geom_vline(xintercept = 0, linetype = 2, colour = "grey55") +
      geom_col(width = .65) +
      facet_grid(cohort ~ fit_metric, scales = "free_x") +
      scale_fill_manual(values = methodCOLOUR, guide = "none") +
      labs(x = "Improvement after adding PRS", y = NULL) +
      theme_minimal(base_size = 11) +
      theme(panel.grid.major.y = element_blank())
    savePLOT(
      "phenotype_incremental_fit",
      fitPLOT,
      11,
      min(14, max(5.5, 2 + .45 * nrow(estimated))),
      "Phenotype",
      "Added model fit",
      "Phenotype incremental fit",
      "Improvement in phenotype-model fit after adding each matching PRS.",
      "Inspect the size of the added fit together with coefficient uncertainty. A positive change is not independent validation.",
      saveFIGUREDATA("phenotype_incremental_fit", estimated)
    )
  }

  partial <- phenotypePLOT[is.finite(adjusted_outcome) & is.finite(adjusted_prs)]
  if (nrow(partial) > 0L) {
    partial <- addMETHODLABEL(partial)
    partialPLOT <- ggplot(
      partial,
      aes(x = adjusted_prs, y = adjusted_outcome, colour = method_label)
    ) +
      geom_hline(yintercept = 0, colour = "grey85") +
      geom_vline(xintercept = 0, colour = "grey85") +
      geom_point(alpha = .75, size = 2) +
      geom_smooth(method = "lm", formula = y ~ x, se = TRUE, linewidth = .7) +
      facet_wrap(vars(cohort, model_id, method_label), scales = "free", ncol = 2) +
      scale_colour_manual(values = methodCOLOUR, guide = "none") +
      labs(
        x = "PRS after covariate adjustment",
        y = "Phenotype after covariate adjustment"
      ) +
      theme_minimal(base_size = 10)
    partialHEIGHT <- max(
      6,
      ceiling(uniqueN(partial[, .(cohort, model_id, method)]) / 2) * 2.7
    )
    savePLOT(
      "adjusted_phenotype_prs_patterns",
      partialPLOT,
      11,
      partialHEIGHT,
      "Phenotype",
      "Covariate-adjusted pattern",
      "Adjusted phenotype PRS patterns",
      "Covariate-adjusted continuous phenotypes and matching PRSs.",
      "Inspect every observation, the direction of the adjusted pattern, and whether one observation appears to control the fitted line.",
      saveFIGUREDATA("adjusted_phenotype_prs_patterns", partial)
    )
  }

  diagnostic <- phenotypePLOT[is.finite(observed) & is.finite(fitted_full)]
  if (nrow(diagnostic) > 0L) {
    diagnostic <- addMETHODLABEL(diagnostic)
    fittedPLOT <- ggplot(
      diagnostic,
      aes(x = observed, y = fitted_full, colour = method_label)
    ) +
      geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey55") +
      geom_point(alpha = .75, size = 2) +
      facet_wrap(vars(cohort, model_id, method_label), scales = "free", ncol = 2) +
      scale_colour_manual(values = methodCOLOUR, guide = "none") +
      labs(x = "Observed phenotype", y = "Full-model fitted value") +
      theme_minimal(base_size = 10)
    diagnosticHEIGHT <- max(
      6,
      ceiling(uniqueN(diagnostic[, .(cohort, model_id, method)]) / 2) * 2.7
    )
    savePLOT(
      "phenotype_observed_fitted",
      fittedPLOT,
      11,
      diagnosticHEIGHT,
      "Phenotype",
      "Model diagnostics",
      "Observed and fitted phenotype",
      "Observed phenotype values and full-model fitted values.",
      "Inspect agreement with the diagonal and any observations that are fitted poorly.",
      saveFIGUREDATA("phenotype_observed_fitted", diagnostic)
    )
  }

  gaussianDIAGNOSTIC <- phenotypePLOT[
    family == "gaussian" & is.finite(residual_full) & is.finite(fitted_full)
  ]
  if (nrow(gaussianDIAGNOSTIC) > 0L) {
    gaussianDIAGNOSTIC <- addMETHODLABEL(gaussianDIAGNOSTIC)
    residualPLOT <- ggplot(
      gaussianDIAGNOSTIC,
      aes(x = fitted_full, y = residual_full, colour = method_label)
    ) +
      geom_hline(yintercept = 0, linetype = 2, colour = "grey55") +
      geom_point(alpha = .75, size = 2) +
      facet_wrap(vars(cohort, model_id, method_label), scales = "free", ncol = 2) +
      scale_colour_manual(values = methodCOLOUR, guide = "none") +
      labs(x = "Full-model fitted value", y = "Residual") +
      theme_minimal(base_size = 10)
    residualHEIGHT <- max(
      6,
      ceiling(uniqueN(gaussianDIAGNOSTIC[, .(cohort, model_id, method)]) / 2) * 2.7
    )
    savePLOT(
      "phenotype_residual_diagnostics",
      residualPLOT,
      11,
      residualHEIGHT,
      "Phenotype",
      "Model diagnostics",
      "Phenotype residual diagnostics",
      "Gaussian model residuals against full-model fitted values.",
      "Inspect residual symmetry around zero, changing spread, and isolated residuals.",
      saveFIGUREDATA("phenotype_residual_diagnostics", gaussianDIAGNOSTIC)
    )
  }

  permutationDIAGNOSTIC <- phenotypePERMUTATION[
    status == "ESTIMATED" & is.finite(permuted_beta) & is.finite(observed_beta)
  ]
  if (nrow(permutationDIAGNOSTIC) > 0L) {
    permutationDIAGNOSTIC <- addMETHODLABEL(permutationDIAGNOSTIC)
    observedPERMUTATION <- unique(permutationDIAGNOSTIC[, .(
      cohort, model_id, method_label, observed_beta, permutation_scheme
    )])
    permutationPLOT <- ggplot(
      permutationDIAGNOSTIC,
      aes(x = permuted_beta, fill = method_label)
    ) +
      geom_histogram(bins = 50, colour = "white", alpha = .88) +
      geom_vline(
        data = observedPERMUTATION,
        aes(xintercept = observed_beta),
        inherit.aes = FALSE,
        colour = DNAPRS_COLOURS[["red"]],
        linewidth = .75
      ) +
      facet_wrap(vars(cohort, model_id, method_label), scales = "free", ncol = 2) +
      scale_fill_manual(values = methodCOLOUR, guide = "none") +
      labs(x = "PRS coefficient under the null", y = "Permutations") +
      theme_minimal(base_size = 10)
    permutationHEIGHT <- max(
      6,
      ceiling(uniqueN(permutationDIAGNOSTIC[, .(cohort, model_id, method)]) / 2) * 3
    )
    savePLOT(
      "phenotype_residual_permutations", permutationPLOT, 11, permutationHEIGHT,
      "Phenotype", "Robustness diagnostics", "Phenotype residual-permutation distributions",
      "Freedman-Lane null distributions for fixed-effects Gaussian PRS coefficients; red lines are the observed coefficients.",
      "Inspect where each observed coefficient falls in its null distribution and use the reported empirical P value with the parametric estimate.",
      saveFIGUREDATA("phenotype_residual_permutations", permutationDIAGNOSTIC)
    )
  }

  influenceDIAGNOSTIC <- phenotypeINFLUENCE[
    status == "ESTIMATED" & is.finite(beta_without) & is.finite(full_beta)
  ]
  if (nrow(influenceDIAGNOSTIC) > 0L) {
    influenceDIAGNOSTIC <- addMETHODLABEL(influenceDIAGNOSTIC)
    influenceDIAGNOSTIC[, deletion_index := seq_len(.N), by = .(cohort, model_id, method)]
    fullINFLUENCE <- unique(influenceDIAGNOSTIC[, .(
      cohort, model_id, method_label, full_beta
    )])
    influencePLOT <- ggplot(
      influenceDIAGNOSTIC,
      aes(x = deletion_index, y = beta_without, colour = method_label)
    ) +
      geom_hline(
        data = fullINFLUENCE,
        aes(yintercept = full_beta),
        inherit.aes = FALSE,
        colour = DNAPRS_COLOURS[["grey"]],
        linetype = 2
      ) +
      geom_line(alpha = .6, linewidth = .5) +
      geom_point(alpha = .82, size = 1.8) +
      facet_wrap(vars(cohort, model_id, method_label), scales = "free", ncol = 2) +
      scale_colour_manual(values = methodCOLOUR, guide = "none") +
      labs(x = "Participant deletion index", y = "PRS coefficient after deletion") +
      theme_minimal(base_size = 10)
    influenceHEIGHT <- max(
      6,
      ceiling(uniqueN(influenceDIAGNOSTIC[, .(cohort, model_id, method)]) / 2) * 3
    )
    savePLOT(
      "phenotype_leave_one_out", influencePLOT, 11, influenceHEIGHT,
      "Phenotype", "Robustness diagnostics", "Leave-one-participant-out phenotype effects",
      "PRS coefficients after deleting each participant from a fixed-effects Gaussian model; dashed lines show the full-model coefficients.",
      "Inspect abrupt coefficient changes or sign reversals. The downloadable source table identifies the omitted participant.",
      saveFIGUREDATA("phenotype_leave_one_out", influenceDIAGNOSTIC)
    )
  }
}

manifest[, `:=`(
  description = vapply(publish_path, describeFILE, character(1)),
  sensitivity = fifelse(
    grepl("^(scores|phenotype)", publish_path),
    "Participant-level research data",
    "Research workflow output"
  )
)]
manifest[, bytes := file.info(file.path(inputROOT, file_name))$size]
manifest[, sha256 := vapply(
  file.path(inputROOT, file_name),
  function(path) {
    connection <- file(path, "rb")
    on.exit(close(connection), add = TRUE)
    as.character(openssl::sha256(connection))
  },
  character(1)
)]
fwrite(manifest, file.path("provenance", "output_files.tsv"), sep = "\t")

reportVERSION <- data.table(
  software = c("R", "Quarto", "data.table", "ggplot2", "knitr", "openssl"),
  version = c(
    paste(R.version$major, R.version$minor, sep = "."),
    Sys.getenv("QUARTO_VERSION", unset = "Recorded by the execution environment"),
    as.character(packageVersion("data.table")),
    as.character(packageVersion("ggplot2")),
    as.character(packageVersion("knitr")),
    as.character(packageVersion("openssl"))
  )
)
fwrite(
  reportVERSION,
  file.path("provenance", "report_software_versions.tsv"),
  sep = "\t"
)

columnMEANING <- c(
  cohort = "Target cohort identifier.",
  IID = "Participant identifier.",
  trait_id = "GWAS trait identifier.",
  prs_name = "Analysis-ready PRS name.",
  method = "PRS generation method.",
  raw_prs = "Weighted allele or dosage sum.",
  prs_z = "PRS standardised within the declared cohort.",
  used_variants = "Variants contributing to the participant score.",
  participants = "Number of participants.",
  status = "Pipeline or model result status.",
  model_id = "Prespecified phenotype-model identifier.",
  beta = "Estimated matching PRS coefficient.",
  se = "Standard error of the PRS coefficient.",
  ci_low = "Lower 95 percent confidence limit.",
  ci_high = "Upper 95 percent confidence limit.",
  parametric_p = "Parametric P value for the PRS coefficient.",
  permutation_p = "Freedman-Lane residual-permutation P value.",
  permutation_holm = "Permutation P value after Holm correction within cohort and method.",
  r2_base = "R-squared for the model without PRS.",
  r2_full = "R-squared for the model with PRS.",
  delta_r2 = "Difference between full-model and base-model R-squared.",
  partial_r2 = "Delta R-squared divided by one minus base-model R-squared.",
  p_value = "PRS coefficient P value.",
  incremental_fit = "Change in the declared model-fit measure after adding PRS.",
  fit_metric = "Model-fit measure used for incremental_fit.",
  stage = "PRS generation stage.",
  variant_count = "Variants retained at the stage.",
  percent_of_source = "Percentage of source GWAS variants retained.",
  path = "Source input path.",
  checksum = "Cryptographic input checksum.",
  algorithm = "Checksum algorithm."
)

dictionaryLIST <- lapply(seq_len(nrow(manifest)), function(row) {
  path <- manifest$relative_path[row]
  if (!grepl("\\.(tsv|csv|ma|txt)(\\.gz)?$", path, ignore.case = TRUE)) {
    return(NULL)
  }
  value <- tryCatch(
    fread(path, nrows = 100L, showProgress = FALSE),
    error = function(error) NULL
  )
  if (is.null(value) || !length(names(value))) return(NULL)
  data.table(
    result_file = manifest$file_name[row],
    column = names(value),
    data_type = vapply(value, function(column) class(column)[1L], character(1)),
    meaning = vapply(names(value), function(column) {
      description <- unname(columnMEANING[column])
      if (is.na(description)) {
        "Field supplied by the named pipeline result."
      } else {
        description
      }
    }, character(1)),
    missing_value = "NA"
  )
})
dataDICTIONARY <- rbindlist(dictionaryLIST, fill = TRUE)
fwrite(dataDICTIONARY, file.path("provenance", "data_dictionary.tsv"), sep = "\t")
fwrite(figureMANIFEST, file.path("provenance", "figure_manifest.tsv"), sep = "\t")

# Build one workbook containing every moderate report table and exact figure-source table.
workbook <- createWorkbook()
addWorksheet(workbook, "Contents")
workbookCONTENT <- data.table(
  worksheet = character(),
  report_page = character(),
  table_title = character(),
  description = character(),
  source_process = character(),
  rows = integer(),
  columns = integer(),
  sensitivity = character()
)
workbookDICTIONARY <- data.table(
  worksheet = character(),
  column = character(),
  definition = character(),
  data_type = character(),
  units_or_scale = character(),
  expected_values = character(),
  missing_value = character(),
  source_stage = character(),
  sensitivity = character()
)
usedSHEETS <- "Contents"

safeSHEET <- function(value) {
  value <- gsub("[\\\\/:?*\\[\\]]", "_", value)
  value <- substr(value, 1L, 31L)
  candidate <- value
  suffix <- 1L
  while (candidate %in% usedSHEETS) {
    suffix <- suffix + 1L
    candidate <- paste0(substr(value, 1L, 28L), "_", suffix)
  }
  usedSHEETS <<- c(usedSHEETS, candidate)
  candidate
}

addWORKBOOKTABLE <- function(sheet, value, page, title, description, source, sensitivity) {
  if (is.null(value) || !ncol(value) || !nrow(value)) return(invisible(NULL))
  if (nrow(value) > 1048575L || ncol(value) > 16384L) {
    workbookCONTENT <<- rbind(
      workbookCONTENT,
      data.table(
        worksheet = "Native download only", report_page = page, table_title = title,
        description = paste(description, "The table exceeds an Excel worksheet limit."),
        source_process = source, rows = nrow(value), columns = ncol(value),
        sensitivity = sensitivity
      )
    )
    return(invisible(NULL))
  }
  sheet <- safeSHEET(sheet)
  addWorksheet(workbook, sheet)
  writeDataTable(workbook, sheet, as.data.frame(value), withFilter = TRUE, tableStyle = "TableStyleLight9")
  freezePane(workbook, sheet, firstRow = TRUE)
  setColWidths(workbook, sheet, cols = seq_len(ncol(value)), widths = 18)
  narrativeCOLUMN <- which(names(value) %in% c(
    "description", "reason", "definition", "meaning", "missing_value",
    "expected_values", "source_genotype", "source_sample", "source_keep"
  ))
  if (length(narrativeCOLUMN)) setColWidths(workbook, sheet, cols = narrativeCOLUMN, widths = 42)
  workbookCONTENT <<- rbind(
    workbookCONTENT,
    data.table(
      worksheet = sheet, report_page = page, table_title = title,
      description = description, source_process = source, rows = nrow(value),
      columns = ncol(value), sensitivity = sensitivity
    )
  )
  workbookDICTIONARY <<- rbind(
    workbookDICTIONARY,
    data.table(
      worksheet = sheet,
      column = names(value),
      definition = vapply(names(value), function(column) {
        meaning <- unname(columnMEANING[column])
        if (is.na(meaning)) paste("Field supplied by", title) else meaning
      }, character(1)),
      data_type = vapply(value, function(column) class(column)[1L], character(1)),
      units_or_scale = "See definition and source result",
      expected_values = "Defined by the source stage",
      missing_value = "NA means unavailable or not applicable unless the source definition states otherwise",
      source_stage = source,
      sensitivity = sensitivity
    ),
    fill = TRUE
  )
  invisible(sheet)
}

addWORKBOOKTABLE("Run", inputCHECK, "Overview", "Run checks", "Validated run-level checks.", "Input checks", "Research workflow output")
addWORKBOOKTABLE("Target_inputs", targetMANIFEST, "Overview", "Target inputs", "Validated generated target records.", "Input checks", "Restricted input identity")
addWORKBOOKTABLE("GWAS_inputs", gwasMANIFEST, "Overview", "GWAS inputs", "Validated generated GWAS records.", "Input checks", "Research workflow output")
addWORKBOOKTABLE("Reference_inputs", referenceMANIFEST, "Overview", "Reference resources", "Validated generated reference records.", "Input checks", "Research workflow output")
addWORKBOOKTABLE("Phenotype_models", phenotypeMANIFEST, "Overview", "Phenotype models", "Prespecified generated phenotype model records.", "Input checks", "Restricted analysis specification")
addWORKBOOKTABLE("Genotype_EDA", genotypeEDA, "Genotype EDA", "Target composition", "Supplied target dimensions and composition at the declared input stage.", "Genotype EDA", "Restricted cohort summary")
addWORKBOOKTABLE("Genotype_EDA_checks", genotypeEDACHECK, "Genotype EDA", "EDA checks", "Exploratory status and reasons.", "Genotype EDA", "Restricted cohort summary")
addWORKBOOKTABLE("Target_PREP", targetPREP, "Target PREP", "Preparation summary", "Target marker preparation summary.", "Target PREP", "Restricted cohort summary")
addWORKBOOKTABLE("Target_QC", targetQC, "Target QC", "QC summary", "Target scoring-readiness summary.", "Target QC", "Restricted cohort summary")
addWORKBOOKTABLE("Target_sample_decisions", targetSAMPLEDECISION, "Target QC", "Participant decisions", "Participant-level technical QC decisions.", "Target QC", "Restricted participant result")
addWORKBOOKTABLE("Target_variant_decisions", targetVARIANTDECISION, "Target QC", "Variant decisions", "Variant-level technical QC decisions.", "Target QC", "Restricted cohort summary")
addWORKBOOKTABLE("Participant_decisions", participantDECISION, "Target QC", "Integrated participant decisions", "Technical-QC, relatedness, reference-ancestry, score-eligibility, and primary-analysis decisions.", "Target QC", "Restricted participant result")
addWORKBOOKTABLE("Target_ancestry", targetANCESTRY, "Target QC", "Target ancestry classification", "Reference-projected target PCs, distance, threshold, and classification.", "Target QC", "Restricted participant result")
addWORKBOOKTABLE("Reference_projection", referencePROJECTION, "Target QC", "Reference ancestry projection", "Population-reference projection used to classify target participants.", "Target QC", "Restricted reference result")
addWORKBOOKTABLE("Ancestry_summary", ancestrySUMMARY, "Target QC", "Ancestry summary", "Variant, participant, and classification counts for reference-anchored ancestry QC.", "Target QC", "Restricted cohort summary")
addWORKBOOKTABLE("Imputation_QC", imputationQC, "Target Imputation", "Imputation QC", "Per-chromosome imputation counts and DR2 threshold.", "Target imputation", "Restricted cohort summary")
addWORKBOOKTABLE("Imputation_manifest", imputationMANIFEST, "Target Imputation", "Imputation handoff", "Reusable chromosome-level imputed target manifest.", "Target imputation", "Restricted input identity")
addWORKBOOKTABLE("Imputation_DR2", imputationDR2, "Target Imputation", "Imputation DR2 bins", "Per-chromosome imputed-variant counts by Beagle dosage R-squared bin.", "Target imputation", "Restricted cohort summary")
addWORKBOOKTABLE("PLINK_sensitivity", plinkSENSITIVITY, "PLINK PRS", "Imputation sensitivity", "Participant-level comparison of imputed-target and direct-genotype PLINK scores.", "PLINK scoring", "Restricted participant result")
addWORKBOOKTABLE("PLINK_sensitivity_QC", plinkSENSITIVITYQC, "PLINK PRS", "Sensitivity and coverage summary", "Typed/imputed variant coverage and score agreement.", "PLINK scoring", "Restricted cohort summary")
addWORKBOOKTABLE("GWAS_QC", gwasSUMMARY, "GWAS QC", "GWAS QC summary", "GWAS source and retained counts.", "GWAS QC", "Research workflow output")
addWORKBOOKTABLE("Variant_flow", variantFLOW, "Overview", "Variant flow", "Variant counts through each PRS method.", "PRS generation", "Research workflow output")
addWORKBOOKTABLE("Score_QC", scoreQC, "PLINK PRS / SBayesRC PRS", "Score QC", "Participant and score coverage checks.", "PRS generation", "Restricted cohort summary")
addWORKBOOKTABLE("Associations", association, "Phenotype", "Phenotype associations", "Prespecified phenotype association estimates.", "Phenotype", "Restricted research result")
addWORKBOOKTABLE("Fitted_models", fittedMODEL, "Phenotype", "Fitted model specifications", "Exact fitted formulas and estimators.", "Phenotype", "Restricted analysis specification")
addWORKBOOKTABLE("Permutations", phenotypePERMUTATION, "Phenotype", "Residual permutations", "Freedman-Lane coefficient null distributions and explicit unsupported-model records.", "Phenotype", "Restricted research result")
addWORKBOOKTABLE("Influence", phenotypeINFLUENCE, "Phenotype", "Case-deletion diagnostics", "Participant-level Gaussian-model case-deletion coefficients and explicit unsupported-model records.", "Phenotype", "Restricted research result")

figureSOURCE <- unique(figureMANIFEST[nzchar(source_table), .(figure_id, page, title, description, source_table)])
if (nrow(figureSOURCE)) {
  for (row in seq_len(nrow(figureSOURCE))) {
    sourcePATH <- figureSOURCE$source_table[row]
    sourceVALUE <- if (file.exists(sourcePATH)) fread(sourcePATH, showProgress = FALSE) else data.table()
    addWORKBOOKTABLE(
      sprintf("Fig_%03d", row), sourceVALUE, figureSOURCE$page[row],
      paste("Figure source:", figureSOURCE$title[row]), figureSOURCE$description[row],
      "Report assembly", "Matches the source pipeline result"
    )
  }
}

addWorksheet(workbook, "Dictionary")
usedSHEETS <- c(usedSHEETS, "Dictionary")
writeDataTable(workbook, "Dictionary", as.data.frame(workbookDICTIONARY), withFilter = TRUE, tableStyle = "TableStyleLight9")
freezePane(workbook, "Dictionary", firstRow = TRUE)
setColWidths(workbook, "Dictionary", cols = seq_len(ncol(workbookDICTIONARY)), widths = 22)
setColWidths(workbook, "Dictionary", cols = c(3L, 5L, 6L, 7L), widths = 42)
workbookCONTENT <- rbind(
  data.table(
    worksheet = "Contents", report_page = "Logs", table_title = "Workbook contents",
    description = "Index of every worksheet and sensitivity classification.",
    source_process = "Report assembly", rows = nrow(workbookCONTENT) + 2L,
    columns = 8L, sensitivity = "Research workflow output"
  ),
  workbookCONTENT,
  data.table(
    worksheet = "Dictionary", report_page = "Logs", table_title = "Column dictionary",
    description = "Plain-language definition and provenance for every workbook column.",
    source_process = "Report assembly", rows = nrow(workbookDICTIONARY),
    columns = ncol(workbookDICTIONARY), sensitivity = "Research workflow output"
  )
)
writeDataTable(workbook, "Contents", as.data.frame(workbookCONTENT), withFilter = TRUE, tableStyle = "TableStyleLight9")
freezePane(workbook, "Contents", firstRow = TRUE)
setColWidths(workbook, "Contents", cols = seq_len(ncol(workbookCONTENT)), widths = 20)
setColWidths(workbook, "Contents", cols = c(3L, 4L), widths = 42)
saveWorkbook(workbook, file.path("downloads", "dnaprs_report_tables.xlsx"), overwrite = TRUE)

reportSTATE <- data.table(
  setting = c("has_plink", "has_sbayesrc", "has_phenotype"),
  value = c(hasPLINK, hasSBAYESRC, hasPHENOTYPE)
)
fwrite(reportSTATE, file.path("provenance", "report_state.tsv"), sep = "\t")
