#!/usr/bin/env Rscript

suppressWarnings(suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
}))

options(stringsAsFactors = FALSE, scipen = 999)

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
#   Save one quantitative plot in the three report formats and record its meaning.
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
#   The three figure files and one row added to figureMANIFEST.
savePLOT <- function(id, plot, width, height, page, section, title,
                     description, inspection) {
  pathTIFF <- file.path("figures", "tiff", paste0(id, ".tiff"))
  pathPNG <- file.path("figures", "png", paste0(id, ".png"))
  pathJPEG <- file.path("figures", "jpeg", paste0(id, ".jpeg"))
  ggsave(pathTIFF, plot, width = width, height = height, dpi = 300,
         compression = "lzw", bg = "white")
  ggsave(pathPNG, plot, width = width, height = height, dpi = 300,
         bg = "white")
  ggsave(pathJPEG, plot, width = width, height = height, dpi = 300,
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
      dpi = 300L,
      tiff = gsub("\\\\", "/", pathTIFF),
      png = gsub("\\\\", "/", pathPNG),
      jpeg = gsub("\\\\", "/", pathJPEG)
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
    inputs = "inputs.qmd",
    target_gwas = "target-gwas.qmd",
    plink = "plink.qmd",
    sbayesrc = "sbayesrc.qmd",
    prs = "prs.qmd",
    phenotype = "phenotype.qmd",
    dictionary = "dictionary.qmd",
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

targetMANIFEST <- readRESULT("target_manifest.validated.tsv")
gwasMANIFEST <- readRESULT("gwas_manifest.validated.tsv")
referenceMANIFEST <- readRESULT("reference_manifest.validated.tsv")
phenotypeMANIFEST <- readRESULT("phenotype_models.validated.tsv")
score <- readRESULT("prs_scores_long.tsv")
scoreQC <- readRESULT("score_qc.tsv")
concordance <- readRESULT("method_concordance.tsv")
variantFLOW <- readRESULT("variant_flow.tsv")
association <- readRESULT("phenotype_associations.tsv")
phenotypePLOT <- readRESULT("phenotype_plot_data.tsv")

methodLABEL <- c(plink_ct = "PLINK C+T", sbayesrc = "SBayesRC")
methodCOLOUR <- c("PLINK C+T" = "#111111", "SBayesRC" = "#777777")
score <- addMETHODLABEL(score)
scoreQC <- addMETHODLABEL(scoreQC)
variantFLOW <- addMETHODLABEL(variantFLOW)

hasPLINK <- nrow(score) > 0L && any(score$method == "plink_ct")
hasSBAYESRC <- nrow(score) > 0L && any(score$method == "sbayesrc")
hasPHENOTYPE <- nrow(association) > 0L
pageVALUE <- c(
  overview = "Overview",
  inputs = "Inputs",
  target_gwas = "Target and GWAS"
)
if (hasPLINK) pageVALUE <- c(pageVALUE, plink = "PLINK")
if (hasSBAYESRC) pageVALUE <- c(pageVALUE, sbayesrc = "SBayesRC")
pageVALUE <- c(pageVALUE, prs = "PRS")
if (hasPHENOTYPE) pageVALUE <- c(pageVALUE, phenotype = "Phenotype")
pageVALUE <- c(pageVALUE, dictionary = "Dictionary", logs = "Logs")
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
  tiff = character(),
  png = character(),
  jpeg = character()
)

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
  "Target and GWAS",
  "Variant flow",
  "Variant generation flow",
  "Variant counts across GWAS harmonisation, method-specific modelling, and final scoring.",
  "Inspect large unexpected losses, differences between traits, and the final number of variants that contributed to scoring."
)

coveragePLOT <- ggplot(scoreQC, aes(x = prs_name, y = used_variants, fill = method_label)) +
  geom_col(width = .68, na.rm = TRUE) +
  geom_text(aes(label = comma(used_variants)), hjust = -.08, size = 3.2, na.rm = TRUE) +
  coord_flip(clip = "off") +
  facet_grid(cohort ~ method_label, scales = "free_x") +
  scale_fill_manual(values = methodCOLOUR, guide = "none") +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, .24))) +
  labs(x = NULL, y = "Variants contributing to the score") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.major.y = element_blank())
savePLOT(
  "scoring_variant_counts",
  coveragePLOT,
  10,
  max(5, uniqueN(scoreQC$cohort) * 2.6),
  "PRS",
  "Variants used",
  "Scoring variant counts",
  "Variants contributing to each participant score.",
  "Compare coverage across traits, cohorts, and methods. A low count should be interpreted with the recorded variant-flow and score-QC tables."
)

distributionPLOT <- ggplot(score, aes(x = prs_z, y = method_label, colour = method_label)) +
  geom_vline(xintercept = 0, colour = "grey75", linewidth = .4) +
  geom_boxplot(
    aes(group = method_label),
    width = .32,
    outlier.shape = NA,
    alpha = .12,
    colour = "grey45"
  ) +
  geom_jitter(height = .12, width = 0, alpha = .75, size = 2) +
  facet_grid(cohort ~ prs_name, scales = "free_y", space = "free_y") +
  scale_colour_manual(values = methodCOLOUR) +
  labs(x = "PRS (within-cohort standard deviations)", y = NULL, colour = NULL) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "top", panel.grid.major.y = element_blank())
savePLOT(
  "participant_prs_distributions",
  distributionPLOT,
  11,
  max(6, uniqueN(score$cohort) * 3.2),
  "PRS",
  "Participant distributions",
  "Participant PRS distributions",
  "Participant-level standardised PRS values.",
  "Inspect every participant for missing values, isolated observations, or method-specific compression of the score range."
)

score[, score_label := paste(prs_name, method_label, sep = " | ")]
correlationLIST <- lapply(unique(score$cohort), function(cohortVALUE) {
  wide <- dcast(score[cohort == cohortVALUE], IID ~ score_label, value.var = "prs_z")
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
  scale_fill_gradient2(
    low = "#202020",
    mid = "white",
    high = "#8a8a8a",
    midpoint = 0,
    limits = c(-1, 1)
  ) +
  coord_equal() +
  labs(x = NULL, y = NULL, fill = "Pearson r") +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 38, hjust = 1), panel.grid = element_blank())
savePLOT(
  "cross_score_correlations",
  correlationPLOT,
  10,
  8,
  "PRS",
  "Cross-score correlation",
  "Cross-score correlations",
  "Pearson correlations among standardised trait and method scores.",
  "Inspect shared information among scores. Correlation describes similarity and does not by itself validate a PRS."
)

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
    "PRS",
    "Method agreement",
    "Method agreement",
    "PLINK C+T and SBayesRC participant scores.",
    "Inspect whether the two methods rank participants similarly and whether any participant departs strongly from the general relationship."
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
    "PRS",
    "Method agreement",
    "Participant method ranks",
    "Participant rank under each PRS method.",
    "Inspect large rank changes between methods. Rank agreement supports consistency but does not establish prediction."
  )
}

if (hasPHENOTYPE) {
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
      "Inspect effect direction, confidence-interval width, and the participant count. Small-sample uncertainty should remain visible."
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
      "Inspect the size of the added fit together with coefficient uncertainty. A positive change is not independent validation."
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
      "Inspect every observation, the direction of the adjusted pattern, and whether one observation appears to control the fitted line."
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
      "Inspect agreement with the diagonal and any observations that are fitted poorly."
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
      "Inspect residual symmetry around zero, changing spread, and isolated residuals."
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

reportSTATE <- data.table(
  setting = c("has_plink", "has_sbayesrc", "has_phenotype"),
  value = c(hasPLINK, hasSBAYESRC, hasPHENOTYPE)
)
fwrite(reportSTATE, file.path("provenance", "report_state.tsv"), sep = "\t")
