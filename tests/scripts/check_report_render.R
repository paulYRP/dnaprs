#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))
manifest <- fread("provenance/figure_manifest.tsv")
stopifnot(all(file.info(unlist(manifest[, .(svg, tiff, png, jpeg)]))$size > 0))
for (row in seq_len(nrow(manifest))) {
  value <- fread(manifest$source_table[row])
  stopifnot(nrow(value) == manifest$source_rows[row], ncol(value) == manifest$source_columns[row])
  if (grepl("^gwas_(qq|manhattan)_", manifest$figure_id[row])) {
    stopifnot(uniqueN(value$trait_id) == 1L, file.info(manifest$preview[row])$size > 0)
    figureTYPE <- if (startsWith(manifest$figure_id[row], "gwas_qq_")) "gwas_qq" else "gwas_manhattan"
    selectionROW <- fread("provenance/gwas_plot_selection.tsv")[trait_id == value$trait_id[1L] & figure_id == figureTYPE]
    stopifnot(nrow(selectionROW) == 1L, selectionROW$displayed_records == nrow(value))
    # Count vector points in bounded chunks, including overlapping strong signals.
    connection <- file(manifest$svg[row], "r")
    points <- 0
    repeat {
      lines <- readLines(connection, n = 10000L, warn = FALSE)
      if (!length(lines)) break
      points <- points + sum(grepl("<circle ", lines, fixed = TRUE))
    }
    close(connection)
    stopifnot(points == nrow(value))
  }
}
selection <- fread("provenance/gwas_plot_selection.tsv")
stopifnot(!anyDuplicated(manifest$figure_id), !any(manifest$figure_id %in% c("gwas_qq", "gwas_manhattan")))
stopifnot(sum(grepl("^gwas_(qq|manhattan)_", manifest$figure_id)) == sum(selection$displayed_records > 0))
stopifnot(all(selection[figure_id == "gwas_qq", displayed_records] <= 51000L))
stopifnot(all(selection$displayed_records <= selection$eligible_records))
eligibility <- fread("data/participant_analysis_eligibility.tsv")
stopifnot(eligibility[analysis_set == "Sensitivity only", N] == 2L,
          eligibility[analysis_set == "Unknown eligibility", N] == 2L,
          nrow(eligibility) == 4L, sum(eligibility$N) == 6L)
stopifnot(identical(unname(tools::md5sum("inputs/phenoPRS.csv")),
                    unname(tools::md5sum("_site/downloads/phenotype/phenoPRS.csv"))))
stopifnot(file.info("_site/index.html")$size > 0, file.info("_site/gwas-qc.html")$size > 0)
stopifnot(file.info("_site/downloads/dnaprs_report_tables.xlsx")$size > 0)
contents <- openxlsx::read.xlsx("_site/downloads/dnaprs_report_tables.xlsx", sheet = "Contents")
large <- manifest[source_rows > 1048575L | source_columns > 16384L]
for (title in large$title) {
  stopifnot(any(contents$worksheet == "Native download only" & contents$table_title == paste("Figure source:", title)))
}
message("Rendered report passed: figures, dimensions, workbook, downloads and unchanged PRS CSV.")

for (page in c("target-prep", "target-imputation", "phenotype")) {
  html <- paste(readLines(file.path("_site", paste0(page, ".html")), warn = FALSE), collapse = "\n")
  stopifnot(grepl("data-paged-table", html, fixed = TRUE), !grepl("```", html, fixed = TRUE))
  stopifnot(!grepl("Browse complete table", html, fixed = TRUE), !grepl("Download complete", html, fixed = TRUE))
  if (page != "phenotype") {
    stage <- if (page == "target-prep") "corrected/SYNTHETIC" else "imputed/SYNTHETIC.imputed"
    for (extension in c("pvar", "pgen", "psam")) {
      stopifnot(grepl(paste0("../data/checkpoints/", stage, "/SYNTHETIC.", extension), html, fixed = TRUE))
      stopifnot(length(regmatches(html, gregexpr(paste0('download="SYNTHETIC.', extension, '"'), html, fixed = TRUE))[[1L]]) == 1L,
        grepl(paste0("Download SYNTHETIC.", extension), html, fixed = TRUE))
    }
  } else stopifnot(length(regmatches(html, gregexpr('download="phenoPRS.csv"', html, fixed = TRUE))[[1L]]) == 1L)
}
stopifnot(all(c("target_preparation_summary", "target_preparation_exclusions") %in% manifest$figure_id))
headings <- list(index = c("Summary", "Runs", "Scores"), `genotype-eda` = c("Target", "Runs", "Results"),
  `target-prep` = c("Figures", "Target", "Results"), `target-qc` = c("Figures", "Participants", "Ancestry", "Variants", "Results"),
  `target-imputation` = c("Figures and Tables", "Target", "Results"), `gwas-qc` = c("Traits", "Figures", "Results"),
  `plink-prs` = c("Variants", "Figures", "Imputation", "Results"),
  `sbayesrc-prs` = c("Alignment and Imputation", "Figures", "Results"),
  phenotype = c("Figures", "Scores", "Estimates", "Models", "Results"), logs = c("Summary", "Results", "Nextflow"))
for (page in names(headings)) {
  path <- file.path("_site", paste0(page, ".html"))
  if (!file.exists(path)) next
  html <- paste(readLines(path, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  found <- regmatches(html, gregexpr('(?s)<h[1-6][^>]*>.*?</h[1-6]>', html, perl = TRUE))[[1L]]
  found <- trimws(gsub('<[^>]+>', '', gsub('<a\\b[^>]*>.*?</a>', '', found, perl = TRUE)))
  for (heading in headings[[page]]) {
    if (!heading %in% found) stop(sprintf("Missing heading '%s' in %s", heading, page))
  }
  anchors <- regmatches(html, gregexpr('<a\\s[^>]*>', html, perl = TRUE))[[1L]]
  downloads <- anchors[grepl(' download="[^"]+"', anchors) & grepl(' href="', anchors)]
  for (link in downloads) {
    href <- utils::URLdecode(sub('[?].*$', '', gsub('&amp;', '&', sub('.* href="([^"]+)".*', '\\1', link), fixed = TRUE)))
    filename <- gsub('&amp;', '&', sub('.* download="([^"]+)".*', '\\1', link), fixed = TRUE)
    expected <- c(input_checks.tsv = "runs.tsv", score_qc.tsv = "scores.tsv")
    expected <- if (page == "index" && basename(href) %in% names(expected)) unname(expected[basename(href)]) else basename(href)
    stopifnot(file.exists(file.path("_site", href)), identical(expected, filename))
  }
}
overviewHTML <- paste(readLines("_site/index.html", warn = FALSE), collapse = "\n")
stopifnot(!grepl("PRS Report", overviewHTML, fixed = TRUE),
          !grepl("Analysis summary", overviewHTML, fixed = TRUE),
          !grepl("dnaprs-summary-grid", overviewHTML, fixed = TRUE),
          grepl('download="runs.tsv">Download runs.tsv', overviewHTML, fixed = TRUE),
          grepl('download="scores.tsv">Download scores.tsv', overviewHTML, fixed = TRUE))
gwasHTML <- paste(readLines("_site/gwas-qc.html", warn = FALSE), collapse = "\n")
stopifnot(!grepl("Discovery studies", gwasHTML, fixed = TRUE),
          grepl("Download TSV", gwasHTML, fixed = TRUE))
imputedHTML <- paste(readLines("_site/target-imputation.html", warn = FALSE), collapse = "\n")
stopifnot(!grepl("The PVAR describes the final filtered imputed variants", imputedHTML, fixed = TRUE))
stopifnot(!file.exists("_site/assets/serve-report.py"), !file.exists("_site/assets/open-report.cmd"),
  file.exists("_site/assets/download-index.js"), file.exists("_site/assets/dnaprs-downloads.js"))
workbookPATH <- "_site/downloads/dnaprs_report_tables.xlsx"
sheets <- openxlsx::getSheetNames(workbookPATH)
stopifnot(identical(sheets[1:2], c("Contents", "Dictionary")))
dictionary <- openxlsx::read.xlsx(workbookPATH, sheet = "Dictionary")
for (sheet in sheets) {
  columns <- names(openxlsx::read.xlsx(workbookPATH, sheet = sheet, rows = 1L, check.names = FALSE))
  stopifnot(setequal(columns, dictionary$column[dictionary$worksheet == sheet]))
}
for (i in seq_along(sheets)) local({
  connection <- unz(workbookPATH, sprintf("xl/worksheets/sheet%d.xml", i), open = "rb")
  on.exit(close(connection))
  tail <- ""
  repeat {
    bytes <- readBin(connection, "raw", n = 65536L)
    if (!length(bytes)) break
    block <- paste0(tail, rawToChar(bytes))
    stopifnot(!grepl("<(autoFilter|pane|cols)\\b", block, perl = TRUE, useBytes = TRUE))
    tail <- rawToChar(tail(charToRaw(block), 128L))
  }
})
rawHTML <- paste(readLines("_site/genotype-eda.html", warn = FALSE), collapse = "\n")
correctedHTML <- paste(readLines("_site/target-qc.html", warn = FALSE), collapse = "\n")
stopifnot(grepl("Invalid raw marker", rawHTML, fixed = TRUE),
          !grepl("Valid corrected calculation", rawHTML, fixed = TRUE),
          grepl("Valid corrected calculation", correctedHTML, fixed = TRUE),
          !grepl("Invalid raw marker", correctedHTML, fixed = TRUE))
message("Rendered report passed: stage separation, native checkpoint viewers and phenotype download links.")
