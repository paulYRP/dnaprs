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
