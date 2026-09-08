#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))
manifest <- fread("provenance/figure_manifest.tsv")
stopifnot(all(file.info(unlist(manifest[, .(svg, tiff, png, jpeg)]))$size > 0))
for (row in seq_len(nrow(manifest))) {
  value <- fread(manifest$source_table[row])
  stopifnot(nrow(value) == manifest$source_rows[row], ncol(value) == manifest$source_columns[row])
  if (manifest$figure_id[row] %in% c("gwas_qq", "gwas_manhattan")) {
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
stopifnot(all(selection[figure_id == "gwas_qq", displayed_records] <= 51000L))
stopifnot(all(selection$displayed_records <= selection$eligible_records))
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
