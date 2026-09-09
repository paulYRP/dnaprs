#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(bit64)
  library(ggplot2)
})
source("assets/report/report-data.R")

testREPORTDATA <- function() {
  decisions <- data.table(cohort = "TEST", primary_analysis = c(TRUE, FALSE, FALSE, FALSE, NA, TRUE),
                          score_eligible = c(TRUE, TRUE, TRUE, FALSE, TRUE, FALSE),
                          reason = c("Pass", "Related", "Ancestry", "QC", "Missing", "Conflicting"))
  originalDECISIONS <- copy(decisions)
  eligibility <- reportELIGIBILITY(decisions)
  stopifnot(identical(decisions, originalDECISIONS), sum(eligibility$totals$N) == 6L,
            eligibility$totals[analysis_set == "Primary analysis", N] == 1L,
            eligibility$totals[analysis_set == "Sensitivity only", N] == 2L,
            eligibility$totals[analysis_set == "Not eligible for scoring", N] == 1L,
            eligibility$totals[analysis_set == "Unknown eligibility", N] == 2L,
            sum(eligibility$reasons[analysis_set == "Sensitivity only", N]) == 2L)
  stopifnot(reportELIGIBILITY(data.table(cohort = "TEST"))$totals$analysis_set == "Unknown eligibility")
  traits <- c("MDD", "mdd", "A/B", "A_B", "A B", "trait.with.dots", "é")
  ids <- vapply(traits, reportTRAITID, character(1))
  stopifnot(!anyDuplicated(ids), all(grepl("^[0-9a-f]+$", ids)),
            identical(ids, vapply(traits, reportTRAITID, character(1))))
  for (n in c(0L, 1L, 9L, 10L, 11L, 100L, 60000L, 1000000L)) {
    p <- if (n) seq(1e-12, 1, length.out = n) else numeric()
    result <- reportQQ(c(p, NA, Inf, -1, 0, 2))
    stopifnot(nrow(result) <= min(n, 51000L))
    if (n) {
      ranks <- sort(unique(c(round(seq(1, n, length.out = min(50000L, n))), seq_len(min(1000L, n)))))
      stopifnot(isTRUE(all.equal(result$expected, -log10(ppoints(n)[ranks]))))
      stopifnot(identical(result$observed, -log10(p[ranks])))
    }
  }
  stopifnot(identical(reportQQ(rep(.1, 60000L)), reportQQ(rep(.1, 60000L))))

  n <- 310000L
  variant <- data.table(SNP = sprintf("rs%07d", seq_len(n)), CHR = rep(1:2, each = n / 2),
                        BP = seq_len(n), p = rep(.5, n))
  variant[1:1200, p := 1e-9]
  layout <- data.table(chromosome = 1:2, chromosome_offset = c(0, n / 2))
  result <- reportMANHATTAN(variant, layout)
  stopifnot(result$full == n, result$strong == 1200L, nrow(result$data) == 251200L)
  stopifnot(all(variant$SNP[1:1200] %in% result$data$SNP))
  stopifnot(identical(result$data, reportMANHATTAN(variant[n:1], layout)$data))
  stopifnot(all(result$data[chromosome == 2, genomic_position] == result$data[chromosome == 2, BP] + n / 2))

  for (x in list(c(-2, -1, 0, 0, .01, .5, 1, 2), rep(.3, 20), seq(0, .5, .005), c(0, .01 - 1e-10, .01, .01 + 1e-10, .02))) {
    for (fixed in c(FALSE, TRUE)) {
      histogram <- if (fixed) geom_histogram(binwidth = .01, boundary = 0) else geom_histogram(bins = 50)
      expected <- as.data.table(ggplot_build(ggplot(data.table(x), aes(x)) + histogram)$data[[1]])
      observed <- if (fixed) reportHISTOGRAM(x, range(x), binwidth = .01, boundary = 0) else reportHISTOGRAM(x, range(x))
      stopifnot(isTRUE(all.equal(as.numeric(observed$count), as.numeric(expected$count))))
      stopifnot(isTRUE(all.equal(observed$xmin, expected$xmin)), isTRUE(all.equal(observed$xmax, expected$xmax)))
    }
  }
  stopifnot(nrow(reportHISTOGRAM(numeric(), c(Inf, -Inf))) == 0L)
  x <- rep(.3, 20)
  expected <- ggplot_build(ggplot(data.table(x), aes(x)) + geom_histogram(bins = 50) + geom_vline(xintercept = 0))$data[[1]]
  observed <- reportHISTOGRAM(x, range(c(0, x)))
  stopifnot(identical(as.numeric(observed$count), as.numeric(expected$count)))
  stopifnot(reportWORKBOOKFITS(1048575L, 16384L), !reportWORKBOOKFITS(1048576L, 2L), !reportWORKBOOKFITS(1L, 16385L))

  directory <- tempfile("report-data-")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  integerPATH <- file.path(directory, "large-integers.tsv")
  original <- data.table(id = as.integer64(c("9007199254740993", "9223372036854775806")))
  fwrite(original, integerPATH, sep = "\t")
  loaded <- fread(integerPATH)
  stopifnot(identical(as.character(loaded$id), as.character(original$id)))
  excelPATH <- file.path(directory, "large-integers.xlsx")
  openxlsx::write.xlsx(reportEXCELDATA(loaded), excelPATH)
  stopifnot(identical(openxlsx::read.xlsx(excelPATH)$id, as.character(original$id)))
  stopifnot(inherits(loaded$id, "integer64"))

  paths <- file.path(directory, c("A.cojo.ma", "B.cojo.ma"))
  for (i in 1:2) fwrite(data.table(SNP = paste0("rs", 1:100), CHR = 1L, BP = (1:100) * i,
                                  b = seq(-i, i, length.out = 100), freq = .3, p = seq(1e-9, .9, length.out = 100)), paths[i])
  prepared <- prepareGWASPLOTS(paths)
  stopifnot(nrow(prepared$qq) == 200L, prepared$layout$chromosome_length == 200L)
  stopifnot(all(prepared$selection$eligible_records == prepared$selection$displayed_records))
  stopifnot(sum(prepared$effect$count) == 200L, sum(prepared$maf$count) == 200L)
  weightPATH <- file.path(directory, "trait.with.dots.sbayesrc.txt")
  fwrite(data.table(BETA = c(-.1, .1)), weightPATH)
  weightBINS <- prepareREPORTBINS(weightPATH, "BETA", function(x) x$BETA, "trait_id", include = 0)
  stopifnot(identical(unique(weightBINS$trait_id), "trait.with.dots"))
  message("Report data tests passed: ranks, signal retention, bins, workbook limits and integers.")
}
testREPORTDATA()
