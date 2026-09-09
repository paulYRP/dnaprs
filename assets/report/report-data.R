# Report-only reductions. Original result files and scoring inputs are unchanged.
reportTRAITID <- function(trait) {
  # Encoding every UTF-8 byte avoids punctuation, case and sanitisation collisions.
  paste(sprintf("%02x", as.integer(charToRaw(enc2utf8(trait)))), collapse = "")
}

reportELIGIBILITY <- function(value) {
  value <- data.table::copy(value)
  for (column in c("primary_analysis", "score_eligible")) {
    if (!column %in% names(value)) value[, (column) := NA]
    flag <- toupper(trimws(as.character(value[[column]])))
    value[, (column) := match(flag, c("FALSE", "TRUE")) == 2L]
  }
  value[, analysis_set := data.table::fcase(
    primary_analysis == TRUE & score_eligible == TRUE, "Primary analysis",
    primary_analysis == FALSE & score_eligible == TRUE, "Sensitivity only",
    (is.na(primary_analysis) | primary_analysis == FALSE) & score_eligible == FALSE, "Not eligible for scoring",
    default = "Unknown eligibility"
  )]
  if (!"reason" %in% names(value)) value[, reason := "Not recorded"]
  list(totals = value[, .N, by = .(cohort, analysis_set)],
       reasons = value[, .N, by = .(cohort, analysis_set, reason)])
}

reportRANGE <- function(x) {
  x <- x[is.finite(x)]
  if (length(x)) range(x) else c(Inf, -Inf)
}

reportQQ <- function(p) {
  p <- sort(p[is.finite(p) & p > 0 & p <= 1])
  n <- length(p)
  if (!n) return(data.table::data.table(expected = numeric(), observed = numeric()))
  ranks <- sort(unique(c(
    round(seq(1, n, length.out = min(50000L, n))),
    seq_len(min(1000L, n))
  )))
  a <- if (n <= 10L) 3 / 8 else 1 / 2
  data.table::data.table(
    expected = -log10((ranks - a) / (n + 1 - 2 * a)),
    observed = -log10(p[ranks])
  )
}

reportMANHATTAN <- function(value, layout) {
  value <- value[is.finite(suppressWarnings(as.numeric(CHR))) &
    is.finite(BP) & is.finite(p) & p > 0 & p <= 1,
    .(SNP, chromosome = as.integer(CHR), BP, p)]
  data.table::setorder(value, chromosome, BP, SNP, p)
  strong <- which(value$p <= 1e-5)
  background <- which(value$p > 1e-5)
  if (length(background) > 250000L) {
    background <- background[unique(round(seq(1, length(background), length.out = 250000L)))]
  }
  full <- nrow(value)
  value <- value[sort(c(strong, background))]
  value[layout, on = "chromosome", genomic_position := BP + i.chromosome_offset]
  value[, `:=`(log10_p = -log10(p), chromosome_group = factor(chromosome %% 2L))]
  list(data = value, full = full, strong = length(strong))
}

# Reuse the pinned ggplot2 bin boundaries and edge tolerance, without allocating
# an entire ggplot layer for every input row. Regression tests check geom_histogram.
reportHISTOGRAM <- function(x, limits, bins = 50L, binwidth = NULL, boundary = NULL) {
  if (length(limits) != 2L || any(!is.finite(limits))) return(data.table::data.table())
  breaks <- if (is.null(binwidth)) {
    getFromNamespace("bin_breaks_bins", "ggplot2")(limits, bins = bins, boundary = boundary)
  } else {
    getFromNamespace("bin_breaks_width", "ggplot2")(limits, width = binwidth, boundary = boundary)
  }
  index <- getFromNamespace("bin_cut", "ggplot2")(x[is.finite(x)], breaks)
  data.table::data.table(
    xmin = head(breaks$breaks, -1L), xmax = tail(breaks$breaks, -1L),
    count = tabulate(index, nbins = length(breaks$breaks) - 1L)
  )
}

readREPORTCOLUMNS <- function(path, columns) {
  message(sprintf("Report read: %s [%s]", basename(path), paste(columns, collapse = ", ")))
  available <- names(data.table::fread(path, nrows = 0L, showProgress = FALSE))
  data.table::fread(path, select = intersect(columns, available), showProgress = FALSE)
}

prepareGWASPLOTS <- function(paths) {
  ranges <- lapply(paths, function(path) {
    value <- readREPORTCOLUMNS(path, c("CHR", "BP", "b", "freq", "p"))
    coordinates <- value[is.finite(suppressWarnings(as.numeric(CHR))) &
      is.finite(BP) & is.finite(p) & p > 0 & p <= 1,
      .(chromosome = as.integer(CHR), BP)]
    list(
      effect = reportRANGE(value$b),
      maf = reportRANGE(pmin(value$freq[value$freq > 0 & value$freq < 1],
                            1 - value$freq[value$freq > 0 & value$freq < 1])),
      layout = coordinates[, .(chromosome_length = max(BP)), by = chromosome]
    )
  })
  layout <- data.table::rbindlist(lapply(ranges, `[[`, "layout"))
  if (nrow(layout)) {
    layout <- layout[, .(chromosome_length = max(chromosome_length)), by = chromosome][order(chromosome)]
    layout[, chromosome_offset := data.table::shift(cumsum(as.numeric(chromosome_length)), fill = 0)]
    layout[, chromosome_midpoint := chromosome_offset + chromosome_length / 2]
  } else {
    layout <- data.table::data.table(chromosome = integer(), chromosome_offset = numeric(), chromosome_midpoint = numeric())
  }
  effectRange <- reportRANGE(c(0, unlist(lapply(ranges, `[[`, "effect"))))
  mafRange <- reportRANGE(unlist(lapply(ranges, `[[`, "maf")))
  records <- lapply(paths, function(path) {
    trait <- sub("\\.cojo\\.ma$", "", basename(path))
    value <- readREPORTCOLUMNS(path, c("SNP", "CHR", "BP", "b", "freq", "p"))
    qq <- reportQQ(value$p)
    manhattan <- reportMANHATTAN(value, layout)
    effect <- reportHISTOGRAM(value$b, effectRange)
    frequency <- value$freq[is.finite(value$freq) & value$freq > 0 & value$freq < 1]
    maf <- reportHISTOGRAM(pmin(frequency, 1 - frequency), mafRange)
    addTrait <- function(x) { x[, trait_id := rep(trait, .N)]; x }
    list(
      qq = addTrait(qq), manhattan = addTrait(manhattan$data),
      effect = addTrait(effect), maf = addTrait(maf),
      selection = data.table::data.table(
        trait_id = trait, figure_id = c("gwas_qq", "gwas_manhattan"),
        eligible_records = c(sum(is.finite(value$p) & value$p > 0 & value$p <= 1), manhattan$full),
        displayed_records = c(nrow(qq), nrow(manhattan$data)),
        selection_rule = c("50000 evenly spaced full-data ranks plus 1000 smallest P values",
                           "All P <= 1e-5 plus 250000 genomically spaced background records")
      )
    )
  })
  result <- lapply(c("qq", "manhattan", "effect", "maf", "selection"), function(name) {
    data.table::rbindlist(lapply(records, `[[`, name), use.names = TRUE, fill = TRUE)
  })
  names(result) <- c("qq", "manhattan", "effect", "maf", "selection")
  result$layout <- layout
  result
}

prepareREPORTBINS <- function(paths, columns, transform, group, bins = 50L,
                              binwidth = NULL, boundary = NULL, include = numeric()) {
  ranges <- lapply(paths, function(path) {
    reportRANGE(transform(readREPORTCOLUMNS(path, columns)))
  })
  limits <- reportRANGE(c(include, unlist(ranges)))
  data.table::rbindlist(lapply(paths, function(path) {
    value <- transform(readREPORTCOLUMNS(path, columns))
    result <- reportHISTOGRAM(value, limits, bins, binwidth, boundary)
    label <- if (group == "trait_id") sub("\\.sbayesrc\\.txt$", "", basename(path)) else sub("\\..*$", "", basename(path))
    result[, (group) := rep(label, .N)]
    result
  }), fill = TRUE)
}

reportWORKBOOKFITS <- function(rows, columns) rows <= 1048575L && columns <= 16384L

reportEXCELDATA <- function(value) {
  value <- data.table::copy(value)
  for (column in names(value)) {
    if (inherits(value[[column]], "integer64")) data.table::set(value, j = column, value = as.character(value[[column]]))
  }
  as.data.frame(value)
}
