#!/usr/bin/env Rscript

argument <- commandArgs(trailingOnly = TRUE)
name <- sub("^--", "", argument[seq.int(1L, length(argument), 2L)])
value <- argument[seq.int(2L, length(argument), 2L)]
option <- stats::setNames(as.list(value), name)
if (!requireNamespace("data.table", quietly = TRUE)) stop("The data.table package is required.", call. = FALSE)

readSCORE <- function(path) {
  value <- data.table::fread(path, colClasses = "character")
  iid <- intersect(c("IID", "#IID"), names(value))
  fid <- intersect(c("#FID", "FID"), names(value))
  pc <- grep("^PC[0-9]+_AVG$", names(value), value = TRUE)
  if (length(iid) != 1L || length(pc) < 2L) stop(sprintf("Invalid PCA projection: %s", path), call. = FALSE)
  result <- value[, c(fid[1L], iid[1L], pc), with = FALSE]
  data.table::setnames(result, c("FID", "IID", sub("_AVG$", "", pc)))
  for (column in sub("_AVG$", "", pc)) result[[column]] <- suppressWarnings(as.numeric(result[[column]]))
  result
}

reference <- readSCORE(option[["reference"]])
target <- readSCORE(option[["target"]])
population <- data.table::fread(option[["population"]], colClasses = "character")
sample_column <- names(population)[tolower(names(population)) == "sample"]
super_column <- names(population)[tolower(names(population)) == "super_pop"]
pop_column <- names(population)[tolower(names(population)) == "pop"]
if (length(sample_column) != 1L || length(super_column) != 1L) {
  stop("Population metadata requires sample and super_pop columns.", call. = FALSE)
}
population <- population[, .(
  IID = get(sample_column),
  population = if (length(pop_column) == 1L) get(pop_column) else "",
  super_population = get(super_column)
)]
reference <- merge(reference, population, by = "IID", all.x = TRUE, sort = FALSE)
if (anyNA(reference$super_population)) stop("Reference participants are missing population metadata.", call. = FALSE)

pc <- intersect(grep("^PC[0-9]+$", names(reference), value = TRUE), names(target))
european <- reference[super_population == "EUR"]
if (nrow(european) < 20L) stop("At least 20 unrelated European reference participants are required for ancestry distance.", call. = FALSE)
pc <- pc[seq_len(min(length(pc), nrow(european) - 2L))]
if (length(pc) < 2L) stop("At least two usable reference PCs are required.", call. = FALSE)

eur_matrix <- as.matrix(european[, ..pc])
target_matrix <- as.matrix(target[, ..pc])
if (any(!is.finite(eur_matrix)) || any(!is.finite(target_matrix))) stop("PCA projections contain non-finite values.", call. = FALSE)
centre <- colMeans(eur_matrix)
covariance <- stats::cov(eur_matrix)
ridge <- max(mean(diag(covariance)), .Machine$double.eps) * 1e-8
inverse <- solve(covariance + diag(ridge, nrow(covariance)))
distanceVALUE <- function(matrix) {
  centred <- sweep(matrix, 2L, centre, "-")
  sqrt(rowSums((centred %*% inverse) * centred))
}
eur_distance <- distanceVALUE(eur_matrix)
percentile <- as.numeric(option[["percentile"]])
if (!is.finite(percentile) || percentile <= 0 || percentile >= 1) stop("Ancestry percentile must be between 0 and 1.", call. = FALSE)
threshold <- as.numeric(stats::quantile(eur_distance, probs = percentile, names = FALSE, type = 8))
target_distance <- distanceVALUE(target_matrix)

reference[, ancestry_distance := distanceVALUE(as.matrix(reference[, ..pc]))]
reference[, ancestry_flag := ifelse(super_population == "EUR" & ancestry_distance <= threshold, "EUR_REFERENCE", "REFERENCE")]
data.table::setorder(reference, super_population, population, IID)
data.table::fwrite(reference, option[["reference-output"]], sep = "\t", na = "NA")

result <- data.table::copy(target)
result[, `:=`(
  cohort = option[["cohort"]],
  ancestry_distance = target_distance,
  ancestry_threshold = threshold,
  ancestry_percentile = percentile,
  ancestry_flag = ifelse(target_distance <= threshold, "PASS", "OUTLIER")
)]
data.table::setcolorder(result, c("cohort", "FID", "IID", pc, "ancestry_distance", "ancestry_threshold", "ancestry_percentile", "ancestry_flag"))
data.table::setorder(result, FID, IID)
data.table::fwrite(result, option[["target-output"]], sep = "\t", na = "NA")

summary <- data.table::data.table(
  cohort = option[["cohort"]],
  matched_variants = as.integer(option[["matched-variants"]]),
  pruned_variants = as.integer(option[["pruned-variants"]]),
  reference_participants = nrow(reference),
  european_reference_participants = nrow(european),
  target_participants = nrow(target),
  pcs = length(pc),
  percentile = percentile,
  distance_threshold = threshold,
  ancestry_pass = sum(result$ancestry_flag == "PASS"),
  ancestry_outlier = sum(result$ancestry_flag != "PASS"),
  status = "PASS"
)
data.table::fwrite(summary, option[["summary-output"]], sep = "\t", na = "NA")
