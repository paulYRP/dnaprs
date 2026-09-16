# Check report sources before copying downloads or calculating file checksums.
validateREPORTINPUTS <- function(inputROOT, manifest) {
  if (!dir.exists(inputROOT)) {
    stop(sprintf("Report input directory does not exist: %s", inputROOT), call. = FALSE)
  }
  paths <- reportINPUTPATHS(manifest, inputROOT)
  for (path in paths) {
    if (dir.exists(path)) {
      stop(
        sprintf("Report input is a directory: %s. Supply individual report files and publish genotype directories separately.", path),
        call. = FALSE
      )
    }
    if (!file_test("-f", path)) {
      stop(
        sprintf("Report input is missing or is not a regular file: %s. Check the source file and any symbolic-link target.", path),
        call. = FALSE
      )
    }
    if (file.access(path, mode = 4L) != 0L) {
      stop(sprintf("Report input is not readable: %s. Check file permissions.", path), call. = FALSE)
    }
  }
  paths
}

# The relative publication path distinguishes identical filenames at different stages.
reportINPUTPATHS <- function(files, root = inputROOT) {
  staged <- if ("staged_path" %in% names(files)) files$staged_path else files$file_name
  file.path(root, staged)
}

reportRESULTPATH <- function(name, section = NULL) {
  rows <- manifest[manifest$file_name == name]
  if (!is.null(section)) rows <- rows[grepl(section, publish_path)]
  if (nrow(rows) > 1L) stop(sprintf("Report result '%s' is ambiguous; select its stage path.", name), call. = FALSE)
  if (!nrow(rows)) return("")
  reportINPUTPATHS(rows)
}

reportREADRESULTS <- function(pattern, section = NULL) {
  files <- manifest[grepl(pattern, file_name)]
  if (!is.null(section)) files <- files[grepl(section, publish_path)]
  data.table::rbindlist(lapply(seq_len(nrow(files)), function(i) {
    value <- data.table::fread(reportINPUTPATHS(files[i]), showProgress = FALSE)
    if (!nrow(value)) return(NULL)
    if (!"cohort" %in% names(value)) value[, cohort := sub("\\..*$", "", files$file_name[i])]
    value[, `:=`(source_file = files$file_name[i], source_path = files$publish_path[i])]
    if (grepl("^(genotype_eda|target_qc)/", files$publish_path[i])) {
      value[, `:=`(
        analysis_stage = if (grepl("^genotype_eda/", files$publish_path[i])) "Input EDA" else "Corrected Target QC",
        source_cohort = sub("^[^/]+/([^/]+).*", "\\1", files$publish_path[i])
      )]
    }
    value
  }), use.names = TRUE, fill = TRUE)
}

# Preserve older status labels while distinguishing failed attempts from skipped checks.
reportEDACHECKS <- function(checks) {
  checks <- data.table::copy(checks)
  if (!nrow(checks)) return(checks)
  failed <- checks$status == "NOT_RUN" & (
    (checks$check == "reported_sex" & checks$reason == "PLINK could not complete the sex check; inspect the stage log.") |
    (checks$check == "internal_pca" & checks$reason == "PLINK could not calculate at least two internal PCs; inspect marker count and the stage log.")
  )
  failed[is.na(failed)] <- FALSE
  if (any(failed)) {
    checks[, recorded_status := status]
    checks[failed, status := "FAIL"]
  }
  checks
}

# Derive coverage from the recorded checks, including reports from earlier runs.
reportEDASUMMARY <- function(summary, checks) {
  summary <- data.table::copy(summary)
  if (!nrow(summary) || !nrow(checks)) return(summary)
  original <- summary$status
  for (i in seq_len(nrow(summary))) {
    rows <- checks[checks$source_path == summary$source_path[i] & checks$cohort == summary$cohort[i]]
    if (!nrow(rows)) next
    count <- table(factor(rows$status, levels = c("PASS", "REVIEW", "FAIL", "NOT_RUN")))
    summary[i, `:=`(pass_items = as.integer(count[1]), review_items = as.integer(count[2]),
      fail_items = as.integer(count[3]), not_run_items = as.integer(count[4]),
      completion = if (count[4] > 0L) "PARTIAL" else "COMPLETE",
      status = if (count[3] > 0L) "FAIL" else if (count[2] > 0L) "REVIEW" else "PASS")]
  }
  if (any(original != summary$status, na.rm = TRUE)) summary[, recorded_status := original]
  summary
}

# Recover historical preparation counts only when the marker audit and final PVAR agree.
reportPREPARATION <- function(summary) {
  summary <- data.table::copy(summary)
  removals <- list()
  if (!nrow(summary)) return(list(summary = summary, removals = data.table::data.table()))
  summary[, count_source := "Recorded preparation summary"]
  for (cohort_id in unique(summary$cohort)) {
    audit <- manifest[file_name == paste0(cohort_id, ".marker_decisions.tsv") &
      publish_path == paste0("target_prep/", cohort_id)]
    if (nrow(audit) != 1L) next
    markers <- data.table::fread(reportINPUTPATHS(audit),
      select = c("decision", "reason", "final_id"), showProgress = FALSE)
    removed <- markers[grepl("^EXCLUDED_", decision), .(variants = .N), by = .(decision, reason)]
    removed[, cohort := cohort_id]
    removals[[cohort_id]] <- removed
    row <- summary[cohort == cohort_id]
    if (nrow(row) != 1L || row$input_stage != "raw" ||
        ("stage_order" %in% names(row) && !is.na(row$stage_order))) next
    retained <- markers[grepl("^RETAINED_", decision)]
    resolved <- sum(grepl("^RETAINED_|^EXCLUDED_.*DUPLICATE", markers$decision))
    checkpoint <- manifest[publish_path == paste0("checkpoints/corrected/", cohort_id) &
      file_name == paste0(cohort_id, ".pvar")]
    if (nrow(checkpoint) != 1L) next
    final_ids <- data.table::fread(reportINPUTPATHS(checkpoint), skip = "#CHROM", select = "ID", showProgress = FALSE)$ID
    raw <- reportREADRESULTS("[.]genotype_eda_summary[.]tsv$", paste0("^genotype_eda/", cohort_id, "$"))
    valid <- all(grepl("^RETAINED_|^EXCLUDED_UNRESOLVED_MARKER$|^EXCLUDED_.*DUPLICATE", markers$decision)) &&
      nrow(retained) == row$variants && data.table::uniqueN(retained$final_id) == row$variants &&
      length(final_ids) == row$variants && setequal(final_ids, retained$final_id) &&
      (nrow(raw) == 0L || (nrow(raw) == 1L && raw$variants == nrow(markers)))
    if (!isTRUE(valid)) {
      summary[cohort == cohort_id, count_source := "Recorded summary; intermediate counts could not be reconciled"]
      next
    }
    stages <- row[rep(1L, 5L)]
    stages[, `:=`(
      step = c("Imported", "Marker resolution", "Duplicate handling", "Allele orientation", "Prepared checkpoint"),
      stage_order = 1:5,
      variants = c(nrow(markers), resolved, row$variants, row$variants, row$variants),
      removed_variants = c(0L, nrow(markers) - resolved, resolved - row$variants, 0L, 0L),
      count_source = "Marker decisions reconciled with the final PVAR"
    )]
    summary <- data.table::rbindlist(list(summary[cohort != cohort_id], stages), use.names = TRUE, fill = TRUE)
  }
  list(summary = summary, removals = data.table::rbindlist(removals, use.names = TRUE, fill = TRUE))
}
