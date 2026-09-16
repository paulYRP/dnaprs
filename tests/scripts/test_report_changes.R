#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))
directory <- commandArgs(trailingOnly = TRUE)[1L]
stopifnot(dir.exists(directory))
setwd(directory)
manifest <- fread("output_files.tsv")
manifest[, staged_path := file_name]
addINPUT <- function(name, value, section) {
  relative <- file.path(section, name)
  path <- file.path("inputs", relative)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  fwrite(value, path, sep = "\t")
  manifest <<- rbind(manifest, data.table(publish_path = section, file_name = name, staged_path = relative))
}
addINPUT("SYNTHETIC.genotype_eda_summary.tsv", data.table(cohort = "SYNTHETIC", role = "target",
  input_stage = "raw", format = "bed", participants = 6, variants = 12, chromosomes = 1,
  autosomal_variants = 12, x_variants = 0, status = "FAIL", review_items = 0, fail_items = 1,
  pass_items = 1, not_run_items = 1, completion = "PARTIAL"), "genotype_eda/SYNTHETIC")
addINPUT("SYNTHETIC.genotype_eda_checks.tsv", data.table(cohort = "SYNTHETIC",
  check = c("heterozygosity", "reported_sex"), status = c("FAIL", "NOT_RUN"),
  value = c("6 participants", "0 X markers"), reason = c("Invalid raw marker", "No X markers")), "genotype_eda/SYNTHETIC")
addINPUT("SYNTHETIC_qc_review.genotype_eda_checks.tsv", data.table(cohort = "SYNTHETIC_qc_review",
  check = "heterozygosity", status = "PASS", value = "6 participants", reason = "Valid corrected calculation"), "target_qc/SYNTHETIC/sample_review")
addINPUT("SYNTHETIC.target_prep_summary.tsv", data.table(cohort = "SYNTHETIC", input_stage = "raw",
  step = c("Imported", "Marker resolution", "Duplicate handling", "Allele orientation", "Prepared checkpoint"),
  participants = 6, variants = c(12, 11, 10, 10, 10), status = "PASS", stage_order = 1:5,
  removed_variants = c(0, 1, 1, 0, 0), reason = "Synthetic preparation count"), "target_prep/SYNTHETIC")
addINPUT("SYNTHETIC.marker_decisions.tsv", data.table(
  decision = c(rep("RETAINED_UNIQUE", 10), "EXCLUDED_UNRESOLVED_MARKER", "EXCLUDED_REDUNDANT_DUPLICATE_PROBE"),
  reason = c(rep("Retained", 10), "Unresolved marker", "Redundant probe"),
  final_id = c(paste0("rs", 1:10), "", "")), "target_prep/SYNTHETIC")
addINPUT("SYNTHETIC.target_qc.tsv", data.table(cohort = "SYNTHETIC", input_stage = "raw", status = "PASS",
  source_participants = 6, retained_participants = 6, source_variants = 10, retained_variants = 9, chromosomes = 1), "target_qc/SYNTHETIC")
addINPUT("SYNTHETIC.imputation_qc.tsv", data.table(cohort = "SYNTHETIC", chromosome = 1,
  corrected_typed_variants = 9, reference_matched_typed_variants = 9, unmatched_typed_variants = 0,
  retained_variants = 20, imputed_variants = 11, dr2_threshold = 0.8, status = "PASS"), "target_imputation/SYNTHETIC")
addINPUT("phenotype_associations.tsv", data.table(model_id = "synthetic", cohort = "SYNTHETIC",
  method = "plink_ct", participants = 6, beta = NA_real_, ci_low = NA_real_, ci_high = NA_real_,
  parametric_p = NA_real_, permutation_p = NA_real_, permutation_holm = NA_real_,
  permutation_scheme = "not_run", permutations = 0, delta_r2 = NA_real_, partial_r2 = NA_real_,
  fit_metric = "delta_r2", status = "NOT_RUN"), "phenotype")
addINPUT("phenotype_models_fitted.tsv", data.table(model_id = character(), cohort = character(),
  method = character(), estimator = character(), formula = character(), null_formula = character(),
  primary = logical()), "phenotype")
addINPUT("phenotype_plot_data.tsv", data.table(adjusted_outcome = numeric(), adjusted_prs = numeric(),
  observed = numeric(), fitted_full = numeric(), family = character(), residual_full = numeric()), "phenotype")
addINPUT("phenotype_permutations.tsv", data.table(status = character(), permuted_beta = numeric(),
  observed_beta = numeric()), "phenotype")
addINPUT("phenotype_influence.tsv", data.table(status = character(), beta_without = numeric(),
  full_beta = numeric()), "phenotype")
for (stage in c("corrected", "imputed")) {
  section <- paste0("checkpoints/", stage, "/SYNTHETIC", if (stage == "imputed") ".imputed" else "")
  for (extension in c("pvar", "psam", "pgen")) {
    name <- paste0("SYNTHETIC.", extension)
    path <- file.path("inputs", section, name)
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    lines <- switch(extension,
      pvar = c("##fileformat=VCFv4.2", "#CHROM\tPOS\tID\tREF\tALT",
        if (stage == "corrected") "1\t100\trs1\tA\tG" else "1\t200\t1:200:C:T\tC\tT"),
      psam = c("#FID\tIID", "S1\tS1"), pgen = "Synthetic download fixture, not an analytical PGEN")
    writeLines(lines, path)
    published <- file.path("data", section, name)
    dir.create(dirname(published), recursive = TRUE, showWarnings = FALSE)
    stopifnot(file.copy(path, published))
    manifest <- rbind(manifest, data.table(publish_path = section, file_name = name,
      staged_path = file.path(section, name)))
  }
}
fwrite(manifest, "output_files.tsv", sep = "\t")
cat("Added synthetic stage-separated QC and native checkpoint report inputs.\n")
