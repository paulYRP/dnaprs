#!/usr/bin/env Rscript

library(data.table)
repo <- normalizePath(".")
task <- tempfile("sbayesrc-test-", tmpdir = Sys.getenv("DNAPRS_TEST_RECORDS", tempdir()))
dir.create(task)
Sys.chmod(task, "0755")
setwd(task)
prepare <- file.path(repo, "bin/prepare_sbayesrc_genotypes.R")
check <- file.path(repo, "bin/check_sbayesrc_scoring.R")
scorer <- file.path(repo, "bin/run_sbayesrc.R")
fixture <- readLines(file.path(repo, "tests/data/sbayesrc/chr1.vcf"))
keep <- file.path(repo, "tests/data/sbayesrc/keep_original.tsv")
equalKEEP <- file.path(repo, "tests/data/sbayesrc/keep.tsv")
keepCHECKSUM <- tools::md5sum(keep)
cohort <- "SYNTHETIC"
passed <- 0L
run <- function(program, args, log, expected = 0L, contains = NULL) {
  status <- system2(program, shQuote(args), stdout = log, stderr = log)
  if ((expected == 0L && status != 0L) || (expected != 0L && status == 0L)) {
    stop(paste("Unexpected exit status", status, "in", log, "\n", paste(readLines(log), collapse = "\n")))
  }
  if (!is.null(contains) && !any(grepl(contains, readLines(log), fixed = TRUE))) stop(paste("Missing expected error", contains, "in", log))
}
ok <- function(label) { passed <<- passed + 1L; message("PASS: ", label) }
writeVCF <- function(lines, path) {
  stream <- gzfile(path, "wt")
  writeLines(lines, stream)
  close(stream)
}
prepARGS <- function(vcf, keepPATH = keep) c(prepare, "--action", "prepare", "--cohort", cohort,
  "--chromosome", "1", "--vcf", vcf, "--keep", keepPATH, "--threads", "1", "--memory", "640")
writeVCF(fixture, "original.vcf.gz")
importARGS <- c("--vcf", "original.vcf.gz", "dosage=DS", "--double-id", "--keep", keep,
  "--set-missing-var-ids", "@:#:$r:$a", "--make-pgen", "--threads", "1", "--memory", "640", "--out", "original_error")
run("plink2", importARGS, "original_samples.log", 1L, "No samples remaining")
ok("reproduce zero retained samples when original family IDs are used with double-id")
importARGS[which(importARGS == "--keep") + 1L] <- equalKEEP
importARGS[which(importARGS == "--out") + 1L] <- "import_control"
run("plink2", importARGS, "import_control.log")
weights <- data.table(SNP = c(as.vector(rbind(paste0("rs", 1:22, "_A"), paste0("rs", 1:22, "_B"))), "rsABSENT"),
  A1 = c(rep(c("G", "C"), 22), "A"), BETA = c(rep(c(0.2, 0.3), 22), 0.8))
fwrite(weights, "weights.tsv", sep = "\t")
manifest <- data.table(chromosome = 1:22,
  directory = paste0(cohort, ".chr", 1:22, ".sbayesrc_genotypes"),
  qc = paste0(cohort, ".chr", 1:22, ".sbayesrc.genotype_qc.tsv"))
for (chromosome in 1:22) {
  lines <- sub("^1\t", paste0(chromosome, "\t"), fixture)
  lines <- gsub("rs1_", paste0("rs", chromosome, "_"), lines, fixed = TRUE)
  lines <- sub("ID=1>", paste0("ID=", chromosome, ">"), lines, fixed = TRUE)
  vcf <- paste0("chr", chromosome, ".vcf.gz")
  writeVCF(lines, vcf)
  args <- prepARGS(vcf)
  args[which(args == "--chromosome") + 1L] <- as.character(chromosome)
  run("Rscript", args, paste0("prepare", chromosome, ".log"))
  psam <- fread(file.path(manifest$directory[chromosome], paste0(cohort, "_chr", chromosome, ".psam")), colClasses = "character", na.strings = NULL)
  stopifnot(identical(psam[["#FID"]], c("FAM01", "FAM02")), identical(psam$IID, c("TEST01", "TEST02")))
}
stopifnot(identical(tools::md5sum(keep), keepCHECKSUM))
controlPSAM <- fread("import_control.psam", colClasses = "character", na.strings = NULL)
controlPSAM[, `#FID` := c("FAM01", "FAM02")]
preparedPSAM <- fread(file.path(manifest$directory[1L], paste0(cohort, "_chr1.psam")), colClasses = "character", na.strings = NULL)
stopifnot(identical(as.data.frame(preparedPSAM), as.data.frame(controlPSAM)),
  unname(tools::md5sum("import_control.pgen")) == unname(tools::md5sum(file.path(manifest$directory[1L], paste0(cohort, "_chr1.pgen")))))
ok("import all autosomes and restore original family IDs without changing sample order, dosages or the keep file")
fwrite(manifest[22:1], "chromosomes.tsv", sep = "\t")
assembleARGS <- c(prepare, "--action", "assemble", "--cohort", cohort, "--manifest", file.path(task, "chromosomes.tsv"), "--keep", keep)
run("Rscript", assembleARGS, "assemble.log")
directory <- file.path(task, paste0(cohort, ".sbayesrc_genotypes"))
qc <- fread(paste0(cohort, ".sbayesrc.genotype_qc.tsv"))
stopifnot(identical(qc$chromosome, 1:22), all(qc$variants == 3L), all(qc$preserved_ids == 2L), all(qc$assigned_ids == 1L), all(qc$participants == 2L))
ok("collect reversed chromosome completion order without changing participants")
checkARGS <- function(weight = file.path(task, "weights.tsv"), target = directory) c(check,
  "--cohort", cohort, "--trait-id", "TRAIT", "--weight", weight, "--target-dir", target, "--keep", keep)
run("Rscript", checkARGS(), "check.log")
matches <- fread(paste0(cohort, ".TRAIT.sbayesrc.match_qc.tsv"))
stopifnot(all(matches$usable_variants == 2L), all(matches$unmatched_weights == 1L))
ok("accept partial model coverage with compatible effect alleles")
for (multiplier in 1:2) {
  trait <- paste0("TRAIT", multiplier)
  traitWEIGHT <- file.path(task, paste0(trait, ".tsv"))
  fwrite(copy(weights)[, BETA := BETA * multiplier], traitWEIGHT, sep = "\t")
  run("Rscript", c(scorer, "--action", "score", "--input", traitWEIGHT,
    "--target-dir", directory, "--cohort", cohort, "--role", "target", "--trait-id", trait,
    "--prs-name", trait, "--keep", keep, "--plink", Sys.which("plink2")), paste0(trait, ".log"))
  score <- fread(paste0(cohort, ".", trait, ".sbayesrc.score.tsv"))
  setorder(score, IID)
  stopifnot(identical(score$FID, c("FAM01", "FAM02")), identical(score$IID, c("TEST01", "TEST02")), all(score$used_variants == 44),
    max(abs(score$raw_prs - c(4.4, 13.2) * multiplier)) < 1e-4)
}
ok("official SBayesRC scores equal known weighted dosages for two traits")

decisions <- data.table(cohort = cohort, FID = c("FAM02", "FAM01"), IID = c("TEST02", "TEST01"),
  technical_pass = TRUE, score_eligible = TRUE, related_flag = c(FALSE, TRUE),
  sample_missingness_pass = TRUE, heterozygosity_z = c(0.2, 0.1), heterozygosity_pass = TRUE,
  sex_check_pass = TRUE, ancestry_flag = "PASS", ancestry_distance = c(2, 1), primary_analysis = c(TRUE, FALSE))
fwrite(decisions, "decisions.tsv", sep = "\t")
run("Rscript", c(file.path(repo, "bin/combine_scores.R"),
  "--scores", paste(paste0(cohort, ".TRAIT", 1:2, ".sbayesrc.score.tsv"), collapse = ","),
  "--participant-decisions", "decisions.tsv"), "combine.log")
combined <- fread("prs_scores_long.tsv")
stopifnot(nrow(combined) == 4L, identical(combined$FID, rep(c("FAM01", "FAM02"), 2)),
  identical(combined$IID, rep(c("TEST01", "TEST02"), 2)),
  identical(combined$primary_analysis, rep(c(FALSE, TRUE), 2)),
  identical(combined$ancestry_distance, rep(1:2, 2)),
  max(abs(combined$raw_prs - rep(c(4.4, 13.2), 2) * rep(1:2, each = 2))) < 1e-4,
  max(abs(combined$prs_z - rep(c(-1, 1) / sqrt(2), 2))) < 1e-10,
  nrow(fread("prs_scores_wide.tsv")) == 2L, all(fread("score_qc.tsv")$status == "PASS"))
ok("join original participant decisions and standardise scores after real SBayesRC scoring")
dir.create("mismatch")
for (chromosome in 2:22) {
  stopifnot(all(file.copy(file.path(directory, paste0(cohort, "_chr", chromosome, c(".pgen", ".pvar", ".psam"))), "mismatch")))
}
run("plink2", c("--pfile", file.path(directory, paste0(cohort, "_chr1")), "--set-all-var-ids", "@:#:$r:$a",
  "--make-pgen", "--memory", "640", "--threads", "1", "--out", file.path("mismatch", paste0(cohort, "_chr1"))), "rename.log")
run("Rscript", checkARGS(target = file.path(task, "mismatch")), "mismatch.log", 1L, "Chromosome 1 has no usable variants: 0 matched IDs")
fwrite(weights[, .(SNP)], "extract.txt", col.names = FALSE)
run("plink2", c("--pfile", file.path("mismatch", paste0(cohort, "_chr1")), "--extract", "extract.txt",
  "--score", "weights.tsv", "1", "2", "3", "header", "no-mean-imputation", "--memory", "640", "--threads", "1", "--out", "original_error"),
  "original_error.log", 1L, "No variants remaining")
ok("reproduce the original zero-overlap failure and detect it before scoring")
bad <- copy(weights)[, A1 := "N"]
fwrite(bad, "bad_alleles.tsv", sep = "\t")
run("Rscript", checkARGS(file.path(task, "bad_alleles.tsv")), "bad_alleles.log", 1L, "2 incompatible alleles")
bad <- copy(weights)
bad[1L, A1 := "T"]
fwrite(bad, "partial_alleles.tsv", sep = "\t")
run("Rscript", checkARGS(file.path(task, "partial_alleles.tsv")), "partial_alleles.log")
stopifnot(fread(paste0(cohort, ".TRAIT.sbayesrc.match_qc.tsv"))[1L, status] == "REVIEW")
ok("report incompatible alleles without silently changing weights")
caseDIR <- function(name) { dir.create(file.path(task, name)); setwd(file.path(task, name)) }
caseDIR("equal_ids")
run("Rscript", prepARGS(file.path(task, "chr1.vcf.gz"), equalKEEP), "test.log")
equalPSAM <- fread(file.path(manifest$directory[1L], paste0(cohort, "_chr1.psam")), colClasses = "character")
stopifnot(identical(equalPSAM[["#FID"]], c("TEST01", "TEST02")), identical(equalPSAM$IID, c("TEST01", "TEST02")))
ok("retain support for equal family and individual IDs")
caseDIR("character_ids")
writeVCF(gsub("TEST02", "SUB_02", gsub("TEST01", "0007", fixture, fixed = TRUE), fixed = TRUE), "character.vcf.gz")
writeLines(c("FAM_02\tSUB_02", "0001\t0007"), "keep.tsv")
run("Rscript", prepARGS("character.vcf.gz", "keep.tsv"), "test.log")
characterPSAM <- fread(file.path(manifest$directory[1L], paste0(cohort, "_chr1.psam")), colClasses = "character")
stopifnot(identical(characterPSAM[["#FID"]], c("0001", "FAM_02")), identical(characterPSAM$IID, c("0007", "SUB_02")))
ok("preserve leading zeros and underscores in participant identifiers")
caseDIR("ambiguous_ids")
writeLines(c("FAM01\tTEST01", "FAM02\tTEST01"), "keep.tsv")
run("Rscript", prepARGS(file.path(task, "chr1.vcf.gz"), "keep.tsv"), "test.log", 1L, "ambiguous IID-to-FID mappings")
ok("reject an individual ID mapped to more than one family")
caseDIR("duplicate_samples")
writeVCF(gsub("TEST02", "TEST01", fixture, fixed = TRUE), "duplicate_samples.vcf.gz")
run("Rscript", prepARGS("duplicate_samples.vcf.gz"), "test.log", 1L, "duplicate sample IDs")
ok("reject duplicate VCF sample names before import")
caseDIR("empty_sample")
writeVCF(sub("\tTEST03$", "\t", fixture), "empty_sample.vcf.gz")
run("Rscript", prepARGS("empty_sample.vcf.gz"), "test.log", 1L, "non-empty sample IDs")
ok("reject an empty VCF sample name before import")
caseDIR("duplicates")
writeVCF(sub("rs1_B", "rs1_A", fixture, fixed = TRUE), "duplicate.vcf.gz")
run("Rscript", prepARGS("duplicate.vcf.gz"), "test.log", 1L, "duplicate variant IDs")
ok("reject duplicate variant IDs without removing records")
caseDIR("missing_participant")
writeLines(c(readLines(keep), "MISSING\tMISSING"), "keep.tsv")
run("Rscript", prepARGS(file.path(task, "chr1.vcf.gz"), "keep.tsv"), "test.log", 1L, "missing eligible IIDs")
ok("reject a missing eligible participant")
caseDIR("missing_chromosome")
manifest[, directory := file.path(task, directory)]
manifest[, qc := file.path(task, qc)]
fwrite(manifest[1:21], "chromosomes.tsv", sep = "\t")
args <- assembleARGS
args[which(args == "--manifest") + 1L] <- "chromosomes.tsv"
run("Rscript", args, "test.log", 1L, "exactly one prepared bundle")
ok("reject incomplete autosomal collection")
caseDIR("sample_order")
dir.create("chr2")
stopifnot(all(file.copy(list.files(manifest$directory[2], full.names = TRUE), "chr2")))
psamPATH <- file.path("chr2", paste0(cohort, "_chr2.psam"))
psam <- fread(psamPATH)
fwrite(psam[2:1], psamPATH, sep = "\t")
manifest[2L, directory := normalizePath("chr2")]
fwrite(manifest, "chromosomes.tsv", sep = "\t")
run("Rscript", args, "test.log", 1L, "different participant order")
ok("reject mixed participant order during collection")
message(sprintf("Passed %s SBayesRC genotype and scoring regression checks. Records: %s", passed, task))
