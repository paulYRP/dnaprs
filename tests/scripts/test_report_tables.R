suppressPackageStartupMessages(library(data.table))
source("assets/report/report-tables.R")
viewerPath <- normalizePath("assets/report/dnaprs-tables.js", winslash = "/")
tablePath <- normalizePath("assets/report/report-tables.R", winslash = "/")
postPath <- normalizePath("assets/report/report-post-render.R", winslash = "/")
arguments <- commandArgs(trailingOnly = TRUE)
destination <- if (length(arguments)) arguments[1L] else tempfile("dnaprs-report-tables-")
stopifnot(!dir.exists(destination))
dir.create(destination, recursive = TRUE)
setwd(destination)
stopifnot(file.copy(tablePath, "report-tables.R"), file.copy(postPath, "report-post-render.R"))
options(dnaprs.report.chunk_rows = 17L)
htmlESCAPE <- function(x) {
  x <- gsub("&", "&amp;", as.character(x), fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  gsub('"', "&quot;", x, fixed = TRUE)
}
downloadBUTTON <- function(path, label) sprintf('<a href="%s" download>%s</a>', htmlESCAPE(path), htmlESCAPE(label))
tableN <- 0L
inputROOT <- "inputs"
dir.create(inputROOT, showWarnings = FALSE)
value <- data.table(SNP = rep(c("rs10", "rs2", "rs10"), 29),
  value = c(NA, seq_len(86)), text = rep(c("quoted\"value", "two\nlines", "</script><b>literal</b>"), 29))
for (i in 1:30) value[, (paste0("field", i)) := paste0("value", i)]
fwrite(value, "inputs/test.tsv", sep = "\t", na = "NA")
before <- tools::md5sum("inputs/test.tsv")
manifest <- data.table(file_name = "test.tsv", relative_path = "inputs/test.tsv")
html <- as.character(reportPAGEDTABLE(reportTABLESOURCE("test.tsv"), "Synthetic complete table"))
stopifnot(identical(before, tools::md5sum("inputs/test.tsv")), !grepl("<td", html, fixed = TRUE))
stopifnot(!grepl("Browse complete", html, fixed = TRUE), !grepl("<details", html, fixed = TRUE),
  grepl("Download test.tsv", html, fixed = TRUE), !grepl("Download complete", html, fixed = TRUE))
chunks <- list.files("assets/tables/-1", pattern = "[.]js$", full.names = TRUE)
stopifnot(length(chunks) == 6L)
values <- lapply(chunks, function(path) {
  text <- paste(readLines(path, warn = FALSE), collapse = "\n")
  text <- sub('^window[.]dnaprsTableChunk\\("[^"]+",', "", text)
  jsonlite::fromJSON(sub("\\);$", "", text))
})
restored <- do.call(rbind, values)
expected <- copy(value)
expected[] <- lapply(expected, function(column) { column <- as.character(column); column[is.na(column)] <- "NA"; column })
stopifnot(nrow(restored) == nrow(value), ncol(restored) == ncol(value) + 2L)
if (!identical(unname(restored[, seq_len(ncol(value))]), unname(as.matrix(expected)))) {
  print(restored[1:3, 1:3]); print(expected[1:3, 1:3]);
  stop("Chunk text differs from the source")
}
file.copy(viewerPath, "viewer.js", overwrite = TRUE)
single <- data.frame(ID = c("0001", "0002"), stringsAsFactors = FALSE)
singleHTML <- as.character(reportPAGEDTABLE(single, "Single-column table"))
stopifnot(file.exists("assets/tables/-2/complete.tsv"))
singleChunk <- readLines("assets/tables/-2/chunk-00001.js", warn = FALSE)
stopifnot(grepl('[["0001"],["0002"]]', singleChunk, fixed = TRUE))
dump(c("htmlESCAPE", "downloadBUTTON"), file = "table-helpers.R")
writeLines(c("---", 'title: "Report table test"', "format: html", "---", "", html,
             "```{r, echo=FALSE, results='asis'}", "library(data.table)",
             'source("report-tables.R")', 'source("table-helpers.R")', "tableN <- 0L",
             'reportPAGEDTABLE(data.frame(ID = sprintf("%04d", 1:220)), "Created during rendering")',
             "```", '<script src="viewer.js"></script>'), "index.qmd")
writeLines(c("project:", "  type: website", "  post-render: report-post-render.R", "  resources:", "    - assets/", "    - inputs/", "    - viewer.js",
             "website:", "  search: false"), "_quarto.yml")
message("Small table test passed: six chunks, repeated IDs, quoted/multiline text, wide rows and unchanged download.")

dir.create("inputs/corrected")
dir.create("inputs/imputed")
writeLines(c("##fileformat=VCFv4.2", "##contig=<ID=1,length=249250621>",
  "#CHROM\tPOS\tID\tREF\tALT", "1\t100\trs1\tA\tG"), "inputs/corrected/TEST.pvar")
writeLines(c("##fileformat=VCFv4.2", "#CHROM\tPOS\tID\tREF\tALT",
  "1\t200\t1:200:C:T\tC\tT"), "inputs/imputed/TEST.pvar")
csvVALUE <- data.table(UID = c("0001", "0002"), note = c('comma, and "quote"', "two\nlines"))
fwrite(csvVALUE, "inputs/phenoPRS.csv", sep = ",")
manifest <- data.table(file_name = c("TEST.pvar", "TEST.pvar", "phenoPRS.csv"),
  publish_path = c("checkpoints/corrected/TEST", "checkpoints/imputed/TEST.imputed", "phenotype"),
  staged_path = c("corrected/TEST.pvar", "imputed/TEST.pvar", "phenoPRS.csv"),
  relative_path = c("data/checkpoints/corrected/TEST/TEST.pvar", "data/checkpoints/imputed/TEST.imputed/TEST.pvar", "inputs/phenoPRS.csv"))
for (stage in c("corrected", "imputed")) {
  html <- as.character(reportPAGEDTABLE(reportTABLESOURCE("[.]pvar$", paste0("/", stage, "/")), stage))
  stopifnot(grepl(paste0("data/checkpoints/", stage, "/"), html, fixed = TRUE),
    grepl('"#CHROM"', html, fixed = TRUE), !grepl('"##fileformat', html, fixed = TRUE))
}
invisible(reportPAGEDTABLE(reportTABLESOURCE("^phenoPRS[.]csv$"), "Phenotype"))
text <- paste(readLines(sprintf("assets/tables/-%s/chunk-00001.js", tableN), warn = FALSE), collapse = "\n")
parsed <- jsonlite::fromJSON(sub("\\);$", "", sub('^window[.]dnaprsTableChunk\\("[^"]+",', "", text)))
stopifnot(identical(unname(parsed[, 1:2]), unname(as.matrix(csvVALUE))))
cat("Stage-specific PVAR headers and quoted, multiline CSV downloads passed.\n")
