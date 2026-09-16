# Write byte-preserving download assets without loading a complete file into R.
reportDOWNLOADFILE <- function(source, href, output, chunk_bytes = 4L * 1024L * 1024L) {
  stopifnot(file_test("-f", source), chunk_bytes > 0L)
  key <- as.character(openssl::md5(charToRaw(enc2utf8(href))))
  directory <- file.path("assets", "downloads", key)
  dir.create(file.path(output, directory), recursive = TRUE, showWarnings = FALSE)
  connection <- file(source, "rb")
  on.exit(close(connection))
  chunks <- 0L
  bytes <- 0
  repeat {
    value <- readBin(connection, "raw", n = chunk_bytes)
    if (!length(value)) break
    chunks <- chunks + 1L
    bytes <- bytes + length(value)
    writeLines(paste0('document.currentScript.downloadBytes="',
      openssl::base64_encode(value), '";'),
      file.path(output, directory, sprintf("%06d.js", chunks)), useBytes = TRUE)
  }
  stopifnot(bytes == file.info(source)$size)
  list(base = paste0(directory, "/"), chunks = chunks, bytes = bytes, chunkBytes = chunk_bytes)
}

# Only files linked from the rendered pages or figure manifest become download assets.
reportDOWNLOADS <- function(output) {
  links <- character()
  for (page in list.files(output, pattern = "[.]html$", full.names = TRUE)) {
    html <- paste(readLines(page, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
    anchors <- regmatches(html, gregexpr('<a\\s[^>]*>', html, perl = TRUE))[[1L]]
    anchors <- anchors[grepl('\\sdownload(?:=|\\s|>)', anchors, perl = TRUE)]
    anchors <- anchors[grepl(' href="', anchors, fixed = TRUE)]
    links <- c(links, sub('.* href="([^"]+)".*', '\\1', anchors))
  }
  figures <- data.table::fread("provenance/figure_manifest.tsv")
  for (column in intersect(c("svg", "png", "tiff", "jpeg", "source_table"), names(figures))) {
    links <- c(links, figures[[column]])
  }
  links <- unique(utils::URLdecode(gsub("&amp;", "&", links, fixed = TRUE)))
  links <- links[!is.na(links) & nzchar(links) & !grepl("^[a-zA-Z]+:|^#|^/", links)]
  manifest <- data.table::fread("provenance/output_files.tsv")
  index <- list()
  deferred <- c("provenance/execution_trace.txt", "provenance/execution_report.html",
    "provenance/execution_timeline.html", "provenance/pipeline_dag.html")
  for (href in sort(links)) {
    if (startsWith(href, "../data/checkpoints/")) {
      row <- manifest[relative_path == href]
      if (nrow(row) != 1L) stop(sprintf("Unlisted checkpoint download: %s", href))
      source <- reportINPUTPATHS(row, Sys.getenv("DNAPRS_REPORT_INPUTS"))
    } else {
      if (".." %in% strsplit(href, "/", fixed = TRUE)[[1L]]) stop(sprintf("Invalid download path: %s", href))
      source <- file.path(output, href)
    }
    if (!file.exists(source) && href %in% deferred) next
    if (!file_test("-f", source)) stop(sprintf("Missing report download: %s", href))
    index[[href]] <- reportDOWNLOADFILE(source, href, output)
  }
  writeLines(paste0("window.dnaprsDownloadIndex=", jsonlite::toJSON(index, auto_unbox = TRUE), ";"),
    file.path(output, "assets", "download-index.js"), useBytes = TRUE)
  message(sprintf("Prepared %d uncompressed offline downloads.", length(index)))
  invisible(index)
}
