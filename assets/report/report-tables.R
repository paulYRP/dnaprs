# Detailed tables are external assets, never complete HTML tables.
reportTABLESOURCE <- function(pattern) {
  files <- manifest[grepl(pattern, file_name)]
  structure(list(files = files), class = "report_table_source")
}

reportJSON <- function(value) {
  text <- jsonlite::toJSON(value, auto_unbox = TRUE, dataframe = "values", na = "null", null = "null", digits = NA)
  gsub("<", "\\u003c", as.character(text), fixed = TRUE)
}

reportPAGEDTABLE <- function(value, caption, page = 25L) {
  page <- if (page %in% c(25L, 50L, 100L)) page else 25L
  tableN <<- tableN + 1L
  document <- gsub("[^A-Za-z0-9_-]", "_", knitr::current_input())
  key <- paste0(document, "-", tableN)
  directory <- file.path("assets", "tables", key)
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  source <- inherits(value, "report_table_source")
  files <- if (source) value$files else NULL
  columns <- if (source) unique(unlist(lapply(files$file_name, function(name) {
    c(names(data.table::fread(file.path(inputROOT, name), nrows = 0L)), "cohort", "source_file")
  }))) else names(value)
  if (!length(columns)) return(knitr::asis_output("<p>No records were produced for this section.</p>"))
  chunkSize <- max(1L, min(as.integer(getOption("dnaprs.report.chunk_rows", 5000L)),
                            floor(100000L / length(columns))))
  count <- 0L
  descriptors <- list()
  identifier <- intersect(c("SNP", "ID", "IID", "UID", "Participant"), columns)
  identifier <- if (length(identifier)) identifier[1L] else ""
  writeChunk <- function(rows) {
    number <- length(descriptors) + 1L
    rows <- as.data.frame(rows)
    for (column in setdiff(columns, names(rows))) rows[[column]] <- NA_character_
    rows <- rows[, columns, drop = FALSE]
    # Preserve source text and avoid rounding large identifiers in JavaScript.
    rows[] <- lapply(rows, as.character)
    values <- rows
    file <- sprintf("chunk-%05d.js", number)
    writeLines(paste0("window.dnaprsTableChunk(", reportJSON(paste0(key, ":", number)),
                     ",", reportJSON(values), ");"), file.path(directory, file), useBytes = TRUE)
    ids <- if (nzchar(identifier)) rows[[identifier]] else character()
    ids <- ids[!is.na(ids)]
    # ASCII bounds use the same ordering in R and JavaScript; other IDs are scanned.
    bounds <- if (length(ids) && all(grepl("^[ -~]*$", ids))) sort(ids, method = "radix") else character()
    descriptors[[number]] <<- list(file = file, start = count, rows = nrow(rows),
      idMin = if (length(bounds)) bounds[1L] else NULL,
      idMax = if (length(bounds)) tail(bounds, 1L) else NULL)
    count <<- count + nrow(rows)
  }
  if (source) {
    for (i in seq_len(nrow(files))) {
      path <- file.path(inputROOT, files$file_name[i])
      connection <- file(path, open = "rt", encoding = "UTF-8")
      tryCatch({
        header <- names(data.table::fread(path, nrows = 0L))
        readLines(connection, n = 1L, warn = FALSE)
        repeat {
          first <- readLines(connection, n = 1L, warn = FALSE)
          if (!length(first)) break
          pushBack(first, connection)
          rows <- utils::read.table(connection, header = FALSE, sep = "\t", quote = '"',
            nrows = chunkSize, col.names = header, colClasses = "character", row.names = NULL,
            comment.char = "", check.names = FALSE, na.strings = NULL, fill = FALSE,
            blank.lines.skip = TRUE)
          if (!nrow(rows)) break
          if (!"cohort" %in% names(rows)) rows$cohort <- sub("\\..*$", "", files$file_name[i])
          rows$source_file <- files$file_name[i]
          writeChunk(rows)
          if (nrow(rows) < chunkSize) break
        }
      }, finally = close(connection))
    }
  } else {
    if (nrow(value)) for (start in seq.int(1L, nrow(value), by = chunkSize)) {
      writeChunk(value[start:min(nrow(value), start + chunkSize - 1L), , drop = FALSE])
    }
  }
  downloads <- if (source) paste(vapply(seq_len(nrow(files)), function(i) {
    downloadBUTTON(files$relative_path[i], paste("Download complete", files$file_name[i]))
  }, character(1)), collapse = " ") else {
    path <- file.path(directory, "complete.tsv")
    data.table::fwrite(value, path, sep = "\t", na = "NA")
    downloadBUTTON(gsub("\\\\", "/", path), "Download complete table (TSV)")
  }
  metadata <- list(key = key, columns = unname(as.list(columns)), total = count, identifier = identifier,
                   chunks = descriptors, base = paste0(gsub("\\\\", "/", directory), "/"))
  knitr::asis_output(paste0(
    '<section class="dnaprs-paged-table" data-paged-table data-page-size="', page, '">',
    '<p>', htmlESCAPE(caption), '</p><p>', downloads, '</p>',
    '<script type="application/json" data-table-metadata>', reportJSON(metadata), '</script>',
    '<div class="dnaprs-table-controls">',
    '<label>Rows <select data-rows><option>25</option><option>50</option><option>100</option></select></label>',
    '<button type="button" data-first>First</button><button type="button" data-previous>Previous</button>',
    '<label>Page <input data-page type="number" min="1" value="1" style="width:6em"></label>',
    '<button type="button" data-next>Next</button><button type="button" data-last>Last</button>',
    '<label>Column <select data-column></select></label>',
    '<label>Condition <select data-operator><option value="contains">contains</option>',
    '<option value="equals">equals</option><option value="lt">less than</option>',
    '<option value="lte">less than or equal</option><option value="gt">greater than</option>',
    '<option value="gte">greater than or equal</option></select></label>',
    '<label>Value <input data-query type="search"></label>',
    '<button type="button" data-apply>Apply</button><button type="button" data-find>Find identifier</button>',
    '<button type="button" data-clear>Clear</button>',
    '<button type="button" data-cancel disabled>Cancel</button></div>',
    '<p data-status role="status" aria-live="polite">', count, ' records. Open the table to browse.</p>',
    '<details data-open><summary>Browse complete table</summary>',
    '<div style="overflow:auto;max-height:65vh"><table class="table table-striped table-sm">',
    '<thead></thead><tbody></tbody></table></div></details>',
    '<noscript>JavaScript is required for browsing. Use the complete-file downloads above.</noscript></section>'
  ))
}
