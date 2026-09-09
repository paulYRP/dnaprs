arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments)) setwd(arguments[1L])
count <- 0L
for (page in list.files("_site", pattern = "[.]html$", full.names = TRUE)) {
  html <- paste(readLines(page, warn = FALSE), collapse = "\n")
  blocks <- regmatches(html, gregexpr('data-table-metadata[^>]*>.*?</script>', html, perl = TRUE))[[1L]]
  for (block in blocks) {
    metadata <- jsonlite::fromJSON(sub('</script>$', '', sub('^[^>]*>', '', block)))
    for (chunk in metadata$chunks$file) {
      source <- file.path(metadata$base, chunk)
      published <- file.path("_site", source)
      stopifnot(file.exists(published),
                identical(unname(tools::md5sum(source)), unname(tools::md5sum(published))))
      count <- count + 1L
    }
  }
}
stopifnot(count > 0L)
message(sprintf("Published report passed: %d table chunks present and unchanged.", count))
