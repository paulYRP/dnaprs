# Quarto discovers resources before rendering; table chunks are created during knitting.
output <- Sys.getenv("QUARTO_PROJECT_OUTPUT_DIR", "_site")
files <- list.files("assets/tables", recursive = TRUE, full.names = TRUE)
for (path in files) {
  destination <- file.path(output, path)
  dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
  if (!file.copy(path, destination, overwrite = TRUE)) {
    stop(sprintf("Cannot publish report table asset: %s", path))
  }
}
