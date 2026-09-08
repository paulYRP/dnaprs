# Check report sources before copying downloads or calculating file checksums.
validateREPORTINPUTS <- function(inputROOT, manifest) {
  if (!dir.exists(inputROOT)) {
    stop(sprintf("Report input directory does not exist: %s", inputROOT), call. = FALSE)
  }
  paths <- file.path(inputROOT, manifest$file_name)
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
