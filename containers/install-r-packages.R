args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
    stop("Usage: install-r-packages.R <repository> <manifest.tsv> [manifest.tsv ...]")
}

repository <- args[[1]]
manifest_files <- args[-1]
manifest <- do.call(
    rbind,
    lapply(manifest_files, read.delim, check.names = FALSE, stringsAsFactors = FALSE)
)
manifest <- manifest[!duplicated(manifest$package, fromLast = TRUE), ]

options(timeout = max(600, getOption("timeout")), repos = c(CRAN = repository))

for (row in seq_len(nrow(manifest))) {
    package <- manifest$package[[row]]
    install_version <- manifest$install_version[[row]]
    expected_version <- manifest$expected_version[[row]]

    for (attempt in seq_len(3)) {
        try(
            remotes::install_version(
                package,
                version = install_version,
                repos = repository,
                upgrade = "never"
            ),
            silent = FALSE
        )

        installed <- requireNamespace(package, quietly = TRUE)
        correct_version <- installed &&
            identical(as.character(packageVersion(package)), expected_version)

        if (correct_version) {
            break
        }

        if (attempt == 3) {
            stop(
                sprintf(
                    "Failed to install %s %s after %d attempts",
                    package,
                    expected_version,
                    attempt
                )
            )
        }

        unlink(Sys.glob(file.path(.libPaths()[[1]], "00LOCK*")), recursive = TRUE)
        Sys.sleep(5 * attempt)
    }
}

observed <- vapply(
    manifest$package,
    function(package) as.character(packageVersion(package)),
    character(1)
)

if (any(observed != manifest$expected_version)) {
    stop(
        paste(
            "Package version mismatch:",
            paste(manifest$package[observed != manifest$expected_version], collapse = ", ")
        )
    )
}
