ARG ANALYSIS_IMAGE=ghcr.io/paulyrp/dnaprs-analysis:1.0.0
FROM ${ANALYSIS_IMAGE}

ARG QUARTO_VERSION=1.7.32
ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        curl \
        libfontconfig1-dev \
        libtiff-dev \
        libuv1-dev \
    && curl --fail --location --retry 5 \
        "https://github.com/quarto-dev/quarto-cli/releases/download/v${QUARTO_VERSION}/quarto-${QUARTO_VERSION}-linux-amd64.deb" \
        --output /tmp/quarto.deb \
    && apt-get install -y --no-install-recommends /tmp/quarto.deb \
    && rm -f /tmp/quarto.deb \
    && rm -rf /var/lib/apt/lists/*

COPY containers/report/R-packages.tsv /tmp/dnaprs-report-R-packages.tsv

RUN Rscript -e "manifest <- read.delim('/tmp/dnaprs-report-R-packages.tsv', check.names=FALSE); for (row in seq_len(nrow(manifest))) remotes::install_version(manifest\$package[row], version=manifest\$install_version[row], repos='https://cloud.r-project.org', upgrade='never'); observed <- vapply(manifest\$package, function(package) as.character(packageVersion(package)), character(1)); if (any(observed != manifest\$expected_version)) stop(paste('Package version mismatch:', paste(manifest\$package[observed != manifest\$expected_version], collapse=', ')))" \
    && rm -f /tmp/dnaprs-report-R-packages.tsv

LABEL org.opencontainers.image.title="dnaprs report"
LABEL org.opencontainers.image.description="Pinned R and Quarto environment for the dnaprs report"
LABEL org.opencontainers.image.source="https://github.com/paulYRP/dnaprs"
LABEL org.opencontainers.image.version="1.0.0"
WORKDIR /work
ENTRYPOINT []
CMD ["/bin/bash"]
