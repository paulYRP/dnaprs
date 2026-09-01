FROM rocker/r-ver:4.4.1@sha256:78cb94ce2db23aaaf7b546450fcf70b5a3f2ace5a9b5fa1f87217da329211312

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        build-essential \
        ca-certificates \
        gfortran \
        libcurl4-openssl-dev \
        libnlopt-dev \
        libssl-dev \
        libxml2-dev \
        pkg-config \
        zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

COPY containers/analysis/R-packages.tsv /tmp/dnaprs-R-packages.tsv

RUN Rscript -e "install.packages('https://cloud.r-project.org/src/contrib/remotes_2.5.0.tar.gz', repos=NULL, type='source'); manifest <- read.delim('/tmp/dnaprs-R-packages.tsv', check.names=FALSE); for (row in seq_len(nrow(manifest))) remotes::install_version(manifest\$package[row], version=manifest\$install_version[row], repos='https://cloud.r-project.org', upgrade='never'); observed <- vapply(manifest\$package, function(package) as.character(packageVersion(package)), character(1)); if (any(observed != manifest\$expected_version)) stop(paste('Package version mismatch:', paste(manifest\$package[observed != manifest\$expected_version], collapse=', ')))" \
    && rm -f /tmp/dnaprs-R-packages.tsv

LABEL org.opencontainers.image.title="dnaprs analysis"
LABEL org.opencontainers.image.description="Pinned R environment for dnaprs validation, QC, score processing, and phenotype models"
LABEL org.opencontainers.image.source="https://github.com/paulYRP/dnaprs"
LABEL org.opencontainers.image.version="1.0.0"
WORKDIR /work
ENTRYPOINT []
CMD ["/bin/bash"]
