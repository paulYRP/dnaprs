FROM rocker/r-ver:4.4.1@sha256:78cb94ce2db23aaaf7b546450fcf70b5a3f2ace5a9b5fa1f87217da329211312

ARG DEBIAN_FRONTEND=noninteractive
ARG CRAN_SNAPSHOT=https://packagemanager.posit.co/cran/2026-08-01
ARG REMOTES_URL=https://cloud.r-project.org/src/contrib/remotes_2.5.0.tar.gz
ARG REMOTES_SHA256=4d663f1426cd88d42f4070f23d969305c575e0499ed1397be6607b0770d2850c
ARG VCS_REF=unknown
ARG BUILD_DATE=unknown

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        build-essential \
        ca-certificates \
        curl \
        gfortran \
        libcurl4-openssl-dev \
        libnlopt-dev \
        libssl-dev \
        libxml2-dev \
        pkg-config \
        zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

COPY containers/analysis/R-packages.tsv /tmp/dnaprs-R-packages.tsv

RUN curl --fail --location --retry 5 "${REMOTES_URL}" --output /tmp/remotes.tar.gz \
    && echo "${REMOTES_SHA256}  /tmp/remotes.tar.gz" | sha256sum --check --strict \
    && Rscript -e "install.packages('/tmp/remotes.tar.gz', repos=NULL, type='source'); manifest <- read.delim('/tmp/dnaprs-R-packages.tsv', check.names=FALSE); for (row in seq_len(nrow(manifest))) remotes::install_version(manifest\$package[row], version=manifest\$install_version[row], repos='${CRAN_SNAPSHOT}', upgrade='never'); observed <- vapply(manifest\$package, function(package) as.character(packageVersion(package)), character(1)); if (any(observed != manifest\$expected_version)) stop(paste('Package version mismatch:', paste(manifest\$package[observed != manifest\$expected_version], collapse=', ')))" \
    && rm -f /tmp/remotes.tar.gz \
    && rm -f /tmp/dnaprs-R-packages.tsv

LABEL org.opencontainers.image.title="dnaprs analysis"
LABEL org.opencontainers.image.description="Pinned R environment for dnaprs validation, QC, score processing, and phenotype models"
LABEL org.opencontainers.image.source="https://github.com/paulYRP/dnaprs"
LABEL org.opencontainers.image.version="1.0.0"
LABEL org.opencontainers.image.revision="${VCS_REF}"
LABEL org.opencontainers.image.created="${BUILD_DATE}"
LABEL org.opencontainers.image.licenses="MIT"
WORKDIR /work
ENTRYPOINT []
CMD ["/bin/bash"]
