FROM zhiliz/sbayesrc:0.2.6

ARG QUARTO_VERSION=1.7.32
ARG DEBIAN_FRONTEND=noninteractive

USER root

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        build-essential \
        ca-certificates \
        curl \
        gfortran \
        libcurl4-openssl-dev \
        libfontconfig1-dev \
        libnlopt-dev \
        libssl-dev \
        libtiff-dev \
        libuv1-dev \
        libxml2-dev \
    && curl --fail --location --retry 5 \
        "https://github.com/quarto-dev/quarto-cli/releases/download/v${QUARTO_VERSION}/quarto-${QUARTO_VERSION}-linux-amd64.deb" \
        --output /tmp/quarto.deb \
    && apt-get install -y --no-install-recommends /tmp/quarto.deb \
    && rm -f /tmp/quarto.deb \
    && rm -rf /var/lib/apt/lists/*

COPY containers/R-packages.tsv /tmp/dnaprs-R-packages.tsv

RUN Rscript -e "install.packages('remotes', repos='https://cloud.r-project.org'); manifest <- read.delim('/tmp/dnaprs-R-packages.tsv', check.names=FALSE); install <- manifest[manifest\$source == 'CRAN', ]; for (row in seq_len(nrow(install))) { package <- install\$package[row]; expected <- install\$expected_version[row]; current <- if (requireNamespace(package, quietly=TRUE)) as.character(packageVersion(package)) else ''; if (current != expected) remotes::install_version(package, version=install\$install_version[row], repos='https://cloud.r-project.org', upgrade='never') }; observed <- vapply(manifest\$package, function(package) if (requireNamespace(package, quietly=TRUE)) as.character(packageVersion(package)) else '', character(1)); if (any(observed != manifest\$expected_version)) stop(paste('Package version mismatch:', paste(manifest\$package[observed != manifest\$expected_version], collapse=', ')))" \
    && rm -f /tmp/dnaprs-R-packages.tsv

LABEL org.opencontainers.image.title="dnaprs"
LABEL org.opencontainers.image.description="Pinned PLINK, SBayesRC, R modelling, and Quarto environment for dnaprs"
LABEL org.opencontainers.image.version="1.0.0"
WORKDIR /work
ENTRYPOINT []
CMD ["/bin/bash"]
