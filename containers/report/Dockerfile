FROM rocker/r-ver:4.4.1@sha256:78cb94ce2db23aaaf7b546450fcf70b5a3f2ace5a9b5fa1f87217da329211312

ARG QUARTO_VERSION=1.7.32
ARG QUARTO_SHA256=9abee72a070a0735fe30c6824e1ca99dae3c118370217357c16ad4f5d6056071
ARG CRAN_SNAPSHOT=https://packagemanager.posit.co/cran/2026-08-01
ARG REMOTES_URL=https://cloud.r-project.org/src/contrib/remotes_2.5.0.tar.gz
ARG REMOTES_SHA256=4d663f1426cd88d42f4070f23d969305c575e0499ed1397be6607b0770d2850c
ARG DEBIAN_FRONTEND=noninteractive
ARG VCS_REF=unknown
ARG BUILD_DATE=unknown

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        build-essential ca-certificates curl gfortran libcurl4-openssl-dev \
        libfontconfig1-dev libnlopt-dev libssl-dev libtiff-dev libuv1-dev \
        libxml2-dev pkg-config zlib1g-dev \
    && curl --fail --location --retry 5 \
        "https://github.com/quarto-dev/quarto-cli/releases/download/v${QUARTO_VERSION}/quarto-${QUARTO_VERSION}-linux-amd64.deb" \
        --output /tmp/quarto.deb \
    && echo "${QUARTO_SHA256}  /tmp/quarto.deb" | sha256sum --check --strict \
    && apt-get install -y --no-install-recommends /tmp/quarto.deb \
    && curl --fail --location --retry 5 "${REMOTES_URL}" --output /tmp/remotes.tar.gz \
    && echo "${REMOTES_SHA256}  /tmp/remotes.tar.gz" | sha256sum --check --strict \
    && rm -f /tmp/quarto.deb \
    && rm -rf /var/lib/apt/lists/*

COPY containers/analysis/R-packages.tsv /tmp/dnaprs-analysis-R-packages.tsv
COPY containers/report/R-packages.tsv /tmp/dnaprs-report-R-packages.tsv
COPY containers/install-r-packages.R /tmp/install-r-packages.R

RUN Rscript -e "install.packages('/tmp/remotes.tar.gz', repos=NULL, type='source')" \
    && Rscript /tmp/install-r-packages.R "${CRAN_SNAPSHOT}" \
        /tmp/dnaprs-analysis-R-packages.tsv /tmp/dnaprs-report-R-packages.tsv \
    && rm -f /tmp/remotes.tar.gz /tmp/dnaprs-analysis-R-packages.tsv \
        /tmp/dnaprs-report-R-packages.tsv /tmp/install-r-packages.R

LABEL org.opencontainers.image.title="dnaprs report"
LABEL org.opencontainers.image.description="Pinned R and Quarto environment for the dnaprs report"
LABEL org.opencontainers.image.source="https://github.com/paulYRP/dnaprs"
LABEL org.opencontainers.image.version="1.0.0"
LABEL org.opencontainers.image.revision="${VCS_REF}"
LABEL org.opencontainers.image.created="${BUILD_DATE}"
LABEL org.opencontainers.image.licenses="MIT"
WORKDIR /work
ENTRYPOINT []
CMD ["/bin/bash"]
