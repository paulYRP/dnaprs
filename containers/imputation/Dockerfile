FROM rocker/r-ver:4.4.1@sha256:78cb94ce2db23aaaf7b546450fcf70b5a3f2ace5a9b5fa1f87217da329211312

ARG DEBIAN_FRONTEND=noninteractive
ARG PLINK2_URL=https://s3.amazonaws.com/plink2-assets/alpha6/plink2_linux_avx2_20250420.zip
ARG PLINK2_SHA256=4032c845da80ecaf55603861bf8a11424b87248146663ca7e11958c3ccff7d1c
ARG BEAGLE_URL=https://faculty.washington.edu/browning/beagle/beagle.27Feb25.75f.jar
ARG BEAGLE_SHA256=7319f4af9638be05c18dcc1bfb8fb41a58a09293507ebf0d54617d0e40df5a70
ARG UNBREF3_URL=https://faculty.washington.edu/browning/beagle/unbref3.27Feb25.75f.jar
ARG UNBREF3_SHA256=db8998ff212b19f72fd0e36e8cd9ef9b4e414d6a6c295fe89ee77858970cfb4a
ARG DATA_TABLE_URL=https://packagemanager.posit.co/cran/2026-08-01/src/contrib/Archive/data.table/data.table_1.18.0.tar.gz
ARG DATA_TABLE_SHA256=da2cde0c5b9718fc093371f71aa06b717ffe3c6569b83e03d9a9d2f78e80f66a
ARG VCS_REF=unknown
ARG BUILD_DATE=unknown

RUN apt-get update \
    && apt-get install -y --no-install-recommends bcftools build-essential ca-certificates curl openjdk-17-jre-headless tabix unzip \
    && curl --fail --location --retry 5 "${PLINK2_URL}" --output /tmp/plink2.zip \
    && echo "${PLINK2_SHA256}  /tmp/plink2.zip" | sha256sum --check --strict \
    && unzip /tmp/plink2.zip plink2 -d /usr/local/bin \
    && chmod 0755 /usr/local/bin/plink2 \
    && mkdir -p /opt/beagle \
    && curl --fail --location --retry 5 "${BEAGLE_URL}" --output /opt/beagle/beagle.jar \
    && curl --fail --location --retry 5 "${UNBREF3_URL}" --output /opt/beagle/unbref3.jar \
    && echo "${BEAGLE_SHA256}  /opt/beagle/beagle.jar" | sha256sum --check --strict \
    && echo "${UNBREF3_SHA256}  /opt/beagle/unbref3.jar" | sha256sum --check --strict \
    && curl --fail --location --retry 5 "${DATA_TABLE_URL}" --output /tmp/data.table.tar.gz \
    && echo "${DATA_TABLE_SHA256}  /tmp/data.table.tar.gz" | sha256sum --check --strict \
    && R CMD INSTALL /tmp/data.table.tar.gz \
    && Rscript -e "stopifnot(as.character(packageVersion('data.table')) == '1.18.0')" \
    && rm -f /tmp/plink2.zip /tmp/data.table.tar.gz \
    && rm -rf /var/lib/apt/lists/* \
    && plink2 --version \
    && bcftools --version | head -n 1 \
    && java -jar /opt/beagle/beagle.jar 2>&1 | head -n 2 \
    && java -jar /opt/beagle/unbref3.jar help 2>&1 | head -n 2

ENV BEAGLE_JAR=/opt/beagle/beagle.jar
ENV UNBREF3_JAR=/opt/beagle/unbref3.jar
LABEL org.opencontainers.image.title="dnaprs target imputation"
LABEL org.opencontainers.image.description="Pinned Beagle, bcftools, tabix, Java, PLINK 2, R, and data.table environment"
LABEL org.opencontainers.image.source="https://github.com/paulYRP/dnaprs"
LABEL org.opencontainers.image.version="1.1.0"
LABEL org.opencontainers.image.revision="${VCS_REF}"
LABEL org.opencontainers.image.created="${BUILD_DATE}"
LABEL org.opencontainers.image.licenses="MIT"
WORKDIR /work
ENTRYPOINT []
CMD ["/bin/bash"]
