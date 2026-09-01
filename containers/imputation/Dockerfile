FROM ghcr.io/paulyrp/dnaprs-plink2:2.0.0-a.6.12

ARG DEBIAN_FRONTEND=noninteractive
ARG BEAGLE_URL=https://faculty.washington.edu/browning/beagle/beagle.27Feb25.75f.jar
ARG BEAGLE_SHA256=7319f4af9638be05c18dcc1bfb8fb41a58a09293507ebf0d54617d0e40df5a70
ARG UNBREF3_URL=https://faculty.washington.edu/browning/beagle/unbref3.27Feb25.75f.jar
ARG UNBREF3_SHA256=db8998ff212b19f72fd0e36e8cd9ef9b4e414d6a6c295fe89ee77858970cfb4a

RUN apt-get update \
    && apt-get install -y --no-install-recommends bcftools openjdk-17-jre-headless tabix \
    && mkdir -p /opt/beagle \
    && curl --fail --location --retry 5 "${BEAGLE_URL}" --output /opt/beagle/beagle.jar \
    && curl --fail --location --retry 5 "${UNBREF3_URL}" --output /opt/beagle/unbref3.jar \
    && echo "${BEAGLE_SHA256}  /opt/beagle/beagle.jar" | sha256sum --check --strict \
    && echo "${UNBREF3_SHA256}  /opt/beagle/unbref3.jar" | sha256sum --check --strict \
    && rm -rf /var/lib/apt/lists/* \
    && bcftools --version | head -n 1 \
    && java -jar /opt/beagle/beagle.jar 2>&1 | head -n 2 \
    && java -jar /opt/beagle/unbref3.jar help 2>&1 | head -n 2

ENV BEAGLE_JAR=/opt/beagle/beagle.jar
ENV UNBREF3_JAR=/opt/beagle/unbref3.jar
LABEL org.opencontainers.image.title="dnaprs target imputation"
LABEL org.opencontainers.image.description="Pinned Beagle, bcftools, tabix, Java, and PLINK 2 environment"
LABEL org.opencontainers.image.source="https://github.com/paulYRP/dnaprs"
LABEL org.opencontainers.image.version="1.1.0"
WORKDIR /work
ENTRYPOINT []
CMD ["/bin/bash"]
