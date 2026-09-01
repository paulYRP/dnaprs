FROM ubuntu:24.04@sha256:33ceb71981b602c1a7443a53469e4dba065f7503eab3078a2d7a57a2ab987517

ARG DEBIAN_FRONTEND=noninteractive
ARG PLINK2_URL=https://s3.amazonaws.com/plink2-assets/alpha6/plink2_linux_avx2_20250420.zip
ARG PLINK2_SHA256=4032c845da80ecaf55603861bf8a11424b87248146663ca7e11958c3ccff7d1c

RUN apt-get update \
    && apt-get install -y --no-install-recommends ca-certificates curl unzip \
    && curl --fail --location --retry 5 "${PLINK2_URL}" --output /tmp/plink2.zip \
    && echo "${PLINK2_SHA256}  /tmp/plink2.zip" | sha256sum --check --strict \
    && unzip /tmp/plink2.zip plink2 -d /usr/local/bin \
    && chmod 0755 /usr/local/bin/plink2 \
    && rm -f /tmp/plink2.zip \
    && rm -rf /var/lib/apt/lists/* \
    && plink2 --version

LABEL org.opencontainers.image.title="dnaprs PLINK 2"
LABEL org.opencontainers.image.description="PLINK 2.0 alpha 6.12 environment for dnaprs"
LABEL org.opencontainers.image.source="https://github.com/paulYRP/dnaprs"
LABEL org.opencontainers.image.version="2.0.0-a.6.12"
WORKDIR /work
ENTRYPOINT []
CMD ["/bin/bash"]
