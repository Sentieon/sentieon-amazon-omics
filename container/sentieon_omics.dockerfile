FROM amazonlinux:2023.12.20260817.0 AS builder
SHELL ["/bin/bash", "-o", "pipefail", "-c"]

ARG SENTIEON_VERSION
RUN test -n "$SENTIEON_VERSION"

# Build dependencies
RUN dnf -y upgrade && \
    dnf -y install \
      autoconf \
      automake \
      libtool \
      binutils \
      make \
      gcc \
      nasm \
      tar \
      gzip \
      bzip2 \
      xz \
      zlib-devel \
      bzip2-devel \
      xz-devel \
      libcurl-devel \
      openssl-devel \
      ncurses-devel \
      perl \
      perl-Data-Dumper && \
    dnf clean all

# Install igzip
RUN mkdir -p /opt/isa-l && \
    curl -fL --retry 3 "https://github.com/intel/isa-l/archive/refs/tags/v2.30.0.tar.gz" | \
        tar -C /opt/isa-l -zxf - && \
    cd /opt/isa-l/isa-l-2.30.0 && \
    ./autogen.sh && \
    ./configure --prefix=/usr --libdir=/usr/lib64 && \
    make -j"$(nproc)" install

# Install samtools
RUN mkdir -p /opt/samtools/ && \
    curl -fL --retry 3 "https://github.com/samtools/samtools/releases/download/1.16.1/samtools-1.16.1.tar.bz2" | \
      tar -C /opt/samtools/ -jxf - && \
    cd /opt/samtools/samtools-1.16.1 && \
    ./configure && \
    make -j"$(nproc)" install

# Install bcftools
RUN mkdir -p /opt/bcftools/ && \
    curl -fL --retry 3 "https://github.com/samtools/bcftools/releases/download/1.16/bcftools-1.16.tar.bz2" | \
      tar -C /opt/bcftools/ -jxf - && \
    cd /opt/bcftools/bcftools-1.16/ && \
    ./configure && \
    make -j"$(nproc)" install

# Install bedtools
RUN curl -fL --retry 3 -o /usr/local/bin/bedtools-2.30.0 \
      "https://github.com/arq5x/bedtools2/releases/download/v2.30.0/bedtools.static.binary"

# Install jemalloc
RUN mkdir -p /opt/jemalloc/ && \
    curl -fL --retry 3 "https://github.com/jemalloc/jemalloc/releases/download/5.3.0/jemalloc-5.3.0.tar.bz2" | \
      tar -C /opt/jemalloc --no-same-owner -jxf - && \
    cd /opt/jemalloc/jemalloc-5.3.0/ && \
    ./configure && \
    make -j"$(nproc)" && \
    make install

# Install sentieon
RUN mkdir -p /opt/sentieon/ && \
    curl -fL --retry 3 "https://s3.amazonaws.com/sentieon-release/software/sentieon-genomics-${SENTIEON_VERSION}.tar.gz" | \
      tar -C /opt/sentieon -zxf -

# Build the container
FROM amazonlinux:2023.12.20260817.0
SHELL ["/bin/bash", "-o", "pipefail", "-c"]

ARG SENTIEON_VERSION
ENV SENTIEON_VERSION=$SENTIEON_VERSION

LABEL container.base.image="amazonlinux:2023.12.20260817.0" \
      software.version="${SENTIEON_VERSION}" \
      software.website="https://www.sentieon.com/"

# Install dependencies
RUN dnf -y upgrade && \
    dnf -y install \
      awscli-2 \
      jq \
      util-linux \
      perl \
      tar \
      gzip \
      which \
      procps-ng \
      findutils \
      zlib \
      bzip2-libs \
      xz-libs \
      libcurl-minimal \
      openssl-libs \
      ncurses-libs && \
    dnf clean all && rm -rf /var/cache/dnf

# Copy dependencies from the first stage
COPY --from=builder /opt/sentieon/sentieon-genomics-${SENTIEON_VERSION} /opt/sentieon/sentieon-genomics-${SENTIEON_VERSION}
COPY --from=builder /usr/bin/igzip /usr/bin/igzip
COPY --from=builder /usr/lib64/libisal.so.2.0.30 /usr/lib64/libisal.so.2.0.30
COPY --from=builder /usr/local/bin/samtools /usr/local/bin/samtools
COPY --from=builder /usr/local/bin/bcftools /usr/local/bin/bcftools
COPY --from=builder /usr/local/bin/bedtools-2.30.0 /usr/local/bin/bedtools-2.30.0
COPY --from=builder /usr/local/lib/libjemalloc.so.2 /usr/local/lib/libjemalloc.so.2

# Create links
RUN ldconfig && \
    ln -sf libisal.so.2 /usr/lib64/libisal.so && \
    ln -sf libjemalloc.so.2 /usr/local/lib/libjemalloc.so && \
    ln -sf bedtools-2.30.0 /usr/local/bin/bedtools && \
    chmod ugo+x /usr/local/bin/bedtools-2.30.0

ENV SENTIEON_INSTALL_DIR=/opt/sentieon/sentieon-genomics-${SENTIEON_VERSION}
ENV PATH=$SENTIEON_INSTALL_DIR/bin:$PATH
ENV LD_PRELOAD=/usr/local/lib/libjemalloc.so.2

# A default jemalloc configuration that should work well for most use-cases, see http://jemalloc.net/jemalloc.3.html
ENV MALLOC_CONF=metadata_thp:auto,background_thread:true,dirty_decay_ms:30000,muzzy_decay_ms:30000

# Test the container
RUN sentieon driver --help && \
    igzip --help && \
    samtools --version && \
    samtools samples < /dev/null && \
    bcftools --version && \
    bedtools --version && \
    aws --version && \
    curl --version && \
    jq --version && \
    perl -MFcntl -e 1 && \
    taskset --version && \
    lscpu > /dev/null && \
    ps --version && \
    sh -c 'cat /proc/self/maps' > /tmp/jemalloc_preload_check && \
    grep -q jemalloc /tmp/jemalloc_preload_check && \
    rm -f /tmp/jemalloc_preload_check

CMD ["/bin/bash"]
