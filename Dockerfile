# Multi-stage Dockerfile for genome-format-converters.
# Build:   docker build -t gfc:0.1.3 .
# Run:     docker run --rm -v $(pwd):/data gfc:0.1.3 gfc <subcommand> ...

FROM python:3.11-slim AS build

WORKDIR /src
COPY . /src

# pysam needs libz / libbz2 / liblzma / libcurl / libssl headers to compile.
# Install them only in the build stage and throw away in the runtime stage.
RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential \
        zlib1g-dev \
        libbz2-dev \
        liblzma-dev \
        libcurl4-openssl-dev \
        libssl-dev \
    && rm -rf /var/lib/apt/lists/*

RUN python -m pip install --upgrade pip build \
 && python -m build --wheel \
 && python -m pip install --prefix=/install dist/*.whl


FROM python:3.11-slim AS runtime

# htslib runtime dependencies for pysam.
RUN apt-get update && apt-get install -y --no-install-recommends \
        zlib1g \
        libbz2-1.0 \
        liblzma5 \
        libcurl4 \
        libssl3 \
    && rm -rf /var/lib/apt/lists/*

COPY --from=build /install /usr/local

WORKDIR /data
ENTRYPOINT ["gfc"]
CMD ["--help"]
