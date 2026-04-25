# Reproducible benchmarks runner.
# Pinned to micromamba + bioconda so reviewers don't have to install
# competitor tools by hand.
#
# Build:
#   docker build -t gfc-bench -f benchmarks/Dockerfile .
# Run:
#   docker run --rm -v $(pwd)/benchmarks/results:/work/benchmarks/results gfc-bench

FROM mambaorg/micromamba:1.5.8

USER root
WORKDIR /work
COPY --chown=mambauser:mambauser . /work

USER mambauser
RUN micromamba install -y -n base -f /work/benchmarks/environment.yml \
 && micromamba clean --all --yes

# Activate the env by default.
SHELL ["/bin/bash", "-c"]
ENV MAMBA_DOCKERFILE_ACTIVATE=1 \
    PYTHONUNBUFFERED=1

ENTRYPOINT ["bash", "/work/benchmarks/run_bench.sh"]
