FROM rocker/rstudio:4.5.2

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
    ca-certificates \
    gcc \
    g++ \
    gfortran \
    make \
    libblas-dev \
    liblapack-dev \
    libbz2-dev \
    liblzma-dev \
    zlib1g-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libcairo2-dev \
    libxt-dev \
    libfontconfig1-dev \
    libfreetype6-dev \
    libfribidi-dev \
    libharfbuzz-dev \
    libpng-dev \
    libjpeg-dev \
    libtiff5-dev \
    libgit2-dev \
    libglpk-dev \
    libgmp-dev \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /opt/hcocena
COPY . /opt/hcocena

# Install hcocena runtime dependencies + optional upstream inference packages.
RUN R -q -e "install.packages(c('remotes','BiocManager'), repos='https://cloud.r-project.org')" && \
    R -q -e "BiocManager::install(c('clusterProfiler','ComplexHeatmap','MultiAssayExperiment','SummarizedExperiment','S4Vectors','decoupleR','dorothea','progeny','RCy3'), ask=FALSE, update=FALSE)" && \
    R -q -e "remotes::install_local('/opt/hcocena/hCoCena-r-package', dependencies=TRUE, upgrade='never')" && \
    R -q -e "packageVersion('hcocena')"

# Provide ready-to-use workflow material in a writable workspace.
RUN mkdir -p /workspace && \
    cp /opt/hcocena/hcocena_main.Rmd /workspace/ && \
    cp /opt/hcocena/hcocena_main_seq_only.Rmd /workspace/ && \
    cp /opt/hcocena/hcocena_satellite.Rmd /workspace/ && \
    cp -r /opt/hcocena/reference_files /workspace/reference_files && \
    cp -r /opt/hcocena/scripts /workspace/scripts && \
    chown -R rstudio:rstudio /workspace

EXPOSE 8787
CMD ["/init"]

