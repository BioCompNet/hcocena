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
    libmagick++-dev \
    libzmq3-dev \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /opt/hcocena
COPY . /opt/hcocena

RUN R -q -e "install.packages(c('remotes', 'BiocManager'), repos='https://cloud.r-project.org')" && \
    R -q -e "remotes::install_local('/opt/hcocena', dependencies = NA, upgrade = 'never')" && \
    R -q -e "BiocManager::install(c('CALIBERrfimpute', 'decoupleR', 'dorothea', 'mice', 'progeny', 'RCy3', 'SpatialExperiment', 'GSVA'), ask = FALSE, update = FALSE)" && \
    R -q -e "status <- vapply(c('CALIBERrfimpute', 'decoupleR', 'dorothea', 'mice', 'progeny', 'RCy3', 'SpatialExperiment', 'GSVA'), requireNamespace, logical(1), quietly = TRUE); print(status); if (!all(status)) stop('Missing optional Docker packages: ', paste(names(status)[!status], collapse = ', '))" && \
    R -q -e "packageVersion('hcocena')"

RUN mkdir -p /home/rstudio/hcocena && \
    mkdir -p /home/rstudio/hcocena/count_data && \
    mkdir -p /home/rstudio/hcocena/annotation_data && \
    mkdir -p /home/rstudio/hcocena/output && \
    cp /opt/hcocena/README.md /home/rstudio/hcocena/ && \
    cp -r /opt/hcocena/docker/reference_files /home/rstudio/hcocena/reference_files && \
    chown -R rstudio:rstudio /home/rstudio

WORKDIR /home/rstudio/hcocena

EXPOSE 8787
CMD ["/init"]
