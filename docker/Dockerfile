FROM rocker/rstudio:4.5.2

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
    ca-certificates \
    cmake \
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
    libnlopt-dev \
    libuv1-dev \
    libzmq3-dev \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /opt/hcocena
COPY . /opt/hcocena

RUN R -q -e "options(timeout = 1200); install.packages(c('remotes', 'BiocManager', 'graphlayouts', 'ellmer'), repos = 'https://cloud.r-project.org')" && \
    R -q -e "options(timeout = 1200); repos <- BiocManager::repositories(); repos['CRAN'] <- 'https://cloud.r-project.org'; options(repos = repos); pkgs <- c('ape', 'treeio', 'ggtree', 'GO.db', 'GOSemSim', 'DOSE', 'enrichplot', 'clusterProfiler', 'CALIBERrfimpute', 'decoupleR', 'dorothea', 'mice', 'progeny', 'RCy3', 'SpatialExperiment', 'GSVA'); for (attempt in seq_len(3)) { missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]; if (!length(missing)) break; message('Docker dependency install attempt ', attempt, ': ', paste(missing, collapse = ', ')); BiocManager::install(missing, ask = FALSE, update = FALSE) }; status <- vapply(pkgs, requireNamespace, logical(1), quietly = TRUE); print(status); if (!all(status)) stop('Missing Docker dependency packages: ', paste(names(status)[!status], collapse = ', '))" && \
    R -q -e "options(timeout = 1200); repos <- BiocManager::repositories(); repos['CRAN'] <- 'https://cloud.r-project.org'; options(repos = repos); remotes::install_local('/opt/hcocena', dependencies = NA, upgrade = 'never')" && \
    R -q -e "status <- vapply(c('graphlayouts', 'ellmer', 'clusterProfiler', 'CALIBERrfimpute', 'decoupleR', 'dorothea', 'mice', 'progeny', 'RCy3', 'SpatialExperiment', 'GSVA', 'hcocena'), requireNamespace, logical(1), quietly = TRUE); print(status); if (!all(status)) stop('Missing Docker packages: ', paste(names(status)[!status], collapse = ', '))" && \
    R -q -e "packageVersion('hcocena')"

RUN mkdir -p /home/rstudio/count_data && \
    mkdir -p /home/rstudio/annotation_data && \
    mkdir -p /home/rstudio/project_data && \
    mkdir -p /home/rstudio/output && \
    mkdir -p /home/rstudio/reference_files && \
    cp /opt/hcocena/docker/README.md /home/rstudio/README_docker.md && \
    cp /opt/hcocena/docker/project_workflows/hcocena_main.Rmd /home/rstudio/01_hcocena_main.Rmd && \
    cp /opt/hcocena/docker/project_workflows/hcocena_satellite.Rmd /home/rstudio/02_hcocena_satellite.Rmd && \
    cp /opt/hcocena/github_workflows/hcocena_main_seq_only.Rmd /home/rstudio/03_hcocena_main_seq_only.Rmd && \
    cp -r /opt/hcocena/docker/reference_files/. /home/rstudio/reference_files/ && \
    cp /opt/hcocena/docker/project_data/README.md /home/rstudio/project_data/README.md && \
    chown -R rstudio:rstudio /home/rstudio

WORKDIR /home/rstudio

EXPOSE 8787
CMD ["/init"]
