# Base image: Bioconductor (R 4.5 + RStudio)
FROM bioconductor/bioconductor_docker:RELEASE_3_23

ENV LANG=C.UTF-8
ENV LC_ALL=C.UTF-8

# Environment variables
ENV CONDA_DIR=/opt/conda
ENV PATH=${CONDA_DIR}/bin:/root/.cargo/bin:$PATH
ENV MAKE="make"
ENV MAKEFLAGS="-j4"

# System dependencies (combined original + extra rpy2/plotting deps)
RUN apt-get update && apt-get install -y --no-install-recommends \
    curl \
    git \
    build-essential \
    gcc \
    g++ \
    gfortran \
    make \
    cmake \
    pkg-config \
    libpcre2-dev \
    libdeflate-dev \
    liblzma-dev \
    zlib1g-dev \
    libbz2-dev \
    libzstd-dev \
    libicu-dev \
    libtirpc-dev \
    libxml2-dev \
    libssl-dev \
    libcurl4-openssl-dev \
    libv8-dev \
    libhdf5-dev \
    libffi-dev \
    libreadline-dev \
    libsqlite3-dev \
    tk-dev \
    libpng-dev \
    libjpeg-dev \
    libfreetype6-dev \
    fontconfig \
    fonts-dejavu-core \
    fonts-dejavu-extra \
    fonts-liberation \
    fonts-freefont-ttf \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

# Refresh font cache
RUN fc-cache -fv

# Switch user so Rust installs into /home/rstudio/.cargo
USER rstudio

# Install Rust
RUN curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y

# Make cargo available in shell
ENV PATH="/home/rstudio/.cargo/bin:${PATH}"

# Make cargo available to R sessions (.Renviron is read by R)
RUN echo 'PATH=${PATH}:/home/rstudio/.cargo/bin' >> /home/rstudio/.Renviron

# Switch back to root
USER root

### Rust for snapATAC2

RUN curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | \
    sh -s -- -y --default-toolchain stable

ENV CARGO_HOME="/root/.cargo"
ENV RUSTUP_HOME="/root/.rustup"
ENV PATH="/root/.cargo/bin:${PATH}"

RUN rustc --version && cargo --version

###

# Install Miniforge3
RUN curl -L https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh -o miniforge.sh \
    && bash miniforge.sh -b -p ${CONDA_DIR} \
    && rm miniforge.sh \
    && conda config --set auto_update_conda False

# Install uv from official image
COPY --from=ghcr.io/astral-sh/uv:latest /uv /uvx /usr/bin/

# Configure R + reticulate to use uv-managed environments
ENV RETICULATE_PYTHON_PREFERENCE=system
ENV RETICULATE_USE_MANAGED_VENV=yes

# Pre-install core R packages
RUN R -e "BiocManager::install(c('scran', 'scater', 'muscat', 'Rhdf5lib', 'edgeR', 'limma', 'glmGamPoi', 'DESeq2'), ask=FALSE)"
RUN R -e "BiocManager::install(c('basilisk', 'SingleCellExperiment', 'SparseArray', 'scrapper', 'sketchR'), ask=FALSE)"

RUN R -e "BiocManager::install(c('alabaster.base', 'alabaster.sce'), ask=FALSE)"

RUN R -e "BiocManager::install(c('lemur', 'MOFA2'), ask=FALSE)"
RUN R -e "BiocManager::install(c('decoupleR', 'clusterProfiler'), ask=FALSE)"
RUN R -e "BiocManager::install(c('ggplot2', 'scattermore', 'patchwork', 'cowplot', 'ComplexHeatmap', 'pheatmap', 'ggrastr'), ask=FALSE)"
RUN R -e "BiocManager::install(c('ggsci'), ask=FALSE)"
RUN R -e "BiocManager::install(c('reticulate', 'anndataR'), ask=FALSE)"

RUN R -e "BiocManager::install(c('tidyverse'), ask=FALSE)"


# Required installs to use R in Jupyter notebooks:

# Update and install the jupyter-client binary
RUN apt-get update
RUN apt-get install -y jupyter-client


RUN R -e "install.packages('IRkernel')"
RUN R -e "IRkernel::installspec(user = FALSE)"


RUN R -e "BiocManager::install(c('airway', 'TENxPBMCData'))"

WORKDIR /home/rstudio/project