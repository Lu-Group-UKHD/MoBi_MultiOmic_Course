# Docker image for the course practicals:
# "Integrative Multi-Omic Data Analysis for Biological and Biomedical Data"
#
# Includes:
# - RStudio Server
# - R, Quarto/R Markdown support
# - CRAN packages used in the Day 1 and Day 2 practicals
# - Bioconductor packages: DESeq2, SummarizedExperiment, MOFA2
# - Python/mofapy2 backend for MOFA2
# - Course Quarto files and example datasets

FROM rocker/rstudio

LABEL maintainer="Junyan Lu"
LABEL description="RStudio Server image for Day 1 and Day 2 multi-omics practicals."

ENV DEBIAN_FRONTEND=noninteractive
ENV RETICULATE_PYTHON=/usr/bin/python3
ENV QUARTO_R=/usr/local/bin/R
ENV COURSE_DIR=/home/rstudio/course

# System libraries commonly needed by R/Bioconductor packages and plotting packages.
RUN apt-get update && apt-get install -y --no-install-recommends \
    python3 \
    python3-pip \
    python3-venv \
    python3-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libhdf5-dev \
    libgsl-dev \
    libgit2-dev \
    libudunits2-dev \
    libfontconfig1-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libcairo2-dev \
    libxt-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Python backend for MOFA2.
# MOFA2 calls the Python package mofapy2 for model training.
RUN python3 -m pip install --no-cache-dir --break-system-packages \
    numpy \
    scipy \
    pandas \
    scikit-learn \
    h5py \
    mofapy2

# CRAN packages used by the practicals.
RUN R -e "install.packages(c('BiocManager', 'data.table', 'tidyverse', 'matrixStats', 'knitr', 'ggpubr', 'randomForest', 'survival', 'survminer', 'reticulate', 'rmarkdown', 'UpSetR','psych'), ask=FALSE)"

# Bioconductor packages used by Day 1 and Day 2.
RUN R -e "options(timeout = 1000); BiocManager::install(c('DESeq2', 'SummarizedExperiment', 'MOFA2'), ask = FALSE, update = FALSE, Ncpus = max(1, parallel::detectCores() - 1))"

# Make sure reticulate sees the system Python with mofapy2 installed.
RUN echo 'RETICULATE_PYTHON=/usr/bin/python3' >> /usr/local/lib/R/etc/Renviron.site

# Course material and datasets.
COPY course/ ${COURSE_DIR}/

RUN chown -R rstudio:rstudio ${COURSE_DIR}

WORKDIR ${COURSE_DIR}

# RStudio Server is exposed on port 8787 by the rocker/rstudio base image.
EXPOSE 8787



