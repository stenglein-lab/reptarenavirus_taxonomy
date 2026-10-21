# ---------------------------------------------------------------------------
#
# This dockerfile defines a docker/singularity image needed to reproduce
# the paper described in this repository.
#
# It handles software dependencies, namely:
#
# - R
# - R packages, like tidyverse, patchwork, etc
# - quarto
# 
# In the interest of reproducibility, it uses:
# - Pinned R version via rocker/r-ver
# - Pinned package versions via a fixed Posit Package Manager (PPM) CRAN
#   snapshot date, installed as fast pre-compiled binaries 
# - Quarto CLI installed directly 
#
# see README.me for build and push instructions
#   
# ---------------------------------------------------------------------------

# the base rocker R image
FROM rocker/r-ver:4.5.0
# --- Pin the CRAN snapshot date for full package-version reproducibility ---
# Every install.packages() call below will resolve against packages exactly
# as they existed on this date. Change this date deliberately when you want
# to bump versions; otherwise builds stay reproducible indefinitely.
ENV CRAN_SNAPSHOT_DATE=2026-08-01
RUN echo "options(repos = c(CRAN = 'https://packagemanager.posit.co/cran/__linux__/noble/${CRAN_SNAPSHOT_DATE}'))" \
    >> /usr/local/lib/R/etc/Rprofile.site

# --- System dependencies commonly needed by tidyverse / plotting packages ---
RUN apt-get update && apt-get install -y --no-install-recommends \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    libuv1-dev \
    libglpk40 \
    pandoc \
    curl \
    gdebi-core \
    && rm -rf /var/lib/apt/lists/*

# --- Install Quarto CLI ---
# Check https://quarto.org/docs/get-started/ for the current release version
ARG QUARTO_VERSION=1.10.18
RUN curl -LO https://github.com/quarto-dev/quarto-cli/releases/download/v${QUARTO_VERSION}/quarto-${QUARTO_VERSION}-linux-amd64.deb \
    && gdebi --non-interactive quarto-${QUARTO_VERSION}-linux-amd64.deb \
    && rm quarto-${QUARTO_VERSION}-linux-amd64.deb

# --- Install R packages (pinned by the CRAN_SNAPSHOT_DATE set above) ---
RUN Rscript -e "install.packages(c( \
    'tidyverse', \
    'rentrez', \
    'readxl', \
    'BiocManager', \
    'phytools', \
    'ape', \
    'seqinr', \
    'RColorBrewer', \
    'patchwork' \
    ))"

# --- Install R packages using Bioconductor, specifying versions to install
RUN Rscript -e "BiocManager::install('ggtree', version = '3.22')"
RUN Rscript -e "BiocManager::install('MSA2dist', version = '3.22')"

# --- Sanity check that everything loads at build time (fail fast) ---
RUN Rscript -e "library(tidyverse); library(rentrez); library(readxl); \
    library(ggtree); library(phytools); library(ape); library(seqinr); \
    library(RColorBrewer); library(MSA2dist); library(patchwork); \
    cat('All packages loaded OK\n')"

WORKDIR /workflow
CMD ["R"]
