FROM r-base

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
    build-essential \
    libcurl4-gnutls-dev \
    libxml2-dev \
    libssl-dev \
    libgit2-dev \
    libcairo2-dev \
    libxt-dev \
    xvfb \
    xauth \
    xfonts-base \
    libssl1.1

# install R packages available in CRAN
RUN R -e 'install.packages(c("dplyr", "devtools", "BiocManager", "plot3D"))'

# install R packages via Bioconductor
RUN R -e 'BiocManager::install(c("biomaRt", "GSVA", "rhdf5", "ComplexHeatmap", "ConsensusClusterPlus", "DESeq2", "tximport", "impute", "limma", "GEOquery"))'

# install NetBID2 from GitHub master branch
RUN xvfb-run R -e 'devtools::install_github("jyyulab/NetBID", ref="master", dependencies="Depends")'

# rstudio stable release version
ENV RSTUDIO_VERSION=1.2.5042
# ENV RSTUDIO_VERSION=1.3.1093

# rstudio daily build
# ENV RSTUDIO_VERSION=1.4.1100

# Install RStudio Server
RUN apt-get install -y --no-install-recommends \
    ca-certificates \
    wget \
    gdebi-core \
    && wget \
    --no-verbose \
    -O rstudio-server.deb \
    "https://download2.rstudio.org/server/trusty/amd64/rstudio-server-${RSTUDIO_VERSION}-amd64.deb" \
    # "https://download2.rstudio.org/server/xenial/amd64/rstudio-server-${RSTUDIO_VERSION}-amd64.deb" \
    # "https://s3.amazonaws.com/rstudio-ide-build/server/bionic/amd64/rstudio-server-${RSTUDIO_VERSION}-amd64.deb" \
    && gdebi -n rstudio-server.deb \
    && rm -f rstudio-server.deb



# Clean up
RUN rm -rf /var/lib/apt/lists/*

# Environment
ENV PATH=/usr/lib/rstudio-server/bin:${PATH}

# install additional lib required to run rserver
RUN apt-get update \
    && apt install -y --no-install-recommends dialog \
    && wget --no-verbose -O multiarch-support.deb http://ftp.us.debian.org/debian/pool/main/g/glibc/multiarch-support_2.19-18+deb8u10_amd64.deb \
    && gdebi -n multiarch-support.deb \
    && rm -f multiarch-support.deb \
    && wget --no-verbose -O libssl1.0.0.deb http://security.debian.org/debian-security/pool/updates/main/o/openssl/libssl1.0.0_1.0.1t-1+deb8u12_amd64.deb \
    && gdebi -n libssl1.0.0.deb \
    && rm -f libssl1.0.0.deb

EXPOSE 8787

# Launch R by default
CMD R

