FROM bioconductor/bioconductor_docker:devel
LABEL Sam Widmayer <sjwidmay@gmail.com>

#RUN apt-get update \
	## Remove packages in '/var/cache/' and 'var/lib'
	## to remove side-effects of apt-get update
#	&& apt-get clean \
#	&& rm -rf /var/lib/apt/lists/*
	
RUN  apt-get --allow-releaseinfo-change update \
    && apt-get install -y procps \
    ssh \
    bash \
    pkg-config \
    libglpk-dev \
    libz-dev \
    tk \
    libxml2 \
    libxml2-dev \
    libbz2-dev \
    liblzma-dev \
    xterm \
    x11-utils \
    libcairo2-dev \
    libblas-dev \
    libssh2-1-dev \
    libgit2-dev
    
RUN R -e "BiocManager::install('VariantAnnotation')"
RUN R -e "install.packages('remotes', repos='http://cran.us.r-project.org')"
RUN R -e "remotes::install_github('rqtl/qtl2')"
RUN R -e "install.packages('qtl2convert')"
RUN R -e "remotes::install_github('rqtl/mmconvert')"
RUN R -e "install.packages('stringr')"
RUN R -e "install.packages('tidyr')"
RUN R -e "install.packages('dplyr')"
RUN R -e "install.packages('ggplot2')"
