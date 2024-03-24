FROM rocker/r-ver:4.2

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

RUN R -e "install.packages('parallel', repos='http://cran.us.r-project.org')"
RUN R -e "install.packages('data.table', repos='http://cran.us.r-project.org')"
RUN R -e "install.packages('vroom', repos='http://cran.us.r-project.org')"
RUN R -e "install.packages('remotes', repos='http://cran.us.r-project.org')"
RUN R -e "install.packages('Rcpp', repos='http://cran.us.r-project.org')"
RUN R -e "remotes::install_github('rqtl/qtl2')"
RUN R -e "remotes::install_github('rqtl/mmconvert')"
RUN R -e "remotes::install_github('kbroman/broman')"
RUN R -e "remotes::install_github('kbroman/qtlcharts')"
RUN R -e "install.packages('qtl2convert')"
RUN R -e "install.packages('qtl2fst',dependencies = T)"
RUN R -e "install.packages('purrr')"
RUN R -e "install.packages('furrr',dependencies = T)"
RUN R -e "install.packages('dplyr')"
RUN R -e "install.packages('tidyr')"
RUN R -e "install.packages('readr')"
RUN R -e "install.packages('magrittr')"
RUN R -e "install.packages('ggplot2')"
RUN R -e "install.packages('ggbeeswarm')"
RUN R -e "install.packages('lsa')"
RUN R -e "install.packages('cowplot')"
RUN R -e "install.packages('tictoc')"
RUN R -e "install.packages('ggrepel')"
