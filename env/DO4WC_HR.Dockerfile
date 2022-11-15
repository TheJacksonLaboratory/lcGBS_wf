FROM ubuntu:latest
LABEL Sam Widmayer <sjwidmay@gmail.com>
RUN apt-get update
RUN apt-get install -y \
    libxml2 \
    libxt6 \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libpcre3-dev \
    libicu-dev \
    libjpeg-dev \
    libpng-dev \
    libxml2-dev \
    libglpk-dev

FROM rocker/r-ver:4.2
RUN Rscript -e "install.packages('parallel', dependencies=TRUE, repos='http://cran.us.r-project.org')"
RUN Rscript -e "install.packages('data.table', dependencies=TRUE, repos='http://cran.us.r-project.org')"
RUN Rscript -e "install.packages('vroom', dependencies=TRUE, repos='http://cran.us.r-project.org')"
# Installing qtl2 utils
RUN Rscript -e "install.packages('remotes', dependencies=TRUE, repos='http://cran.us.r-project.org')"
RUN Rscript -e "install.packages('Rcpp', dependencies=TRUE, repos='http://cran.us.r-project.org')"
RUN Rscript -e "remotes::install_github('rqtl/qtl2')"
RUN Rscript -e "remotes::install_github('rqtl/qtl2convert')"
RUN Rscript -e "remotes::install_github('rqtl/qtl2fst')"
RUN Rscript -e "remotes::install_github('byandell/qtl2ggplot')"
# Installing bioconductor utils
RUN Rscript -e "install.packages('BiocManager')"
RUN Rscript -e "BiocManager::install('VariantAnnotation)"
