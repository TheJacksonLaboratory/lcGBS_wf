# Docker inheritance for bioconductor
FROM bioconductor/bioconductor_docker:devel
LABEL Sam Widmayer <sjwidmay@gmail.com>

RUN apt-get update \
	## Remove packages in '/var/cache/' and 'var/lib'
	## to remove side-effects of apt-get update
	&& apt-get clean \
	&& rm -rf /var/lib/apt/lists/*
RUN Rscript -e "BiocManager::install('VariantAnnotation')"

# Installing other R packages
FROM rocker/r-ver:4.2
RUN Rscript -e "install.packages('parallel', repos='http://cran.us.r-project.org')"
RUN Rscript -e "install.packages('data.table', repos='http://cran.us.r-project.org')"
RUN Rscript -e "install.packages('vroom', repos='http://cran.us.r-project.org')"
RUN Rscript -e "install.packages('remotes', repos='http://cran.us.r-project.org')"
RUN Rscript -e "install.packages('Rcpp', repos='http://cran.us.r-project.org')"
RUN Rscript -e "remotes::install_github('rqtl/qtl2')"
RUN Rscript -e "remotes::install_github('rqtl/qtl2convert')"
RUN Rscript -e "remotes::install_github('rqtl/qtl2fst')"
RUN Rscript -e "remotes::install_github('byandell/qtl2ggplot',force = T)"
RUN Rscript -e "install.packages('purrr')"
RUN Rscript -e "install.packages('dplyr')"
RUN Rscript -e "install.packages('magrittr')"
RUN Rscript -e "install.packages('ggplot2')"
RUN Rscript -e "install.packages('tidyr')"
