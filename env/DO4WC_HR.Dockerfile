# Docker inheritance for bioconductor
LABEL Sam Widmayer <sjwidmay@gmail.com>

# Installing other R packages
FROM rocker/r-ver:4.2
RUN R -e "install.packages('parallel', repos='http://cran.us.r-project.org')"
RUN R -e "install.packages('data.table', repos='http://cran.us.r-project.org')"
RUN R -e "install.packages('vroom', repos='http://cran.us.r-project.org')"
RUN R -e "install.packages('remotes', repos='http://cran.us.r-project.org')"
RUN R -e "install.packages('Rcpp', repos='http://cran.us.r-project.org')"
RUN R -e "remotes::install_github('rqtl/qtl2')"
RUN R -e "remotes::install_github('rqtl/qtl2convert')"
RUN R -e "remotes::install_github('rqtl/qtl2fst')"
RUN R -e "remotes::install_github('byandell/qtl2ggplot',force = T)"
RUN R -e "install.packages('purrr')"
RUN R -e "install.packages('dplyr')"
RUN R -e "install.packages('magrittr')"
RUN R -e "install.packages('ggplot2')"
RUN R -e "install.packages('tidyr')"

FROM bioconductor/bioconductor_docker:devel
RUN apt-get update \
	## Remove packages in '/var/cache/' and 'var/lib'
	## to remove side-effects of apt-get update
	&& apt-get clean \
	&& rm -rf /var/lib/apt/lists/*
RUN R -e "BiocManager::install('VariantAnnotation')"
