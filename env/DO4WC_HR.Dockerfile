FROM rocker/r-ver:4.2
LABEL Sam Widmayer <sjwidmay@gmail.com>
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
RUN Rscript -e "remotes::install_github('Bioconductor/VariantAnnotation', dependencies=TRUE)"