FROM rocker/r-ver:4.2
LABEL Sam Widmayer <sjwidmay@gmail.com>
install.packages('parallel', dependencies=TRUE, repos='http://cran.us.r-project.org')"
RUN Rscript -e "install.packages('data.table', dependencies=TRUE, repos='http://cran.us.r-project.org')"
RUN Rscript -e "install.packages('purrr', dependencies=TRUE, repos='http://cran.us.r-project.org')"
RUN Rscript -e "install.packages('qtl2', dependencies=TRUE, repos='http://cran.us.r-project.org')"
RUN Rscript -e "install.packages('furrr', dependencies=TRUE, repos='http://cran.us.r-project.org')"
RUN Rscript -e "install.packages('magrittr', dependencies=TRUE, repos='http://cran.us.r-project.org')"
RUN Rscript -e "install.packages('DT', dependencies=TRUE, repos='http://cran.us.r-project.org')"
RUN Rscript -e "install.packages('plotly', dependencies=TRUE, repos='http://cran.us.r-project.org')"
RUN Rscript -e "install.packages('progress', dependencies=TRUE, repos='http://cran.us.r-project.org')"
RUN Rscript -e "install.packages('ape', dependencies=TRUE, repos='http://cran.us.r-project.org')"
RUN Rscript -e "install.packages('RColorBrewer', dependencies=TRUE, repos='http://cran.us.r-project.org')"
RUN Rscript -e "install.packages('tictoc', dependencies=TRUE, repos='http://cran.us.r-project.org')"
RUN Rscript -e "install.packages('vroom', dependencies=TRUE, repos='http://cran.us.r-project.org')"
RUN Rscript -e "install.packages('fst', dependencies=TRUE, repos='http://cran.us.r-project.org')"
RUN Rscript -e "install.packages('workflowr', dependencies=TRUE, repos='http://cran.us.r-project.org')"
