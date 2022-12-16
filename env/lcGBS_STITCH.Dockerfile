FROM continuumio/miniconda
LABEL Sam Widmayer <sjwidmay@gmail.com>


COPY lcGBS_STITCH.yml .
RUN \
   conda env update -n root -f lcGBS_STITCH.yml \
&& conda clean -a

RUN R -e "install.packages('parallel', repos='http://cran.us.r-project.org')"
RUN R -e "install.packages('Rcpp', repos='http://cran.us.r-project.org')"
RUN R -e "install.packages('RcppArmadillo', repos='http://cran.us.r-project.org')"
RUN R -e "library('STITCH')"