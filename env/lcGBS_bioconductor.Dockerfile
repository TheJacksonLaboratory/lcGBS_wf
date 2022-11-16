FROM bioconductor/bioconductor_docker:devel
LABEL Sam Widmayer <sjwidmay@gmail.com>
RUN apt-get update \
	## Remove packages in '/var/cache/' and 'var/lib'
	## to remove side-effects of apt-get update
	&& apt-get clean \
	&& rm -rf /var/lib/apt/lists/*
RUN R -e "BiocManager::install('VariantAnnotation')"
