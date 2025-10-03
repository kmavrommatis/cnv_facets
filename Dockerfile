FROM ubuntu:24.04

RUN apt-get update 
RUN DEBIAN_FRONTEND=noninteractive \
    apt-get install -y \
        git \
        wget \
        python3 \
        r-base \
        libssl-dev \
        libcurl4-openssl-dev \
        libxml2-dev \
        libfontconfig1-dev \
        libharfbuzz-dev \
        libfribidi-dev \
        libfreetype6-dev \
        libpng-dev \
        libtiff5-dev \
        libjpeg-dev \
        rsync \
        curl

RUN R -e 'install.packages(c("devtools","BiocManager"),dependencies=TRUE); '
RUN R -e 'BiocManager::install("Rsamtools")'
RUN R -e 'devtools::install_github("mskcc/facets", ref = "f3c93ee")'

RUN R -e 'install.packages(c("testthat","covr","data.table","gridExtra", "ggplot2"),dependencies=TRUE);'
RUN R -e 'devtools::install_github(c("trevorld/argparse"))'

RUN apt-get install -y samtools. 
RUN apt-get install -y bcftools

WORKDIR /usr/local/lib/R/site-library/facets/extcode/
RUN ln -s /usr/local/lib/R/site-library/Rhtslib/usrlib/libhts.a /usr/local/lib/R/site-library/Rhtslib/usrlib/libhts-static.a && \
    g++ -std=c++11 -I /usr/local/lib/R/site-library/Rhtslib/include/ snp-pileup.cpp -L /usr/local/lib/R/site-library/Rhtslib/usrlib/ -lhts-static -o snp-pileup -lcurl -lz -lpthread -lcrypto -llzma -lbz2 &&\
    ln snp-pileup /usr/local/bin/snp-pileup

RUN R -e 'install.packages(c("future","future.apply"))'
RUN R -e "install.packages('box', repos = 'https://klmr.r-universe.dev')"
RUN R -e "install.packages('memuse')"

WORKDIR /opt/cnv_facets/
COPY . /opt/cnv_facets/
RUN chmod uog+x /opt/cnv_facets/bin/cnv_facets.R
# we have installed all the packages in the previous commands - no need to run the setup script.
# RUN /opt/cnv_facets/setup.sh --bin_dir /usr/bin


COPY ./Dockerfile /opt/Dockerfile

