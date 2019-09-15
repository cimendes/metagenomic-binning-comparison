FROM ubuntu:16.04
MAINTAINER I Mendes <cimendes@medicina.ulisboa.pt>

WORKDIR /NGStools/

RUN apt-get update && apt-get -y install perl r-base git wget make build-essential

# Samtools 0.1.19
RUN apt-get install -qq samtools

# Bowtie2.1.0
RUN apt-get install -qq bowtie2

# R packages
RUN Rscript -e "install.packages(c('ggplot2','Rtsne','getopt','dbscan', 'RColorBrewer'), repos = 'http://cran.us.r-project.org')"

# Install Esel
RUN git clone https://github.com/EddyRivasLab/easel && cd easel && autoconf
RUN cd easel && ./configure && make

RUN git clone https://github.com/laxeye/YAMB.git
ENV PATH="/NGStools/YAMB:$PATH"