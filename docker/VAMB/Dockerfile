FROM ubuntu:16.04
MAINTAINER I Mendes <cimendes@medicina.ulisboa.pt>

WORKDIR /NGStools/

RUN apt-get update && apt-get -y install wget unzip make gcc zlib1g-dev libbz2-dev

# INSTALL MINICONDA
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
RUN chmod +x Miniconda3-latest-Linux-x86_64.sh
RUN ./Miniconda3-latest-Linux-x86_64.sh -b -p /NGStools/miniconda
ENV PATH="/NGStools/miniconda/bin:$PATH"
# Configure channels for bioconda
RUN conda config --add channels defaults && conda config --add channels conda-forge && conda config --add channels bioconda

RUN wget https://github.com/RasmussenLab/vamb/archive/v2.0.1.zip && unzip v2.0.1.zip && rm v2.0.1.zip
RUN cd vamb* && pip install -e .
