FROM ubuntu:16.04
MAINTAINER Inês Mendes <cimendes@medicina.ulisboa.pt>

RUN  apt-get update && apt-get install -y wget zip gcc

RUN apt-get install -y python-dev zlib1g-dev

WORKDIR /NGStools

# Dependencies Gallore!

# INSTALL MINICONDA
RUN wget https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh
RUN chmod +x Miniconda2-latest-Linux-x86_64.sh
RUN ./Miniconda2-latest-Linux-x86_64.sh -b -p /NGStools/miniconda
ENV PATH="/NGStools/miniconda/bin:$PATH"
# Configure channels for bioconda
RUN conda config --add channels defaults && conda config --add channels conda-forge && conda config --add channels bioconda
RUN conda install --yes -c bioconda subread bedtools bowtie2 samtools
RUN conda install --yes -c bioconda pysam

# install python packages
RUN apt-get install -y libbz2-dev libhts-dev
RUN pip install Cython
RUN pip install numpy scipy pandas scikit-learn matplotlib biopython pysam dendropy

RUN wget https://github.com/matsen/pplacer/releases/download/v1.1.alpha19/pplacer-linux-v1.1.alpha19.zip && \
    unzip pplacer-linux-v1.1.alpha19.zip && rm  pplacer-linux-v1.1.alpha19.zip
ENV PATH="/NGStools/pplacer-linux-v1.1.alpha19:$PATH"

## INSTALL CheckM

RUN apt-get install -y --no-install-recommends git && \
    git clone https://github.com/Ecogenomics/CheckM.git && \
    cd CheckM && python setup.py install

# configure checkM database
RUN wget https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz && \
    tar zxfv checkm_data_2015_01_16.tar.gz && rm checkm_data_2015_01_16.tar.gz

## INSTALL HMMER
RUN apt-get install -y make build-essential
RUN wget http://eddylab.org/software/hmmer3/3.1b2/hmmer-3.1b2.tar.gz && \
    tar -zxvf hmmer-3.1b2.tar.gz
RUN cd hmmer-3.1b2 && ./configure && make && make install
RUN cd hmmer-3.1b2/easel && make check && make install


# Install BInSanity
RUN pip install BinSanity

WORKDIR /data


