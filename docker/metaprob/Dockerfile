FROM ubuntu:16.04
MAINTAINER I Mendes <cimendes@medicina.ulisboa.pt>

WORKDIR /NGStools/

RUN apt-get update && apt-get -y install libboost-all-dev \
    build-essential \
    wget \
    libeigen3-dev

RUN wget https://bitbucket.org/samu661/metaprob/downloads/MetaProb_v2.tar.gz
RUN tar xf MetaProb_v2.tar.gz && rm MetaProb_v2.tar.gz
RUN cd MetaProb/Release && make all

ENV PATH="NGStools/MetaProb/Release:$PATH"cd

