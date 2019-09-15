FROM ubuntu:16.04
MAINTAINER Inês Mendes <cimendes@medicina.ulisboa.pt>

RUN  apt-get update && apt-get install -y wget curl unzip

WORKDIR /NGStools

# INSTALL MAXBIN
RUN apt-get install -y g++ perl build-essential autoconf
RUN perl -MCPAN -e 'install LWP::Simple'
RUN wget https://sourceforge.net/projects/maxbin/files/MaxBin-2.2.7.tar.gz
RUN tar xf MaxBin-2.2.7.tar.gz && rm MaxBin-2.2.7.tar.gz
RUN cd MaxBin-2.2.7/src && make
RUN cd MaxBin-2.2.7 && ./autobuild_auxiliary

ENV PATH="/NGStools/MaxBin-2.2.7:$PATH"