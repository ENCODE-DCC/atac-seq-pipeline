############################################################
# Dockerfile for ENCODE DCC containerized pipelines
# Based on Ubuntu 14.04
############################################################

# Set the base image to Ubuntu 14.04
FROM ubuntu:14.04

# File Author / Maintainer
MAINTAINER Jin Lee

# Update the repository sources list
# Install base packages: git, python, java
RUN apt-get update && apt-get install -y \
    ant \
    build-essential \
    git \
    libncurses5-dev \
    libncursesw5-dev \
    python-dev \
    python-setuptools \
    software-properties-common \
    wget \
    zlib1g-dev

# Install pip and python packages: common and python-dateutils
RUN easy_install pip
RUN pip install common
RUN pip install python-dateutil

# Install Java 8 for picard
RUN echo debconf shared/accepted-oracle-license-v1-1 select true | \
    sudo debconf-set-selections
RUN echo debconf shared/accepted-oracle-license-v1-1 seen true | \
    sudo debconf-set-selections
RUN add-apt-repository ppa:webupd8team/java
RUN apt-get update && apt-get install -y \
	oracle-java8-installer

RUN mkdir /image_software
WORKDIR /image_software

# Install bedtools 2.26.0
RUN git clone --branch v2.26.0 https://github.com/arq5x/bedtools2.git
RUN cd bedtools2 && make && make install

# Install samtools dependencies
RUN wget http://zlib.net/zlib-1.2.11.tar.gz && tar -xvf zlib-1.2.11.tar.gz
RUN cd zlib-1.2.11 && ./configure && make && make install
RUN wget http://bzip.org/1.0.6/bzip2-1.0.6.tar.gz && tar -xvf bzip2-1.0.6.tar.gz
RUN cd bzip2-1.0.6 && make && make install
RUN wget https://tukaani.org/xz/xz-5.2.3.tar.gz && tar -xvf xz-5.2.3.tar.gz
RUN cd xz-5.2.3 && ./configure && make && make install
# Install samtools 1.5
RUN git clone --branch 1.5 https://github.com/samtools/samtools.git
RUN git clone --branch 1.5 git://github.com/samtools/htslib.git
RUN cd samtools && make && make install

# Install Picard 2.10.2
RUN wget https://github.com/broadinstitute/picard/releases/download/2.10.2/picard.jar
ENV PICARD_HOME /image_software/

# Get ENCODE pipeline container repository
# This COPY asumes the build context is the root of the pipeline-container repo
# and it gets whatever is checked out plus local modifications
RUN mkdir pipeline-container
COPY / pipeline-container

# Set up the user directory
RUN groupadd -r encode && useradd --no-log-init -r -m -d /home/encode/ -g encode encode
USER encode
WORKDIR /home/encode

#ENTRYPOINT ["python", "/image_software/pipeline-container/src/filter_qc.py"]
ENTRYPOINT ["/bin/bash","-c"]
