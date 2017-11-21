############################################################
# Dockerfile for ENCODE DCC atac-seq-pipeline
# Based on Ubuntu 14.04
############################################################

# Set the base image to Ubuntu 14.04
FROM ubuntu:14.04

# File Author / Maintainer
MAINTAINER Jin Lee

# Update the repository sources list
# Install base packages: git, python, java, R
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
    zlib1g-dev \
    unzip \
    libcurl4-openssl-dev \
    r-base

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

RUN mkdir /software
WORKDIR /software
ENV PATH="/software:${PATH}"

# Install cutadapt 1.9.1
RUN pip install cutadapt==1.9.1

# Install samtools dependencies
RUN wget http://zlib.net/zlib-1.2.11.tar.gz && tar -xvf zlib-1.2.11.tar.gz
RUN cd zlib-1.2.11 && ./configure && make && make install
RUN wget http://bzip.org/1.0.6/bzip2-1.0.6.tar.gz && tar -xvf bzip2-1.0.6.tar.gz
RUN cd bzip2-1.0.6 && make && make install
RUN wget https://tukaani.org/xz/xz-5.2.3.tar.gz && tar -xvf xz-5.2.3.tar.gz
RUN cd xz-5.2.3 && ./configure && make && make install

# Install samtools 1.2
RUN git clone --branch 1.2 https://github.com/samtools/samtools.git
RUN git clone --branch 1.2 git://github.com/samtools/htslib.git
RUN cd samtools && make && make install

# Install bedtools 2.26.0
RUN git clone --branch v2.26.0 https://github.com/arq5x/bedtools2.git
RUN cd bedtools2 && make && make install

# Install Picard 2.10.6
RUN wget https://github.com/broadinstitute/picard/releases/download/2.10.6/picard.jar
RUN chmod +x picard.jar

# Install Bowtie2 2.2.6
RUN wget https://github.com/BenLangmead/bowtie2/releases/download/v2.2.6/bowtie2-2.2.6-linux-x86_64.zip
RUN unzip bowtie2-2.2.6-linux-x86_64.zip
ENV PATH="/software/bowtie2-2.2.6:${PATH}"

# Install UCSC tools (Kent utils) latest
RUN wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bedClip
RUN chmod +x bedClip
RUN wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bedGraphToBigWig
RUN chmod +x bedGraphToBigWig

# Install NumPy 1.10.2 (python2)
RUN pip install numpy==1.10.2

# Install MACS2 2.1.0
RUN pip install macs2==2.1.0.20151222

# Install R packages
RUN echo "r <- getOption('repos'); r['CRAN'] <- 'http://cran.us.r-project.org'; options(repos = r);" > ~/.Rprofile
RUN Rscript -e "install.packages('snow')"
RUN Rscript -e "install.packages('snowfall')"
RUN Rscript -e "install.packages('bitops')"
RUN Rscript -e "install.packages('caTools')"
RUN Rscript -e "source('http://bioconductor.org/biocLite.R'); biocLite('Rsamtools')"

# Install boost (required for R-spp)
RUN apt-get install -y libboost-dev

# Install R package spp 1.14 (required for phantompeakqualtools)
#RUN wget https://github.com/hms-dbmi/spp/archive/1.14.tar.gz
#RUN Rscript -e "install.packages('./1.14.tar.gz')"

# Install R package spp 1.13 (required for phantompeakqualtools)
RUN wget https://github.com/hms-dbmi/spp/archive/1.13.tar.gz
RUN Rscript -e "install.packages('./1.13.tar.gz')"

# Install phantompeakqualtools 1.2
RUN wget https://github.com/kundajelab/phantompeakqualtools/archive/1.2.tar.gz && tar -xvf 1.2.tar.gz
ENV PATH="/software/phantompeakqualtools-1.2:${PATH}"

# Install Python3 and packages (for IDR)
RUN apt-get install -y python3-dev python3-numpy python3-scipy python3-setuptools python3-matplotlib

# Install IDR 2.0.4
RUN git clone --branch 2.0.4 https://github.com/kundajelab/idr
RUN cd idr && python3 setup.py install

# Install sambamba 0.6.5
RUN wget https://github.com/lomereiter/sambamba/releases/download/v0.6.5/sambamba_v0.6.5_linux.tar.bz2 && tar -xvjf sambamba_v0.6.5_linux.tar.bz2
RUN mv sambamba_v0.6.5 sambamba

# Install ghostscript for reporting module (pdf2png)
RUN apt-get install -y ghostscript

# Get ENCODE atac-seq-pipeline container repository
# This COPY assumes the build context is the root of the atac-seq-pipeline repo
# and it gets whatever is checked out plus local modifications
RUN mkdir atac-seq-pipeline
COPY / atac-seq-pipeline
ENV PATH="/software/atac-seq-pipeline/src:${PATH}"

ENTRYPOINT ["/bin/bash","-c"]
