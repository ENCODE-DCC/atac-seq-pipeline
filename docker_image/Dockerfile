############################################################
# Dockerfile for ENCODE DCC atac-seq-pipeline
# Based on Ubuntu 16.04
############################################################

# IMPORTANT!
# If you install python2/3 packages using pip/pip3
#  and not sure about math library dependencies like BLAS and numpy,
#  then install with --no-dependencies

# Set the base image to Ubuntu 16.04
FROM ubuntu:16.04

# File Author / Maintainer
MAINTAINER Jin Lee

# Update the repository sources list
# Install base packages: git, python2/3, java, R
RUN apt-get update && apt-get install -y \
    ant \
    build-essential \
    software-properties-common \
    libncurses5-dev \
    libncursesw5-dev \
    libcurl4-openssl-dev \
    libboost-dev \
    libfreetype6-dev \
    zlib1g-dev \
    r-base \
    default-jdk \
    python-dev \
    python-setuptools \
    python-pip \
    python3-dev \
    python3-setuptools \
    python3-pip \
    git \
    wget \
    unzip \
    ghostscript \
    pkg-config

# Make directory for all softwares
RUN mkdir /software
WORKDIR /software
ENV PATH="/software:${PATH}"

# Install OpenBLAS 0.2.19 (NEHALEM architecture, single-threaded)
# to have the same result as of Conda's non-MKL openblas 0.2.19
# setting TARGET=HASWELL will give you different result.
# so let's stick to TARGET=NEHALEM
RUN git clone --branch v0.2.19 --single-branch https://github.com/xianyi/OpenBLAS
RUN cd OpenBLAS && make FC=gfortran TARGET=NEHALEM USE_THREAD=0 && make PREFIX=/opt/openblas install
ENV LD_LIBRARY_PATH="/opt/openblas/lib:${LD_LIBRARY_PATH}"

# Install samtools dependencies
RUN wget http://zlib.net/zlib-1.2.11.tar.gz && tar -xvf zlib-1.2.11.tar.gz
RUN cd zlib-1.2.11 && ./configure && make && make install
RUN wget http://bzip.org/1.0.6/bzip2-1.0.6.tar.gz && tar -xvf bzip2-1.0.6.tar.gz
RUN cd bzip2-1.0.6 && make && make install
RUN wget https://tukaani.org/xz/xz-5.2.3.tar.gz && tar -xvf xz-5.2.3.tar.gz
RUN cd xz-5.2.3 && ./configure && make && make install

# Install basic python2 packages
RUN pip install common python-dateutil cython

# Install basic python3 packages
RUN pip3 install common python-dateutil cython

# Install numpy 1.11.3 (python2/3)
RUN git clone --branch v1.11.3 --single-branch https://github.com/numpy/numpy
COPY /docker_image/site.cfg numpy/
RUN cd numpy && python setup.py install
RUN cd numpy && python3 setup.py install

# Install scipy 1.0.0 (python2/3)
RUN git clone --branch v1.0.0 --single-branch https://github.com/scipy/scipy
RUN cd scipy && python setup.py install
RUN cd scipy && python3 setup.py install

# Install matplotlib 1.5.1 (python2/3)
RUN git clone --branch v1.5.1 --single-branch https://github.com/matplotlib/matplotlib
RUN cd matplotlib && python setup.py install
RUN cd matplotlib && python3 setup.py install

# Install MACS2 2.1.1.20160309 (python2)
RUN pip install --no-dependencies macs2==2.1.1.20160309

# Install IDR 2.0.4.1 (python3)
RUN git clone --branch 2.0.4.1 --single-branch https://github.com/kundajelab/idr
RUN cd idr && python3 setup.py install

# Install UCSC tools (Kent utils) latest
RUN wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bedClip
RUN chmod +x bedClip
RUN wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bedGraphToBigWig
RUN chmod +x bedGraphToBigWig

# Install samtools 1.2
RUN git clone --branch 1.2 --single-branch https://github.com/samtools/samtools.git
RUN git clone --branch 1.2 --single-branch git://github.com/samtools/htslib.git
RUN cd samtools && make && make install

# Install bedtools 2.26.0
RUN git clone --branch v2.26.0 --single-branch https://github.com/arq5x/bedtools2.git
RUN cd bedtools2 && make && make install

# Install Picard 2.10.6
RUN wget https://github.com/broadinstitute/picard/releases/download/2.10.6/picard.jar
RUN chmod +x picard.jar

# Install sambamba 0.6.5
RUN wget https://github.com/lomereiter/sambamba/releases/download/v0.6.5/sambamba_v0.6.5_linux.tar.bz2 && tar -xvjf sambamba_v0.6.5_linux.tar.bz2
RUN mv sambamba_v0.6.5 sambamba

# Install R packages
RUN echo "r <- getOption('repos'); r['CRAN'] <- 'http://cran.us.r-project.org'; options(repos = r);" > ~/.Rprofile
RUN Rscript -e "install.packages('snow')"
RUN Rscript -e "install.packages('snowfall')"
RUN Rscript -e "install.packages('bitops')"
RUN Rscript -e "install.packages('caTools')"
RUN Rscript -e "source('http://bioconductor.org/biocLite.R'); biocLite('Rsamtools')"

# Install R package spp 1.13 (required for phantompeakqualtools)
RUN wget https://github.com/hms-dbmi/spp/archive/1.13.tar.gz
RUN Rscript -e "install.packages('./1.13.tar.gz')"

# Install phantompeakqualtools 1.2
RUN wget https://github.com/kundajelab/phantompeakqualtools/archive/1.2.tar.gz && tar -xvf 1.2.tar.gz
ENV PATH="/software/phantompeakqualtools-1.2:${PATH}"

# Install cutadapt 1.9.1
RUN pip install --no-dependencies cutadapt==1.9.1 

# Install Bowtie2 2.2.6
RUN wget https://github.com/BenLangmead/bowtie2/releases/download/v2.2.6/bowtie2-2.2.6-linux-x86_64.zip
RUN unzip bowtie2-2.2.6-linux-x86_64.zip
ENV PATH="/software/bowtie2-2.2.6:${PATH}"

# Get ENCODE atac-seq-pipeline container repository
# This COPY assumes the build context is the root of the atac-seq-pipeline repo
# and it gets whatever is checked out plus local modifications
# so the buildling command should:
# cd [GIT_REPO_DIR] && docker build -f docker_images/pipeline/Dockerfile .
RUN mkdir -p atac-seq-pipeline/src
COPY /src atac-seq-pipeline/src
ENV PATH="/software/atac-seq-pipeline/src:${PATH}"

ENTRYPOINT ["/bin/bash","-c"]
