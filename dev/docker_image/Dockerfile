############################################################
# Dockerfile for ENCODE DCC atac-seq-pipeline
# Based on Ubuntu 18.04.3
############################################################

# Set the base image to Ubuntu 18.04.3
#FROM ubuntu:18.04
FROM ubuntu@sha256:d1d454df0f579c6be4d8161d227462d69e163a8ff9d20a847533989cf0c94d90

MAINTAINER Jin Lee

# To prevent time zone prompt
ENV DEBIAN_FRONTEND=noninteractive

# Install softwares from apt repo
RUN apt-get update && apt-get install -y \
    libncurses5-dev libcurl4-openssl-dev libfreetype6-dev zlib1g-dev \
    python python-setuptools python-pip python3 python3-setuptools python3-pip \
    git wget unzip ghostscript pkg-config libboost-dev \
    default-jre apt-transport-https tabix \
    r-base \
    && rm -rf /var/lib/apt/lists/*

# Make directory for all softwares
RUN mkdir /software
WORKDIR /software
ENV PATH="/software:${PATH}"

# Downgrade openssl to 1.0.2t (to get the same random number for SPR)
RUN wget https://github.com/openssl/openssl/archive/OpenSSL_1_0_2t.tar.gz && tar zxvf OpenSSL_1_0_2t.tar.gz && cd openssl-OpenSSL_1_0_2t/ && ./config && make && make install && cd ../ && rm -rf openssl-OpenSSL_1_0_2t* && rm /usr/bin/openssl && ln -s /usr/local/ssl/bin/openssl /usr/bin/openssl

# Install system/math python packages (python3)
RUN pip3 install --no-cache-dir jsondiff==1.1.1 common python-dateutil cython pandas==0.25.1 jinja2==2.10.1 matplotlib==3.1.1

# Install genomic python package (python3)
RUN pip3 install --no-cache-dir pyBigwig==0.3.13 cutadapt==2.5 pyfaidx==0.5.5.2 pybedtools==0.8.0 pysam==0.15.3 deeptools==3.3.1

# Install R packages including spp
RUN echo "r <- getOption('repos'); r['CRAN'] <- 'http://cran.r-project.org'; options(repos = r);" > ~/.Rprofile && \
    Rscript -e "install.packages('snow')" && \
    Rscript -e "install.packages('snowfall')" && \
    Rscript -e "install.packages('bitops')" && \
    Rscript -e "install.packages('caTools')" && \
    Rscript -e "install.packages('Rcpp')"

# Install bioconductor and Rsamtools which is required by spp package
RUN Rscript -e "source('http://bioconductor.org/biocLite.R'); biocLite('Rsamtools')"

# Install r-spp 1.15 (1.13 in Conda env)
RUN wget https://github.com/hms-dbmi/spp/archive/1.15.tar.gz && Rscript -e "install.packages('./1.15.tar.gz')" && rm -f https://github.com/hms-dbmi/spp/archive/1.15.tar.gz

# Install samtools 1.9
RUN git clone --branch 1.9 --single-branch https://github.com/samtools/samtools.git && \
    git clone --branch 1.9 --single-branch https://github.com/samtools/htslib.git && \
    cd samtools && make && make install && cd ../ && rm -rf samtools* htslib*

# Install bedtools 2.29.0
RUN git clone --branch v2.29.0 --single-branch https://github.com/arq5x/bedtools2.git && \
    cd bedtools2 && make && make install && cd ../ && rm -rf bedtools2*

# Install Picard 2.20.7
RUN wget https://github.com/broadinstitute/picard/releases/download/2.20.7/picard.jar && chmod +x picard.jar

# Install sambamba 0.6.6
RUN wget https://github.com/lomereiter/sambamba/releases/download/v0.6.6/sambamba_v0.6.6_linux.tar.bz2 && tar -xvjf sambamba_v0.6.6_linux.tar.bz2 && mv sambamba_v0.6.6 sambamba && rm -rf sambamba_*

# Install gsl 1.16
RUN wget http://gnu.mirror.vexxhost.com/gsl/gsl-1.16.tar.gz && tar -zxvf gsl-1.16.tar.gz && cd gsl-1.16 && ./configure && make && make install && cd .. && rm -rf gsl-1.16 gsl-1.16.tar.gz
ENV LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/usr/local/lib

# Install preseq 2.0.3
RUN git clone --branch v2.0.3 --single-branch --recursive https://github.com/smithlabcode/preseq preseq_2.0.3 && cd preseq_2.0.3 && make && cd ../ && mv preseq_2.0.3/preseq . && rm -rf preseq_2.0.3

# Install Bowtie2 2.3.4.3
RUN wget https://github.com/BenLangmead/bowtie2/releases/download/v2.3.4.3/bowtie2-2.3.4.3-linux-x86_64.zip && \
    unzip bowtie2-2.3.4.3-linux-x86_64.zip && mv bowtie2*/bowtie2* . && rm -rf bowtie2-2.3.4.3*

# Install Bwa 0.7.17
RUN git clone --branch v0.7.17 --single-branch https://github.com/lh3/bwa.git && \
    cd bwa && make && cp bwa /usr/local/bin/ && cd ../ && rm -rf bwa*

# Install phantompeakqualtools 1.2.1
RUN wget https://github.com/kundajelab/phantompeakqualtools/archive/1.2.1.tar.gz && tar -xvf 1.2.1.tar.gz && rm -f 1.2.1.tar.gz
ENV PATH="/software/phantompeakqualtools-1.2.1:${PATH}"

# Install SAMstats
RUN pip3 install --no-cache-dir SAMstats==0.2.1

# Install IDR 2.0.4.2
RUN git clone --branch 2.0.4.2 --single-branch https://github.com/kundajelab/idr && \
    cd idr && python3 setup.py install && cd ../ && rm -rf idr*

# Install system/math python packages (python2)
RUN pip2 install --no-cache-dir numpy matplotlib==2.2.4

# Install genomic python packages (python2)
RUN pip2 install --no-cache-dir Cython
RUN pip2 install --no-cache-dir macs2==2.1.3.3 metaseq==0.5.6

# Install UCSC tools (v377)
RUN git clone https://github.com/ENCODE-DCC/kentUtils_bin_v377
ENV PATH=${PATH}:/software/kentUtils_bin_v377/bin
ENV LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/software/kentUtils_bin_v377/lib

# Prevent conflict with locally installed python outside of singularity container
ENV PYTHONNOUSERSITE=True

# Disable multithreading for BLAS
ENV OPENBLAS_NUM_THREADS=1
ENV MKL_NUM_THREADS=1

# Get ENCODE atac-seq-pipeline container repository
# This COPY assumes the build context is the root of the atac-seq-pipeline repo
# and it gets whatever is checked out plus local modifications
# so the buildling command should be:
# cd [GIT_REPO_DIR] && docker build -f dev/docker_image/Dockerfile .
RUN mkdir -p atac-seq-pipeline/src
COPY src atac-seq-pipeline/src/
COPY atac.wdl atac-seq-pipeline/
ENV PATH="/software/atac-seq-pipeline:/software/atac-seq-pipeline/src:${PATH}"

# make some temporary directories
RUN mkdir -p /mnt/ext_{0..9}
