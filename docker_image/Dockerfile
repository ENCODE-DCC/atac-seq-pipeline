############################################################
# Dockerfile for ENCODE DCC atac-seq-pipeline
# Based on Ubuntu 16.04
############################################################

# IMPORTANT!
# If you install python2/3 packages using pip/pip3
#  and not sure about math library dependencies like BLAS and numpy,
#  then install with --no-dependencies

# Set the base image to Ubuntu 16.04
#FROM ubuntu:16.04
FROM ubuntu@sha256:e4a134999bea4abb4a27bc437e6118fdddfb172e1b9d683129b74d254af51675

# File Author / Maintainer
MAINTAINER Jin Lee

# Update the repository sources list
# Install base packages: git, python2/3, java, R
RUN apt-get update && apt-get install -y \
    libncurses5-dev \
    libncursesw5-dev \
    libcurl4-openssl-dev \
    libfreetype6-dev \
    zlib1g-dev \
    python \
    python-setuptools \
    python-pip \
    python3 \
    python3-setuptools \
    python3-pip \
    git \
    wget \
    unzip \
    ghostscript \
    pkg-config \
    libboost-dev \
    r-base-core \
    default-jre \
    apt-transport-https \
    tabix \
&& rm -rf /var/lib/apt/lists/*

# Install Intel MKL for BLAS
RUN wget https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS-2019.PUB && apt-key add GPG-PUB-KEY-INTEL-SW-PRODUCTS-2019.PUB && rm -rf GPG-PUB-KEY-INTEL-SW-PRODUCTS-2019.PUB && sh -c 'echo deb https://apt.repos.intel.com/mkl all main > /etc/apt/sources.list.d/intel-mkl.list' && apt-get update && apt-get install intel-mkl-64bit-2018.0-033 -y && rm -rf /var/lib/apt/lists/*
ENV LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/opt/intel/compilers_and_libraries_2018.0.128/linux/mkl/lib/intel64_lin

# Install basic python2/3 packages
RUN pip install --no-cache-dir common python-dateutil cython==0.27.3 && \
    pip3 install --no-cache-dir common python-dateutil cython==0.27.3

# Make directory for all softwares
RUN mkdir /software
WORKDIR /software
ENV PATH="/software:${PATH}"

# Install numpy 1.11.3 (python2/3, linked with MKL)
RUN git clone --branch v1.11.3 https://github.com/numpy/numpy && cd numpy && \
    /bin/bash -c 'echo -e "[mkl]\nlibrary_dirs = /opt/intel/compilers_and_libraries_2018/linux/mkl/lib/intel64\ninclude_dirs = /opt/intel/compilers_and_libraries_2018/linux/mkl/include\nmkl_libs = mkl_rt\nlapack_libs =" > site.cfg' && python setup.py install && python3 setup.py install && cd ../ && rm -rf numpy*

# Install scipy 1.0.0 (python2/3)
RUN git clone --branch v1.0.0 --single-branch https://github.com/scipy/scipy && \
    cd scipy && python setup.py install && python3 setup.py install && cd ../ && rm -rf scipy*

# Install matplotlib 1.5.1 (python2/3)
RUN git clone --branch v1.5.1 --single-branch https://github.com/matplotlib/matplotlib && \
    cd matplotlib && python setup.py install && python3 setup.py install && cd ../ && rm -rf matplotlib*

# Install MACS2 2.1.1.20160309 (python2)
RUN pip install --no-cache-dir --no-dependencies macs2==2.1.1.20160309

# Install IDR 2.0.4.2 (python3)
RUN git clone --branch 2.0.4.2 --single-branch https://github.com/kundajelab/idr && \
    cd idr && python3 setup.py install && cd ../ && rm -rf idr*

# Install UCSC tools (Kent utils) latest
RUN wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bedClip && chmod +x bedClip
RUN wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bedGraphToBigWig && chmod +x bedGraphToBigWig

# Install samtools 1.2
RUN git clone --branch 1.2 --single-branch https://github.com/samtools/samtools.git && \
    git clone --branch 1.2 --single-branch https://github.com/samtools/htslib.git && \
    cd samtools && make && make install && cd ../ && rm -rf samtools* htslib*

# Install bedtools 2.26.0
RUN git clone --branch v2.26.0 --single-branch https://github.com/arq5x/bedtools2.git && \
    cd bedtools2 && make && make install && cd ../ && rm -rf bedtools2*

# Install Picard 2.10.6
RUN wget https://github.com/broadinstitute/picard/releases/download/2.10.6/picard.jar && chmod +x picard.jar

# Install sambamba 0.6.6
RUN wget https://github.com/lomereiter/sambamba/releases/download/v0.6.6/sambamba_v0.6.6_linux.tar.bz2 && tar -xvjf sambamba_v0.6.6_linux.tar.bz2 && mv sambamba_v0.6.6 sambamba && rm -rf sambamba_*

# Install R packages
RUN echo "r <- getOption('repos'); r['CRAN'] <- 'http://cran.r-project.org'; options(repos = r);" > ~/.Rprofile && \
    Rscript -e "install.packages('snow')" && \
    Rscript -e "install.packages('snowfall')" && \
    Rscript -e "install.packages('bitops')" && \
    Rscript -e "install.packages('caTools')" && \
    Rscript -e "source('http://bioconductor.org/biocLite.R'); biocLite('Rsamtools')"

# Install R package spp 1.13 (required for phantompeakqualtools)
RUN wget https://github.com/hms-dbmi/spp/archive/1.13.tar.gz && Rscript -e "install.packages('./1.13.tar.gz')" && rm -f 1.13.tar.gz

# Install phantompeakqualtools 1.2
RUN wget https://github.com/kundajelab/phantompeakqualtools/archive/1.2.tar.gz && tar -xvf 1.2.tar.gz && rm -f 1.2.tar.gz
ENV PATH="/software/phantompeakqualtools-1.2:${PATH}"

# Install cutadapt 1.9.1
RUN pip install --no-cache-dir --no-dependencies cutadapt==1.9.1 

# Install Bowtie2 2.2.6
RUN wget https://github.com/BenLangmead/bowtie2/releases/download/v2.2.6/bowtie2-2.2.6-linux-x86_64.zip && \
    unzip bowtie2-2.2.6-linux-x86_64.zip && mv bowtie2*/bowtie2* . && rm -rf bowtie2-2.2.6*

# Install pyfaidx (for building genome data)
RUN pip install --no-cache-dir pyfaidx==0.4.7.1

# Install some apt packages (liblzo2, libgsl2) for ATAQC
RUN apt-get update && apt-get install -y \
    liblzo2-dev \
&& rm -rf /var/lib/apt/lists/*

# Install python packages for ATAQC (pysam, pybedtools, metaseq, pandas, jinja2)
RUN pip install --no-cache-dir pysam==0.8.2.1 pybedtools==0.6.9 pandas==0.21.1 metaseq==0.5.6 jinja2==2.10

# Install gsl 1.16
RUN wget http://gnu.mirror.vexxhost.com/gsl/gsl-1.16.tar.gz && tar -zxvf gsl-1.16.tar.gz && cd gsl-1.16 && ./configure && make && make install && cd .. && rm -rf gsl-1.16 gsl-1.16.tar.gz
ENV LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/usr/local/lib

# Install preseq 2.0.3 for ATAQC (build from source)
RUN git clone --branch v2.0.3 --single-branch --recursive https://github.com/smithlabcode/preseq preseq_2.0.3 && cd preseq_2.0.3 && make && cd ../ && mv preseq_2.0.3/preseq . && rm -rf preseq_2.0.3

# Install UCSC tools (Kent utils) for ATAQC
RUN wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bigWigAverageOverBed && chmod +x bigWigAverageOverBed

# Install UCSC tools (Kent utils) latest for peak2bigbed conversion
RUN wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bedToBigBed && chmod +x bedToBigBed

# Install bgzip/tabix for Wash U browser track (hammock type)
#RUN apt-get update && apt-get install -y tabix && rm -rf /var/lib/apt/lists/*

# Install python jsondiff
RUN pip install --no-cache-dir jsondiff==1.1.1

# Prevent conflict with locally installed python outside of singularity container
ENV PYTHONNOUSERSITE=True

# Get ENCODE atac-seq-pipeline container repository
# This COPY assumes the build context is the root of the atac-seq-pipeline repo
# and it gets whatever is checked out plus local modifications
# so the buildling command should:
# cd [GIT_REPO_DIR] && docker build -f docker_images/pipeline/Dockerfile .
RUN mkdir -p atac-seq-pipeline/src
COPY src atac-seq-pipeline/src/
COPY atac.wdl atac-seq-pipeline/
ENV PATH="/software/atac-seq-pipeline:/software/atac-seq-pipeline/src:${PATH}"

# make some temporary directories
RUN for i in $(seq 0 9); do mkdir -p /mnt/ext_$i; done

#ENTRYPOINT ["/bin/bash","-c"]
