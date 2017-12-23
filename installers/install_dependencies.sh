#!/bin/bash
# Stop on error
set -e

echo "=== Installing packages for python3 env..."
conda install nomkl -y
conda install --file requirements_py3.txt -y \
  -c defaults -c bioconda -c r -c bcbio -c daler -c asmeurer

echo "=== Installing packages for python2 env..."
conda create -n encode-atac-seq-pipelinepipeline --file requirements.txt \
  -y -c defaults -c bioconda -c r -c bcbio -c daler -c asmeurer

echo "=== Installing additional packages for python3 env..."
  CONDA_BIN=$(dirname $(which activate))
  cd $CONDA_BIN

  # uninstall IDR 2.0.3 and install the latest one
  conda uninstall idr -y
  rm -rf idr_tmp
  git clone --branch 2.0.4.1 https://github.com/kundajelab/idr idr_tmp
  cd idr_tmp
  $CONDA_BIN/python3 setup.py install
  cd $CONDA_BIN
  rm -rf idr_tmp

echo "=== Installing additional packages for python2 env..."
source activate encode-atac-seq-pipelinepipeline
  CONDA_BIN=$(dirname $(which activate))
  CONDA_LIB="$CONDA_BIN/../lib"
  cd $CONDA_BIN

  # soft-link IDR
  ln -s ../../../bin/idr
  ln -s ../../../bin/python3
  ln -s ../../../bin/pip3

  # graphviz in bioconda has segmentation fault bug
  conda uninstall -y graphviz -y
  conda install -y graphviz -c anaconda
  conda install -y ucsc-bedgraphtobigwig -c bioconda
  conda install -y ucsc-bedtobigbed -c bioconda

  # install Anshul's phantompeakqualtool 1.2
  wget -N -c https://github.com/kundajelab/phantompeakqualtools/archive/1.2.tar.gz
  tar -zxvf 1.2.tar.gz
  mv phantompeakqualtools-1.2/*.R .

  # # install R-spp 1.4
  # Rscript -e "install.packages('phantompeakqualtools-1.2/spp_1.14.tar.gz')"
  
  # install picard 2.10.6
  wget -N -c https://github.com/broadinstitute/picard/releases/download/2.10.6/picard.jar
  chmod +x picard.jar

  # resolve permission issue for python libraries
  if [[ $(find $CONDA_LIB -name '*egg-info*' -not -perm -o+r | wc -l ) > 0 ]]; then
    find $CONDA_LIB -name '*egg-info*' -not -perm -o+r -exec dirname {} \; | xargs chmod o+r -R
  fi

source deactivate

echo "=== All done."