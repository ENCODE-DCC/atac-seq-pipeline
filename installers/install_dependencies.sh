#!/bin/bash
# Stop on error
set -e

echo "=== Installing packages for python3 env..."
conda create -n encode-atac-seq-pipeline-python3 \
  --file requirements_py3.txt -y \
  -c bioconda -c conda-forge -c defaults -c r


echo "=== Installing packages for python2 env..."
conda create -n encode-atac-seq-pipeline --file requirements.txt \
  -y -c bioconda -c conda-forge -c defaults -c r

echo "=== Installing additional packages for python2 env..."
source activate encode-atac-seq-pipeline
  CONDA_BIN=$(dirname $(which bedtools))
  CONDA_LIB="$CONDA_BIN/../lib"

  #hack around the need for both python2 and python3 
  #in the same environment
  rm -f $CONDA_BIN/idr
  ln -s ../../encode-atac-seq-pipeline-python3/bin/idr $CONDA_BIN/idr

  # decompress MACS2 python egg to ensure functionality on Sherlock
  cd $CONDA_LIB/python2.7/site-packages
  unzip -o MACS2-2.1.1.20160309-py2.7-linux-x86_64.egg
  

  # resolve permission issue for python libraries
#  if [[ $(find $CONDA_LIB -name '*egg-info*' -not -perm -o+r | wc -l ) > 0 ]]; then
#    find $CONDA_LIB -name '*egg-info*' -not -perm -o+r -exec dirname {} \; | xargs chmod o+r -R
#  fi

source deactivate

echo "=== All done."
