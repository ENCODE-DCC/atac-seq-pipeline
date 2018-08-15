#!/bin/bash
# Stop on error
set -e

CONDA_ENV=encode-atac-seq-pipeline
CONDA_ENV_PY3=encode-atac-seq-pipeline-python3

SH_SCRIPT_DIR=$(cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd)

REQ_TXT=${SH_SCRIPT_DIR}/requirements.txt
REQ_TXT_PY3=${SH_SCRIPT_DIR}/requirements_py3.txt

if which conda; then
  echo "=== Found Conda ($(conda --version))."
else
  echo "=== Conda does not exist on your system. Please install Conda first."
  echo "https://conda.io/docs/user-guide/install/index.html#regular-installation"
  exit 1
fi

if conda env list | grep -wq ${CONDA_ENV_PY3}; then
  echo "=== Pipeline's py3 Conda env (${CONDA_ENV_PY3}) already exists."
  echo "=== Please remove it first (uninstall_dependencies.sh)"
  exit 2
fi

if conda env list | grep -wq ${CONDA_ENV}; then
  echo "=== Pipeline's Conda env (${CONDA_ENV}) already exists."
  echo "=== Please remove it first (uninstall_dependencies.sh)"
  exit 3
fi

echo "=== Installing packages for python3 env..."
conda create -n ${CONDA_ENV_PY3} --file ${REQ_TXT_PY3} -y -c bioconda -c conda-forge -c defaults -c r

echo "=== Installing packages for python2 env..."
conda create -n ${CONDA_ENV} --file ${REQ_TXT} -y -c bioconda -c conda-forge -c defaults -c r

echo "=== Installing additional packages for python2 env..."
source activate ${CONDA_ENV}
  CONDA_BIN=$(dirname $(which bedtools))
  CONDA_SHARE="${CONDA_BIN}/../share"

  #hack around the need for both python2 and python3 
  #in the same environment
  cd ${CONDA_BIN}
  rm -f idr
  ln -s ../../${CONDA_ENV_PY3}/bin/idr

  # make an executable symlink for cromwell.jar on conda bin dir
  chmod +rx ${CONDA_SHARE}/cromwell/cromwell.jar
  cd ${CONDA_BIN}
  ln -s ../share/cromwell/cromwell.jar

source deactivate

# update pipeline's python scripts (src/*.py) on conda environment
cd ${SH_SCRIPT_DIR}
bash update_conda_env.sh

echo "=== All done."
