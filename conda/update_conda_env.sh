#!/bin/bash
# Stop on error
set -e

CONDA_ENV=encode-atac-seq-pipeline

SH_SCRIPT_DIR=$(cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd)

if which conda; then
  echo "=== Found Conda ($(conda --version))."
else
  echo "=== Conda does not exist on your system. Please install Conda first."
  echo "=== https://conda.io/docs/user-guide/install/index.html#regular-installation"
  exit 1
fi

if conda env list | grep -wq ${CONDA_ENV}; then
  echo "=== Found Pipeline's Conda env (${CONDA_ENV})."
else
  echo "=== Pipeline's Conda env (${CONDA_ENV}) does not exist. Please install it first."
  echo "=== Run install_dependencies.sh"
  exit 2
fi

echo "=== Updating pipeline's source code on Conda env (${CONDA_ENV})..."
source activate ${CONDA_ENV}
  CONDA_BIN="${CONDA_PREFIX}/bin"
  cd ${CONDA_BIN}

  # copy all python scripts in /src into conda env bin dir
  cp -f ${SH_SCRIPT_DIR}/../src/*.py ${CONDA_BIN}/
  chmod +rx ${CONDA_BIN}/*.py
source deactivate

echo "=== All done."