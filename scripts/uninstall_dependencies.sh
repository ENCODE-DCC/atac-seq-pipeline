#!/bin/bash
# Stop on error
set -e

CONDA_ENV=encode-atac-seq-pipeline
CONDA_ENV_PY3=encode-atac-seq-pipeline-python3

if which conda; then
  echo "=== Found Conda ($(conda --version))."
else
  echo "=== Conda does not exist on your system. Please install Conda first."
  echo "https://conda.io/docs/user-guide/install/index.html#regular-installation"
  exit 1
fi

if conda env list | grep -wq ${CONDA_ENV_PY3}; then
  echo "=== Removing pipeline's py3 Conda env (${CONDA_ENV_PY3})..."
  conda env remove -n ${CONDA_ENV_PY3} -y
else
  echo "=== Pipeline's py3 Conda env (${CONDA_ENV_PY3}) does not exist or has already been removed."
fi

if conda env list | grep -wq ${CONDA_ENV}; then
  echo "=== Removing pipeline's Conda env (${CONDA_ENV})..."
  conda env remove -n ${CONDA_ENV} -y
else
  echo "=== Pipeline's Conda env (${CONDA_ENV}) does not exist or has already been removed."
fi
