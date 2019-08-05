#!/bin/bash
# Stop on error
set -e

CONDA_ENV=encode-atac-seq-pipeline
CONDA_ENV_PY3=encode-atac-seq-pipeline-python3

SH_SCRIPT_DIR=$(cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd)

REQ_TXT=${SH_SCRIPT_DIR}/requirements.txt
REQ_TXT_PY3=${SH_SCRIPT_DIR}/requirements_py3.txt

echo "=== Installing packages for python3 env..."
conda create -n ${CONDA_ENV_PY3} --file ${REQ_TXT_PY3} -y -c bioconda -c conda-forge -c defaults -c r

echo "=== Installing packages for python2 env..."
conda create -n ${CONDA_ENV} --file ${REQ_TXT} -y -c bioconda -c conda-forge -c defaults -c r

echo "=== Done."
