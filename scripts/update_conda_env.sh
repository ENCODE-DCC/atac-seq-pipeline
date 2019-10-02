#!/bin/bash
set -e  # Stop on error

SH_SCRIPT_DIR=$(cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd)

CONDA_ENV_PY3=encode-atac-seq-pipeline
SRC_DIR=${SH_SCRIPT_DIR}/../src

conda --version  # check if conda exists

CONDA_PREFIX_PY3=$(conda env list | grep -P "\b${CONDA_ENV_PY3}\s" | awk '{if (NF==3) print $3; else print $2}')

if [ ! "${CONDA_PREFIX_PY3}" ];
then
	echo "Error: Pipeline's Conda environments not found."
	echo "Try to reinstall pipeline's Conda environments."
	echo
	echo "1) $ bash uninstall_conda_env.sh"
	echo "2) $ bash install_conda_env.sh"
	exit 1
fi

CONDA_BIN="${CONDA_PREFIX_PY3}/bin"

echo "=== Updating pipeline's Conda environments ==="
cd ${CONDA_BIN}
chmod u+rx ${SRC_DIR}/*.py
# copy all python scripts in /src into conda env bin dir
cp -f ${SRC_DIR}/*.py ${CONDA_BIN}/
chmod u+rx ${CONDA_BIN}/encode_*.py

echo "=== All done successfully ==="
