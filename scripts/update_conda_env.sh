#!/bin/bash
# Stop on error
set -e

SH_SCRIPT_DIR=$(cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd)

CONDA_BIN="${CONDA_PREFIX}/bin"
cd ${CONDA_BIN}
chmod u+rx ${SH_SCRIPT_DIR}/../src/*.py
# copy all python scripts in /src into conda env bin dir
cp -f ${SH_SCRIPT_DIR}/../src/*.py ${CONDA_BIN}/
chmod u+rx ${CONDA_BIN}/encode_*.py

echo "=== Done."
