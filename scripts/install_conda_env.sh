#!/bin/bash
set -e  # Stop on error

SH_SCRIPT_DIR=$(cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd)

echo "$(date): Installing pipeline's Conda environments..."

conda create -n encode-atac-seq-pipeline --file ${SH_SCRIPT_DIR}/requirements.txt \
  --override-channels -c bioconda -c defaults -y

conda create -n encode-atac-seq-pipeline-macs2 --file ${SH_SCRIPT_DIR}/requirements.macs2.txt \
  --override-channels -c bioconda -c defaults -y

conda create -n encode-atac-seq-pipeline-spp --file ${SH_SCRIPT_DIR}/requirements.spp.txt \
  --override-channels -c r -c bioconda -c defaults -y

# adhoc fix for https://github.com/ENCODE-DCC/chip-seq-pipeline2/issues/259
# force-install readline 6.2 from conda-forge (ignoring dependencies)
conda install -n encode-atac-seq-pipeline-spp --no-deps --no-update-deps -y readline==6.2 -c conda-forge

conda create -n encode-atac-seq-pipeline-python2 --file ${SH_SCRIPT_DIR}/requirements.python2.txt \
  --override-channels -c conda-forge -c bioconda -c defaults -y

echo "$(date): Done successfully."

bash ${SH_SCRIPT_DIR}/update_conda_env.sh
