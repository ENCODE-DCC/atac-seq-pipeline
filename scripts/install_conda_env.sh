#!/bin/bash
set -e  # Stop on error

install_ucsc_tools_369() {
  # takes in conda env name and find conda bin
  CONDA_BIN=$(conda run -n $1 bash -c "echo \$(dirname \$(which python))")
  curl -o "$CONDA_BIN/fetchChromSizes" "https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v369/fetchChromSizes"
  curl -o "$CONDA_BIN/wigToBigWig" "https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v369/wigToBigWig"
  curl -o "$CONDA_BIN/bedGraphToBigWig" "https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v369/bedGraphToBigWig"
  curl -o "$CONDA_BIN/bigWigInfo" "https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v369/bigWigInfo"
  curl -o "$CONDA_BIN/bedClip" "https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v369/bedClip"
  curl -o "$CONDA_BIN/bedToBigBed" "https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v369/bedToBigBed"
  curl -o "$CONDA_BIN/twoBitToFa" "https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v369/twoBitToFa"
  curl -o "$CONDA_BIN/bigWigAverageOverBed" "https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v369/bigWigAverageOverBed"

  chmod +x "$CONDA_BIN/fetchChromSizes"
  chmod +x "$CONDA_BIN/wigToBigWig"
  chmod +x "$CONDA_BIN/bedGraphToBigWig"
  chmod +x "$CONDA_BIN/bigWigInfo"
  chmod +x "$CONDA_BIN/bedClip"
  chmod +x "$CONDA_BIN/bedToBigBed"
  chmod +x "$CONDA_BIN/twoBitToFa"
  chmod +x "$CONDA_BIN/bigWigAverageOverBed"
}

SH_SCRIPT_DIR=$(cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd)

echo "$(date): Installing pipeline's Conda environments..."

conda create -n encd-atac --file ${SH_SCRIPT_DIR}/requirements.txt \
  --override-channels -c bioconda -c defaults -y

conda create -n encd-atac-macs2 --file ${SH_SCRIPT_DIR}/requirements.macs2.txt \
  --override-channels -c bioconda -c defaults -y

conda create -n encd-atac-py2 --file ${SH_SCRIPT_DIR}/requirements.python2.txt \
  --override-channels -c conda-forge -c bioconda -c defaults -y

conda create -n encd-atac-spp --file ${SH_SCRIPT_DIR}/requirements.spp.txt \
  -c r -c bioconda -c defaults -y

# adhoc fix for the following issues:
# - https://github.com/ENCODE-DCC/chip-seq-pipeline2/issues/259
# - https://github.com/ENCODE-DCC/chip-seq-pipeline2/issues/265
# force-install readline 6.2, ncurses 5.9 from conda-forge (ignoring dependencies)
#conda install -n encd-atac-spp --no-deps --no-update-deps -y \
#  readline==6.2 ncurses==5.9 -c conda-forge

CONDA_BIN=$(conda run -n encd-atac-spp bash -c "echo \$(dirname \$(which python))")

echo "$(date): Installing phantompeakqualtools in Conda environments..."
RUN_SPP="https://raw.githubusercontent.com/kundajelab/phantompeakqualtools/1.2.2/run_spp.R"
conda run -n encd-atac-spp bash -c \
  "curl -o $CONDA_BIN/run_spp.R $RUN_SPP && chmod +x $CONDA_BIN/run_spp.R"

echo "$(date): Installing R packages in Conda environments..."
CRAN="https://cran.r-project.org/"
conda run -n encd-atac-spp bash -c \
  "Rscript -e \"install.packages('snow', repos='$CRAN')\""
conda run -n encd-atac-spp bash -c \
  "Rscript -e \"install.packages('snowfall', repos='$CRAN')\""
conda run -n encd-atac-spp bash -c \
  "Rscript -e \"install.packages('bitops', repos='$CRAN')\""
conda run -n encd-atac-spp bash -c \
  "Rscript -e \"install.packages('caTools', repos='$CRAN')\""
conda run -n encd-atac-spp bash -c \
  "Rscript -e \"install.packages('BiocManager', repos='$CRAN')\""
conda run -n encd-atac-spp bash -c \
  "Rscript -e \"require('BiocManager'); BiocManager::install('Rsamtools'); BiocManager::install('Rcpp')\""

echo "$(date): Installing R spp 1.15.5 in Conda environments..."
SPP="https://cran.r-project.org/src/contrib/Archive/spp/spp_1.15.5.tar.gz"
SPP_BASENAME=$(basename $SPP)
curl -o "$CONDA_BIN/$SPP_BASENAME" "$SPP"
conda run -n encd-atac-spp bash -c \
  "Rscript -e \"install.packages('$CONDA_BIN/$SPP_BASENAME')\""

echo "$(date): Installing USCS tools (v369)..."
install_ucsc_tools_369 encd-atac
install_ucsc_tools_369 encd-atac-spp
install_ucsc_tools_369 encd-atac-macs2

echo "$(date): Done successfully."
echo
echo "If you see openssl,readline,ncurses lib errors while running pipelines"
echo "then switch to Singularity method. Conda method will not work on your system."

bash ${SH_SCRIPT_DIR}/update_conda_env.sh
