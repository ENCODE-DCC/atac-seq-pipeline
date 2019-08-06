#!/bin/bash
# Stop on error
set -e

CONDA_ENV_PY3=encode-atac-seq-pipeline-python3

echo "=== Installing additional packages for python2 env..."
# all numpy packages in the conda-forge repo are linked to OpenBLAS,
# which generate slightly different SVD results for finding summits in MACS2
# according to cpu arch, num threads.
# so we need to keep using numpy with MKL in the defaults repo
# until this issue is fixed.
conda install numpy==1.11.3 -c defaults -y

# init script for activation/deactivation.
# limiting number of threads for BLAS is important to ensure the same output for MACS2
# most pipeline users have their own numpy and scipy installed on their home.
# so we also need to set PYTHONNOUSERSITE=True to prevent python from 
# looking for external packages locally installed on user's home.
# see https://github.com/conda/conda/issues/6018 for details
CONDA_LIB="${CONDA_PREFIX}/lib"
CONDA_PY3_LIB="${CONDA_PREFIX}/../${CONDA_ENV_PY3}/lib"
CONDA_ACTIVATE_D="${CONDA_PREFIX}/etc/conda/activate.d"
CONDA_DEACTIVATE_D="${CONDA_PREFIX}/etc/conda/deactivate.d"
CONDA_ACTIVATE_SH="${CONDA_ACTIVATE_D}/env_vars.sh"
CONDA_DEACTIVATE_SH="${CONDA_DEACTIVATE_D}/env_vars.sh"
mkdir -p ${CONDA_ACTIVATE_D}
mkdir -p ${CONDA_DEACTIVATE_D}
echo "export OPENBLAS_NUM_THREADS=1" > ${CONDA_ACTIVATE_SH}
echo "export MKL_NUM_THREADS=1" >> ${CONDA_ACTIVATE_SH}
echo "export PYTHONNOUSERSITE=True" >> ${CONDA_ACTIVATE_SH}
echo "unset OPENBLAS_NUM_THREADS MKL_NUM_THREADS PYTHONNOUSERSITE" > ${CONDA_DEACTIVATE_SH}

# to prevent conflict between Conda's R and global(local) R
echo "export R_HOME=${CONDA_LIB}/R" >> ${CONDA_ACTIVATE_SH}
echo "export R_LIBS=${CONDA_LIB}/R/library" >> ${CONDA_ACTIVATE_SH}
echo "unset R_HOME R_LIBS" > ${CONDA_DEACTIVATE_SH}

# hack around the need for both python2 and python3 in the same environment
CONDA_BIN="${CONDA_PREFIX}/bin"
cd ${CONDA_BIN}
rm -f idr python3
ln -s ../../${CONDA_ENV_PY3}/bin/idr
ln -s ../../${CONDA_ENV_PY3}/bin/python3

# make an executable symlink for cromwell.jar on conda bin dir
CONDA_SHARE="${CONDA_PREFIX}/share"
chmod +rx ${CONDA_SHARE}/cromwell/cromwell.jar
cd ${CONDA_BIN}
rm -f cromwell.jar
ln -s ../share/cromwell/cromwell.jar

# install SAMstats 0.2.0 (float div bug fixed)
pip install --no-dependencies SAMstats==0.2.0

# install deeptools 3.3.0
pip install --no-dependencies deeptools==3.3.0 #deeptoolsintervals==0.1.8 pyBigWig==0.3.11 plotly retrying 

# make a soft link for picard.jar
if [ -f ${CONDA_PREFIX}/bin/picard ]; then
  PICARD_JAR=../share/picard*/picard.jar
  chmod +rx ${PICARD_JAR}
  cd ${CONDA_PREFIX}/bin && ln -s ${PICARD_JAR}
fi

echo "=== Done."
