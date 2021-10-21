#!/bin/bash

PIPELINE_CONDA_ENVS=(
  encode-atac-seq-pipeline
  encode-atac-seq-pipeline-macs2
  encode-atac-seq-pipeline-spp
  encode-atac-seq-pipeline-python2
)
for PIPELINE_CONDA_ENV in "${PIPELINE_CONDA_ENVS[@]}"
do
  conda env remove -n ${PIPELINE_CONDA_ENV} -y
done
