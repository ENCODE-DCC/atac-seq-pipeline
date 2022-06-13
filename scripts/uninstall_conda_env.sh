#!/bin/bash

PIPELINE_CONDA_ENVS=(
  encd-atac
  encd-atac-macs2
  encd-atac-spp
  encd-atac-py2
)
for PIPELINE_CONDA_ENV in "${PIPELINE_CONDA_ENVS[@]}"
do
  conda env remove -n ${PIPELINE_CONDA_ENV} -y
done
