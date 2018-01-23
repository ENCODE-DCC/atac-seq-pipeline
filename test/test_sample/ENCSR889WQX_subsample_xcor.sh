#!/bin/bash

# make sure that chmod +x cromwell-30.1 and it's in $PATH
CROMWELL_JAR=$(which cromwell-30.1.jar)
WDL=../../atac.wdl
INPUT=ENCSR889WQX_subsample_xcor.json
WF_OPT=../../workflow_opts/docker_google.json
BACKEND_CONF=../../backends/backend.conf
BACKEND=google
GC_PROJ=encode-dcc-1016
GC_ROOT=gs://encode-pipeline-test-runs/test_atac_subsampled_sample/ENCSR889WQX_subsample_xcor

java -Dconfig.file=${BACKEND_CONF} -Dbackend.default=${BACKEND} -Dbackend.providers.google.config.project=${GC_PROJ} \
-Dbackend.providers.google.config.root=${GC_ROOT} -jar ${CROMWELL_JAR} run ${WDL} -i ${INPUT} -o ${WF_OPT}