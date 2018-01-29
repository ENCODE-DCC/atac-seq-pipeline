#!/bin/bash

# make sure that chmod +x cromwell-30.1 and it's in $PATH
CROMWELL_JAR=$(which cromwell-30.1.jar)
BACKEND_CONF=../backends/backend_with_db.conf
BACKEND=google
GC_PROJ=encode-dcc-1016
GC_ROOT=gs://encode-pipeline-test-runs

java -Dconfig.file=${BACKEND_CONF} -Dbackend.default=${BACKEND} -Dbackend.providers.google.config.project=${GC_PROJ} \
-Dbackend.providers.google.config.root=${GC_ROOT} -jar ${CROMWELL_JAR} server