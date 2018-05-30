#!/bin/bash

if [ -f "cromwell-32.jar" ]; then
  echo "Skip downloading cromwell."
else
  wget -N -c https://github.com/broadinstitute/cromwell/releases/download/32/cromwell-32.jar
fi
CROMWELL_JAR=cromwell-32.jar
BACKEND_CONF=../backends/backend_with_db.conf
BACKEND=google
GC_PROJ=encode-dcc-1016
GC_ROOT=gs://encode-pipeline-test-runs

java -Dconfig.file=${BACKEND_CONF} -Dbackend.default=${BACKEND} -Dbackend.providers.google.config.project=${GC_PROJ} \
-Dbackend.providers.google.config.root=${GC_ROOT} -jar ${CROMWELL_JAR} server
