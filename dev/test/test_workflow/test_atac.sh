#!/bin/bash
set -e # exit on error

if [ $# -lt 2 ]; then
  echo "Usage: ./test_atac.sh [INPUT_JSON] [GCLOUD_SERVICE_ACCOUNT_SECRET_JSON_FILE] [DOCKER_IMAGE](optional)"
  exit 1
fi
if [ $# -gt 2 ]; then
  DOCKER_PARAM=$3
else
  DOCKER_PARAM="conda"
fi
INPUT=$1
GCLOUD_SERVICE_ACCOUNT_SECRET_JSON_FILE=$2
PREFIX=$(basename $INPUT .json)



CROMWELL_JAR="cromwell-42.jar"
if [ -f ${CROMWELL_JAR} ]; then
  echo "Skip downloading cromwell."
else
  wget -N -c https://storage.googleapis.com/encode-pipeline-test-samples/cromwell_jar/cromwell-42.jar 
fi

METADATA=${PREFIX}.metadata.json # metadata

cp $GCLOUD_SERVICE_ACCOUNT_SECRET_JSON_FILE tmp_secret_key.json
caper run ../../../atac.wdl -i ${INPUT} ${WF_OPT} -m ${METADATA}

java -Dconfig.file=backend_gcp_service_account.conf \
-Dbackend.default=google \
-Dbackend.providers.google.config.project=encode-dcc-1016 \
-Dbackend.providers.google.config.root="gs://encode-pipeline-test-runs/circleci" \
-Dbackend.providers.google.config.genomics.auth=service-account \
-Dbackend.providers.google.config.filesystems.gcs.auth=service-account \
-jar ${CROMWELL_JAR} run \
../../../atac.wdl \
-i ${INPUT} ${WF_OPT} -m ${METADATA}
 
rm -f tmp_secret_key
