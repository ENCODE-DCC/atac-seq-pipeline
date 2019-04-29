#!/bin/bash
set -e # exit on error

if [ $# -lt 2 ]; then
  echo "Usage: ./test_atac.sh [INPUT_JSON] [GCLOUD_SERVICE_ACCOUNT_SECRET_JSON_FILE] [DOCKER_IMAGE](optional)"
  exit 1
fi
if [ $# -gt 2 ]; then
  DOCKER_IMAGE=$3
else
  DOCKER_IMAGE=quay.io/encode-dcc/atac-seq-pipeline:test-v1.3.0
fi
INPUT=$1
GCLOUD_SERVICE_ACCOUNT_SECRET_JSON_FILE=$2
PREFIX=$(basename $INPUT .json)

CROMWELL_JAR="cromwell-34.jar"
if [ -f ${CROMWELL_JAR} ]; then
  echo "Skip downloading cromwell."
else
  wget -N -c https://github.com/broadinstitute/cromwell/releases/download/34/cromwell-34.jar
fi

# Write workflow option JSON file
TMP_WF_OPT=$PREFIX.test_atac_wf_opt.json
cat > $TMP_WF_OPT << EOM
{
    "default_runtime_attributes" : {
        "docker" : "$DOCKER_IMAGE"
    }
}
EOM

METADATA=${PREFIX}.metadata.json # metadata
RESULT=${PREFIX}.result.txt # output

cp $GCLOUD_SERVICE_ACCOUNT_SECRET_JSON_FILE tmp_secret_key.json

java -Dconfig.file=backend_gcp_service_account.conf \
-Dbackend.default=google \
-Dbackend.providers.google.config.project=encode-dcc-1016 \
-Dbackend.providers.google.config.root="gs://encode-pipeline-test-runs/circleci" \
-Dbackend.providers.google.config.genomics.auth=service-account \
-Dbackend.providers.google.config.filesystems.gcs.auth=service-account \
-jar ${CROMWELL_JAR} run \
../../atac.wdl \
-i ${INPUT} -o ${TMP_WF_OPT} -m ${METADATA}
 
rm -f tmp_secret_key ${TMP_WF_OPT}
