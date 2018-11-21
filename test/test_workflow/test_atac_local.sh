#!/bin/bash
set -e # exit on error

if [ $# -lt 1 ]; then
  echo "Usage: ./test_atac_local.sh [INPUT_JSON] [DOCKER_IMAGE](optional)"
  exit 1
fi
if [ $# -gt 1 ]; then
  DOCKER_IMAGE=$2
else
  DOCKER_IMAGE=quay.io/encode-dcc/atac-seq-pipeline:v1.1.2
fi
INPUT=$1
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

WDL=../../atac.wdl
BACKEND_CONF=../../backends/backend.conf
BACKEND=google
PROJECT="encode-dcc-1016"
BUCKET="gs://encode-pipeline-test-runs/circleci"
PREFIX=$(basename ${INPUT} .json)
METADATA=${PREFIX}.metadata.json # metadata
RESULT=${PREFIX}.result.json # output

java -Dconfig.file=${BACKEND_CONF} -Dbackend.default=${BACKEND} \
-Dbackend.providers.google.config.project=${PROJECT} \
-Dbackend.providers.google.config.root=${BUCKET} \
-jar ${CROMWELL_JAR} run ${WDL} \
-i ${INPUT} -o ${TMP_WF_OPT} -m ${METADATA}

# parse output metadata json
cat ${METADATA} | python -c "import json,sys;obj=json.load(sys.stdin);print(obj['outputs']['${PREFIX}.compare_md5sum.json_str'])" > ${RESULT}
cat ${RESULT}

rm -f ${METADATA} ${TMP_WF_OPT}

