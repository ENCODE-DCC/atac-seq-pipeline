#!/bin/bash
set -e # exit on error

if [ $# -lt 2 ]; then
  echo "Usage: ./test.sh [WDL] [INPUT_JSON] [DOCKER_IMAGE](optional)"
  echo "Make sure to have cromwell-30.1.jar in your \$PATH as an executable (chmod +x)."
  exit 1
fi

WDL=$1
INPUT=$2
if [ $# -gt 2 ]; then
  DOCKER_IMAGE=quay.io/encode-dcc/atac-seq-pipeline:latest
else
  DOCKER_IMAGE=$3
fi
CROMWELL_JAR=$(which cromwell-30.1.jar)
BACKEND_CONF=../../backends/backend.conf
BACKEND=Local
WF_OPT=../../workflow_opts/docker.json
PREFIX=$(basename ${WDL} .wdl)
METADATA=${PREFIX}.metadata.json # metadata
RESULT=${PREFIX}.result.json # output

java -Dconfig.file=${BACKEND_CONF} -Dbackend.default=${BACKEND} -jar ${CROMWELL_JAR} run ${WDL} -i ${INPUT} -o ${WF_OPT} -m ${METADATA}

# parse output metadata json
cat ${METADATA} | python -c "import json,sys;obj=json.load(sys.stdin);print(obj['outputs']['${PREFIX}.compare_md5sum.json_str'])" > ${RESULT}
cat ${RESULT}
rm -f ${METADATA}
