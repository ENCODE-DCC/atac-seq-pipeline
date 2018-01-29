#!/bin/bash
set -e # exit on error

CROMWELL_SVR_URL=35.185.235.240:8000
DOCKER_IMAGE=quay.io/encode-dcc/atac-seq-pipeline:latest
WDL=../../atac.wdl

if [ $# -lt 2 ]; then
  echo "Usage: ./test.atac.sh [INPUT_JSON] [QC_JSON_TO_COMPARE] [DOCKER_IMAGE](optional)"
  exit 1
fi
if [ $# -gt 2 ]; then
  DOCKER_IMAGE=$3
fi
INPUT=$1
QC_JSON_TO_COMPARE=$2
PREFIX=$(basename $INPUT .json)

# Write workflow option JSON file
TMP_WF_OPT=$PREFIX.test_atac_wf_opt.json
cat > $TMP_WF_OPT << EOM
{
    "default_runtime_attributes" : {
        "docker" : "$DOCKER_IMAGE",
        "cpu": "1",
        "memory": "4 GB",
        "disks": "local-disk 100 HDD",
        "zones": "us-west1-a us-west1-b us-west1-c us-central1-c us-central1-b",
        "failOnStderr" : false,
        "continueOnReturnCode" : 0,
        "preemptible": "0",
        "bootDiskSizeGb": "10",
        "noAddress": "false"
    }
}
EOM

# Submit workflow (POST)
curl -X POST --header "Accept: application/json" -v "$CROMWELL_SVR_URL/api/workflows/v1" \
-F workflowSource=@$WDL \
-F workflowInputs=@$INPUT \
-F workflowOptions=@$TMP_WF_OPT > $PREFIX.submit.json
rm -f $TMP_WF_OPT

# Get workflow-id
WF_ID=$(cat $PREFIX.submit.json | python -c 'import json,sys;obj=json.load(sys.stdin);print(obj["id"])')
rm -f $PREFIX.submit.json
echo "Workflow ID: $WF_ID"
sleep 30

# Check status of running job every 300 second
ITER=0
ITER_MAX=40
while true; do
  ITER=$(($ITER+1))
  # Get status (GET)
  curl -X GET --header "Accept: application/json" -v "$CROMWELL_SVR_URL/api/workflows/v1/$WF_ID/status" > $PREFIX.status.json
  WF_STATUS=$(cat $PREFIX.status.json | python -c 'import json,sys;obj=json.load(sys.stdin);print(obj["status"])')
  echo "Workflow status: $WF_STATUS, Iter: $ITER"
  if [ $WF_STATUS == Succeeded ]; then
  	echo "Workflow has been done successfully."
    break
  elif [ $WF_STATUS == Failed ]; then
    echo "Workflow has failed. Check out $PREFIX.status.json."
    cat $PREFIX.status.json
    exit 2
  fi
  if [ $ITER -gt $ITER_MAX ]; then
  	echo "Iteration reached limit ($ITER_MAX). Failed."
  	exit 3
    break
  fi
  rm -f $PREFIX.status.json
  sleep 300
done 

# Get output of workflow
curl -X GET --header "Accept: application/json" -v "$CROMWELL_SVR_URL/api/workflows/v1/$WF_ID/outputs" > $PREFIX.result.json
cat $PREFIX.result.json | python -c "import json,sys;obj=json.load(sys.stdin);print(obj['outputs']['atac.qc_report.qc_json_str'])" > $PREFIX.result.qc.json

echo "Done testing successfully."
diff $PREFIX.result.qc.json $QC_JSON_TO_COMPARE > $PREFIX.qc_json_diff.txt