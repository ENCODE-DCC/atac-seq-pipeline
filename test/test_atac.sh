#!/bin/bash
set -e # exit on error

CROMWELL_SVR_URL=35.185.235.240:8000
TEST_WDL=test_atac.wdl
WDL=../atac.wdl
STANDALONE_WDL=test_atac_standalone.wdl
INPUT=test_atac.json
WF_OPT=test_atac_wf_opt.json

# Append task block of WDL to TEST_WDL
LINE_TASK_BEGIN=$(grep -n '@TASK_DEF_BEGIN' $WDL | awk -F":" '{print $1}')
LINE_TASK_END=$(grep -n '@TASK_DEF_END' $WDL | awk -F":" '{print $1}')
cp -f $TEST_WDL $STANDALONE_WDL
echo >> $STANDALONE_WDL
sed -n $LINE_TASK_BEGIN','$LINE_TASK_END'p' $WDL >> $STANDALONE_WDL

# Submit workflow (POST)
curl -X POST --header "Accept: application/json" -v "$CROMWELL_SVR_URL/api/workflows/v1" \
-F workflowSource=@$STANDALONE_WDL \
-F workflowInputs=@$INPUT \
-F workflowOptions=@$WF_OPT > tmp.submit.json
# curl -X POST --header "Accept: application/json" -v "localhost:8000/api/workflows/v1" -F workflowSource=@test_cromwell_server.wdl -F workflowInputs=@test_cromwell_server.json > tmp.submit.json
rm -f $STANDALONE_WDL

# Get workflow-id
WF_ID=$(cat tmp.submit.json | python -c 'import json,sys;obj=json.load(sys.stdin);print(obj["id"])')
rm -f tmp.submit.json
echo "Workflow ID: $WF_ID"
sleep 30

# Check status of running job every 300 second
ITER=0
ITER_MAX=20
while true; do
  ITER=$(($ITER+1))
  # Get status (GET)
  WF_STATUS=$(curl -X GET --header "Accept: application/json" -v "$CROMWELL_SVR_URL/api/workflows/v1/$WF_ID/status" | \
  	python -c 'import json,sys;obj=json.load(sys.stdin);print(obj["status"])')  
  echo "Workflow status: $WF_STATUS, Iter: $ITER"
  if [ $WF_STATUS == Succeeded ]; then
  	echo "Workflow has been done successfully."
    break
  fi
  if [ $ITER -gt $ITER_MAX ]; then
  	echo "Iteration reached limit ($ITER_MAX). Failed."
  	exit 1
    break
  fi
  sleep 300
done 

# Get output of workflow
curl -X GET --header "Accept: application/json" -v "$CROMWELL_SVR_URL/api/workflows/v1/$WF_ID/outputs" > tmp.result.json
RESULT=$(cat tmp.result.json | python -c 'import json,sys;obj=json.load(sys.stdin);print(obj["outputs"]["test_atac.compare_md5sum.match_overall"])')
cat tmp.result.json | python -c 'import json,sys;obj=json.load(sys.stdin);print(obj["outputs"]["test_atac.compare_md5sum.json_str"])' > result.json
rm -f tmp.result.json

echo "Test result (match_overall): $RESULT"
echo $RESULT > result.txt
echo "Test result (detail, JSON file): result.json"