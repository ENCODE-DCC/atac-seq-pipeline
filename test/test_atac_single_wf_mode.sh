#!/bin/bash

# make sure that chmod +x cromwell-30.1 and it's in $PATH
CROMWELL_JAR=$(which cromwell-30.1.jar)
TEST_WDL=test_atac.wdl
WDL=../atac.wdl
STANDALONE_WDL=test_atac_standalone.wdl
INPUT=test_atac.json
WF_OPT=../workflow_opts/docker_google.json
BACKEND_CONF=../backends/backend.conf
BACKEND=google
GC_PROJ=encode-dcc-1016
GC_ROOT=gs://encode-pipeline-test-runs/test_atac

# Append task block of WDL to TEST_WDL
LINE_TASK_BEGIN=$(grep -n '@TASK_DEF_BEGIN' $WDL | awk -F":" '{print $1}')
LINE_TASK_END=$(grep -n '@TASK_DEF_END' $WDL | awk -F":" '{print $1}')
cp -f $TEST_WDL $STANDALONE_WDL
echo >> $STANDALONE_WDL
sed -n $LINE_TASK_BEGIN','$LINE_TASK_END'p' $WDL >> $STANDALONE_WDL

java -Dconfig.file=${BACKEND_CONF} -Dbackend.default=${BACKEND} -Dbackend.providers.google.config.project=${GC_PROJ} \
-Dbackend.providers.google.config.root=${GC_ROOT} -jar ${CROMWELL_JAR} run ${STANDALONE_WDL} -i ${INPUT} -o ${WF_OPT}
