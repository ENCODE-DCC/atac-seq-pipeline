#!/bin/bash
set -e

VER=$(cat atac.wdl | grep "String pipeline_ver = " | awk '{gsub("'"'"'",""); print $4}')
DOCKER=encodedcc/atac-seq-pipeline:$VER
DXWDL=~/dxWDL-v1.46.4.jar

# check if docker image exists on dockerhub
docker pull $DOCKER

# make a copy of WDL
# workflows generated from original WDL cannot start from inputs other than FASTQs
# make all vars non-optional (without ?)
cp atac.wdl atac.dx.wdl
sed -i 's/Array\[File?\] bams = \[\]/Array\[File\] bams = \[\]/g' atac.dx.wdl
sed -i 's/Array\[File?\] nodup_bams = \[\]/Array\[File\] nodup_bams = \[\]/g' atac.dx.wdl
sed -i 's/Array\[File?\] tas = \[\]/Array\[File\] tas = \[\]/g' atac.dx.wdl

# general
java -jar ${DXWDL} compile atac.dx.wdl -project "ENCODE Uniform Processing Pipelines" -extras \
<(echo "{\"default_runtime_attributes\":{\"docker\":\"${DOCKER}\"}}") -f -folder \
/ATAC-seq/workflows/$VER/general -defaults example_input_json/dx/template_general.json

# hg38
java -jar ${DXWDL} compile atac.dx.wdl -project "ENCODE Uniform Processing Pipelines" -extras \
<(echo "{\"default_runtime_attributes\":{\"docker\":\"${DOCKER}\"}}") -f -folder \
/ATAC-seq/workflows/$VER/hg38 -defaults example_input_json/dx/template_hg38.json

# hg19
java -jar ${DXWDL} compile atac.dx.wdl -project "ENCODE Uniform Processing Pipelines" -extras \
<(echo "{\"default_runtime_attributes\":{\"docker\":\"${DOCKER}\"}}") -f -folder \
/ATAC-seq/workflows/$VER/hg19 -defaults example_input_json/dx/template_hg19.json

# mm10
java -jar ${DXWDL} compile atac.dx.wdl -project "ENCODE Uniform Processing Pipelines" -extras \
<(echo "{\"default_runtime_attributes\":{\"docker\":\"${DOCKER}\"}}") -f -folder \
/ATAC-seq/workflows/$VER/mm10 -defaults example_input_json/dx/template_mm10.json

# mm9
java -jar ${DXWDL} compile atac.dx.wdl -project "ENCODE Uniform Processing Pipelines" -extras \
<(echo "{\"default_runtime_attributes\":{\"docker\":\"${DOCKER}\"}}") -f -folder \
/ATAC-seq/workflows/$VER/mm9 -defaults example_input_json/dx/template_mm9.json

# test sample
java -jar ${DXWDL} compile atac.dx.wdl -project "ENCODE Uniform Processing Pipelines" -extras \
<(echo "{\"default_runtime_attributes\":{\"docker\":\"${DOCKER}\"}}") -f -folder \
/ATAC-seq/workflows/$VER/test_ENCSR356KRQ_subsampled -defaults example_input_json/dx/ENCSR356KRQ_subsampled_dx.json

# test sample (single rep)
java -jar ${DXWDL} compile atac.dx.wdl -project "ENCODE Uniform Processing Pipelines" -extras \
<(echo "{\"default_runtime_attributes\":{\"docker\":\"${DOCKER}\"}}") -f -folder \
/ATAC-seq/workflows/$VER/test_ENCSR356KRQ_subsampled_rep1 -defaults example_input_json/dx/ENCSR356KRQ_subsampled_rep1_dx.json

## DX Azure

# general
java -jar ${DXWDL} compile atac.dx.wdl -project "ENCODE Uniform Processing Pipelines Azure" -extras \
<(echo "{\"default_runtime_attributes\":{\"docker\":\"${DOCKER}\"}}") -f -folder \
/ATAC-seq/workflows/$VER/general -defaults example_input_json/dx_azure/template_general.json

# hg38
java -jar ${DXWDL} compile atac.dx.wdl -project "ENCODE Uniform Processing Pipelines Azure" -extras \
<(echo "{\"default_runtime_attributes\":{\"docker\":\"${DOCKER}\"}}") -f -folder \
/ATAC-seq/workflows/$VER/hg38 -defaults example_input_json/dx_azure/template_hg38.json

# hg19
java -jar ${DXWDL} compile atac.dx.wdl -project "ENCODE Uniform Processing Pipelines Azure" -extras \
<(echo "{\"default_runtime_attributes\":{\"docker\":\"${DOCKER}\"}}") -f -folder \
/ATAC-seq/workflows/$VER/hg19 -defaults example_input_json/dx_azure/template_hg19.json

# mm10
java -jar ${DXWDL} compile atac.dx.wdl -project "ENCODE Uniform Processing Pipelines Azure" -extras \
<(echo "{\"default_runtime_attributes\":{\"docker\":\"${DOCKER}\"}}") -f -folder \
/ATAC-seq/workflows/$VER/mm10 -defaults example_input_json/dx_azure/template_mm10.json

# mm9
java -jar ${DXWDL} compile atac.dx.wdl -project "ENCODE Uniform Processing Pipelines Azure" -extras \
<(echo "{\"default_runtime_attributes\":{\"docker\":\"${DOCKER}\"}}") -f -folder \
/ATAC-seq/workflows/$VER/mm9 -defaults example_input_json/dx_azure/template_mm9.json

# test sample
java -jar ${DXWDL} compile atac.dx.wdl -project "ENCODE Uniform Processing Pipelines Azure" -extras \
<(echo "{\"default_runtime_attributes\":{\"docker\":\"${DOCKER}\"}}") -f -folder \
/ATAC-seq/workflows/$VER/test_ENCSR356KRQ_subsampled -defaults example_input_json/dx_azure/ENCSR356KRQ_subsampled_dx_azure.json

rm -f atac.dx.wdl
