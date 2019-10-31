#!/bin/bash
set -e

# version
VER=$(cat atac.wdl | grep "#CAPER docker" | awk 'BEGIN{FS=":"} {print $2}')
DOCKER=quay.io/encode-dcc/atac-seq-pipeline:$VER

# general
java -jar ~/dxWDL-0.79.1.jar compile atac.wdl -project "ENCODE Uniform Processing Pipelines" -extras <(echo "{\"default_runtime_attributes\":{\"docker\":\"${DOCKER}\"}}") -f -folder /ATAC-seq/workflows/$VER/general -defaults example_input_json/dx/template_general.json

# hg38
java -jar ~/dxWDL-0.79.1.jar compile atac.wdl -project "ENCODE Uniform Processing Pipelines" -extras <(echo "{\"default_runtime_attributes\":{\"docker\":\"${DOCKER}\"}}") -f -folder /ATAC-seq/workflows/$VER/hg38 -defaults example_input_json/dx/template_hg38.json

# hg19
java -jar ~/dxWDL-0.79.1.jar compile atac.wdl -project "ENCODE Uniform Processing Pipelines" -extras <(echo "{\"default_runtime_attributes\":{\"docker\":\"${DOCKER}\"}}") -f -folder /ATAC-seq/workflows/$VER/hg19 -defaults example_input_json/dx/template_hg19.json

# mm10
java -jar ~/dxWDL-0.79.1.jar compile atac.wdl -project "ENCODE Uniform Processing Pipelines" -extras <(echo "{\"default_runtime_attributes\":{\"docker\":\"${DOCKER}\"}}") -f -folder /ATAC-seq/workflows/$VER/mm10 -defaults example_input_json/dx/template_mm10.json

# mm9
java -jar ~/dxWDL-0.79.1.jar compile atac.wdl -project "ENCODE Uniform Processing Pipelines" -extras <(echo "{\"default_runtime_attributes\":{\"docker\":\"${DOCKER}\"}}") -f -folder /ATAC-seq/workflows/$VER/mm9 -defaults example_input_json/dx/template_mm9.json

# test sample
java -jar ~/dxWDL-0.79.1.jar compile atac.wdl -project "ENCODE Uniform Processing Pipelines" -extras <(echo "{\"default_runtime_attributes\":{\"docker\":\"${DOCKER}\"}}") -f -folder /ATAC-seq/workflows/$VER/test_ENCSR356KRQ_subsampled -defaults example_input_json/dx/ENCSR356KRQ_subsampled_dx.json

# test sample (single rep)
java -jar ~/dxWDL-0.79.1.jar compile atac.wdl -project "ENCODE Uniform Processing Pipelines" -extras <(echo "{\"default_runtime_attributes\":{\"docker\":\"${DOCKER}\"}}") -f -folder /ATAC-seq/workflows/$VER/test_ENCSR356KRQ_subsampled_rep1 -defaults example_input_json/dx/ENCSR356KRQ_subsampled_rep1_dx.json

## DX Azure

# general
java -jar ~/dxWDL-0.79.1.jar compile atac.wdl -project "ENCODE Uniform Processing Pipelines Azure" -extras <(echo "{\"default_runtime_attributes\":{\"docker\":\"${DOCKER}\"}}") -f -folder /ATAC-seq/workflows/$VER/general -defaults example_input_json/dx_azure/template_general.json

# hg38
java -jar ~/dxWDL-0.79.1.jar compile atac.wdl -project "ENCODE Uniform Processing Pipelines Azure" -extras <(echo "{\"default_runtime_attributes\":{\"docker\":\"${DOCKER}\"}}") -f -folder /ATAC-seq/workflows/$VER/hg38 -defaults example_input_json/dx_azure/template_hg38.json

# hg19
java -jar ~/dxWDL-0.79.1.jar compile atac.wdl -project "ENCODE Uniform Processing Pipelines Azure" -extras <(echo "{\"default_runtime_attributes\":{\"docker\":\"${DOCKER}\"}}") -f -folder /ATAC-seq/workflows/$VER/hg19 -defaults example_input_json/dx_azure/template_hg19.json

# mm10
java -jar ~/dxWDL-0.79.1.jar compile atac.wdl -project "ENCODE Uniform Processing Pipelines Azure" -extras <(echo "{\"default_runtime_attributes\":{\"docker\":\"${DOCKER}\"}}") -f -folder /ATAC-seq/workflows/$VER/mm10 -defaults example_input_json/dx_azure/template_mm10.json

# mm9
java -jar ~/dxWDL-0.79.1.jar compile atac.wdl -project "ENCODE Uniform Processing Pipelines Azure" -extras <(echo "{\"default_runtime_attributes\":{\"docker\":\"${DOCKER}\"}}") -f -folder /ATAC-seq/workflows/$VER/mm9 -defaults example_input_json/dx_azure/template_mm9.json

# test sample
java -jar ~/dxWDL-0.79.1.jar compile atac.wdl -project "ENCODE Uniform Processing Pipelines Azure" -extras <(echo "{\"default_runtime_attributes\":{\"docker\":\"${DOCKER}\"}}") -f -folder /ATAC-seq/workflows/$VER/test_ENCSR356KRQ_subsampled -defaults example_input_json/dx_azure/ENCSR356KRQ_subsampled_dx_azure.json

