# Dev

## Command line for version change
```bash
PREV_VER=v1.1.7
NEW_VER=v1.1.7
for f in $(grep -rl ${PREV_VER} --include=*.{wdl,md,sh,yml})
do
  sed -i "s/${PREV_VER}/${NEW_VER}/g" ${f}
done
cd workflow_opts
for f in $(grep -rl ${PREV_VER} --include=*.json)
do
  sed -i "s/${PREV_VER}/${NEW_VER}/g" ${f}
done
cd ..
```

## Building templates on DX for each genome

Make sure that you have [`dxWDL-0.77.jar`](https://github.com/DNAnexus/dxWDL/releases/download/0.77/dxWDL-0.77.jar) on your `$HOME`. Install [DNAnexus Platform SDK](https://wiki.DNAnexus.com/downloads) with `pip install dxpy`. Log-in on DNAnexus with `dx login` and choose "ENCODE Uniform Processing Pipelines" (name of our official DNAnexus project for pipelines).

Run the following command line locally to build out DX workflows for this pipeline on our official one. This will overwrite (`-f` parameter does it).

```bash
# version
VER=v1.1.7

# general
java -jar ~/dxWDL-0.77.jar compile atac.wdl -project "ENCODE Uniform Processing Pipelines" -extras workflow_opts/docker.json -f -folder /ATAC-seq/workflows/$VER/general -defaults examples/dx/template_general.json

# hg38
java -jar ~/dxWDL-0.77.jar compile atac.wdl -project "ENCODE Uniform Processing Pipelines" -extras workflow_opts/docker.json -f -folder /ATAC-seq/workflows/$VER/hg38 -defaults examples/dx/template_hg38.json

# hg19
java -jar ~/dxWDL-0.77.jar compile atac.wdl -project "ENCODE Uniform Processing Pipelines" -extras workflow_opts/docker.json -f -folder /ATAC-seq/workflows/$VER/hg19 -defaults examples/dx/template_hg19.json

# mm10
java -jar ~/dxWDL-0.77.jar compile atac.wdl -project "ENCODE Uniform Processing Pipelines" -extras workflow_opts/docker.json -f -folder /ATAC-seq/workflows/$VER/mm10 -defaults examples/dx/template_mm10.json

# mm9
java -jar ~/dxWDL-0.77.jar compile atac.wdl -project "ENCODE Uniform Processing Pipelines" -extras workflow_opts/docker.json -f -folder /ATAC-seq/workflows/$VER/mm9 -defaults examples/dx/template_mm9.json

# test sample
java -jar ~/dxWDL-0.77.jar compile atac.wdl -project "ENCODE Uniform Processing Pipelines" -extras workflow_opts/docker.json -f -folder /ATAC-seq/workflows/$VER/test_ENCSR356KRQ_subsampled -defaults examples/dx/ENCSR356KRQ_subsampled_dx.json

## DX Azure

# general
java -jar ~/dxWDL-0.77.jar compile atac.wdl -project "ENCODE Uniform Processing Pipelines Azure" -extras workflow_opts/docker.json -f -folder /ATAC-seq/workflows/$VER/general -defaults examples/dx_azure/template_general.json

# hg38
java -jar ~/dxWDL-0.77.jar compile atac.wdl -project "ENCODE Uniform Processing Pipelines Azure" -extras workflow_opts/docker.json -f -folder /ATAC-seq/workflows/$VER/hg38 -defaults examples/dx_azure/template_hg38.json

# hg19
java -jar ~/dxWDL-0.77.jar compile atac.wdl -project "ENCODE Uniform Processing Pipelines Azure" -extras workflow_opts/docker.json -f -folder /ATAC-seq/workflows/$VER/hg19 -defaults examples/dx_azure/template_hg19.json

# mm10
java -jar ~/dxWDL-0.77.jar compile atac.wdl -project "ENCODE Uniform Processing Pipelines Azure" -extras workflow_opts/docker.json -f -folder /ATAC-seq/workflows/$VER/mm10 -defaults examples/dx_azure/template_mm10.json

# mm9
java -jar ~/dxWDL-0.77.jar compile atac.wdl -project "ENCODE Uniform Processing Pipelines Azure" -extras workflow_opts/docker.json -f -folder /ATAC-seq/workflows/$VER/mm9 -defaults examples/dx_azure/template_mm9.json

# test sample
java -jar ~/dxWDL-0.77.jar compile atac.wdl -project "ENCODE Uniform Processing Pipelines Azure" -extras workflow_opts/docker.json -f -folder /ATAC-seq/workflows/$VER/test_ENCSR356KRQ_subsampled -defaults examples/dx_azure/ENCSR356KRQ_subsampled_dx_azure.json
```
