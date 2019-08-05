# Dev

## Command line for version change
```bash
PREV_VER=v1.5.0
NEW_VER=v1.5.0
for f in $(grep -rl ${PREV_VER} --include=*.{wdl,md,sh})
do
  sed -i "s/${PREV_VER}/${NEW_VER}/g" ${f}
done
cd dev/workflow_opts
for f in $(grep -rl ${PREV_VER} --include=*.json)
do
  sed -i "s/${PREV_VER}/${NEW_VER}/g" ${f}
done
cd ../../
```

## Building templates on DX for each genome

Make sure that you have [`dxWDL-0.79.1.jar`](https://github.com/DNAnexus/dxWDL/releases/download/0.79.1/dxWDL-0.79.1.jar) on your `$HOME`. Install [DNAnexus Platform SDK](https://wiki.DNAnexus.com/downloads) with `pip install dxpy`. Log-in on DNAnexus with `dx login` and choose "ENCODE Uniform Processing Pipelines" (name of our official DNAnexus project for pipelines).

Run the following command line locally to build out DX workflows for this pipeline on our official one. This will overwrite (`-f` parameter does it).

```bash
# version
VER=v1.5.0

# general
java -jar ~/dxWDL-0.79.1.jar compile atac.wdl -project "ENCODE Uniform Processing Pipelines" -extras dev/workflow_opts/docker.json -f -folder /ATAC-seq/workflows/$VER/general -defaults dev/examples/dx/template_general.json

# hg38
java -jar ~/dxWDL-0.79.1.jar compile atac.wdl -project "ENCODE Uniform Processing Pipelines" -extras dev/workflow_opts/docker.json -f -folder /ATAC-seq/workflows/$VER/hg38 -defaults dev/examples/dx/template_hg38.json

# hg19
java -jar ~/dxWDL-0.79.1.jar compile atac.wdl -project "ENCODE Uniform Processing Pipelines" -extras dev/workflow_opts/docker.json -f -folder /ATAC-seq/workflows/$VER/hg19 -defaults dev/examples/dx/template_hg19.json

# mm10
java -jar ~/dxWDL-0.79.1.jar compile atac.wdl -project "ENCODE Uniform Processing Pipelines" -extras dev/workflow_opts/docker.json -f -folder /ATAC-seq/workflows/$VER/mm10 -defaults dev/examples/dx/template_mm10.json

# mm9
java -jar ~/dxWDL-0.79.1.jar compile atac.wdl -project "ENCODE Uniform Processing Pipelines" -extras dev/workflow_opts/docker.json -f -folder /ATAC-seq/workflows/$VER/mm9 -defaults dev/examples/dx/template_mm9.json

# test sample
java -jar ~/dxWDL-0.79.1.jar compile atac.wdl -project "ENCODE Uniform Processing Pipelines" -extras dev/workflow_opts/docker.json -f -folder /ATAC-seq/workflows/$VER/test_ENCSR356KRQ_subsampled -defaults dev/examples/dx/ENCSR356KRQ_subsampled_dx.json

## DX Azure

# general
java -jar ~/dxWDL-0.79.1.jar compile atac.wdl -project "ENCODE Uniform Processing Pipelines Azure" -extras dev/workflow_opts/docker.json -f -folder /ATAC-seq/workflows/$VER/general -defaults dev/examples/dx_azure/template_general.json

# hg38
java -jar ~/dxWDL-0.79.1.jar compile atac.wdl -project "ENCODE Uniform Processing Pipelines Azure" -extras dev/workflow_opts/docker.json -f -folder /ATAC-seq/workflows/$VER/hg38 -defaults dev/examples/dx_azure/template_hg38.json

# hg19
java -jar ~/dxWDL-0.79.1.jar compile atac.wdl -project "ENCODE Uniform Processing Pipelines Azure" -extras dev/workflow_opts/docker.json -f -folder /ATAC-seq/workflows/$VER/hg19 -defaults dev/examples/dx_azure/template_hg19.json

# mm10
java -jar ~/dxWDL-0.79.1.jar compile atac.wdl -project "ENCODE Uniform Processing Pipelines Azure" -extras dev/workflow_opts/docker.json -f -folder /ATAC-seq/workflows/$VER/mm10 -defaults dev/examples/dx_azure/template_mm10.json

# mm9
java -jar ~/dxWDL-0.79.1.jar compile atac.wdl -project "ENCODE Uniform Processing Pipelines Azure" -extras dev/workflow_opts/docker.json -f -folder /ATAC-seq/workflows/$VER/mm9 -defaults dev/examples/dx_azure/template_mm9.json

# test sample
java -jar ~/dxWDL-0.79.1.jar compile atac.wdl -project "ENCODE Uniform Processing Pipelines Azure" -extras dev/workflow_opts/docker.json -f -folder /ATAC-seq/workflows/$VER/test_ENCSR356KRQ_subsampled -defaults dev/examples/dx_azure/ENCSR356KRQ_subsampled_dx_azure.json
```

## Building genome data (for mito-only index)

```
conda/build_genome_data.sh hg38 /mnt/lab_data/pipeline_genome_data/hg38
conda/build_genome_data.sh hg38_chr19_chrM /mnt/lab_data/pipeline_genome_data/hg38_chr19_chrM
conda/build_genome_data.sh mm10 /mnt/lab_data/pipeline_genome_data/mm10
conda/build_genome_data.sh mm10_chr19_chrM /mnt/lab_data/pipeline_genome_data/mm10_chr19_chrM
conda/build_genome_data.sh hg19 /mnt/lab_data/pipeline_genome_data/hg19
conda/build_genome_data.sh mm9 /mnt/lab_data/pipeline_genome_data/mm9

cd /mnt/lab_data/pipeline_genome_data

find -name *.fasta -delete
find -name *.fa -delete

mv hg19/hg19.tsv hg19_klab.tsv
mv hg38_chr19_chrM/hg38_chr19_chrM.tsv hg38_chr19_chrM_klab.tsv
mv hg38/hg38.tsv hg38_klab.tsv
mv mm10_chr19_chrM/mm10_chr19_chrM.tsv mm10_chr19_chrM_klab.tsv
mv mm10/mm10.tsv mm10_klab.tsv 
mv mm9/mm9.tsv mm9_klab.tsv

DIR_SRC=/mnt/lab_data/pipeline_genome_data

PLATFORM=caper
DIR_TARGET=https://storage.googleapis.com/encode-pipeline-genome-data
for G in "hg19" "mm9" "hg38" "mm10" "hg38_chr19_chrM" "mm10_chr19_chrM";
do
  TSV1=${G}_klab.tsv
  TSV2=${G}_$PLATFORM.tsv
  cp $TSV1 $TSV2
  sed -i "s~$DIR_SRC~$DIR_TARGET~g" $TSV2
done

PLATFORM=gcp
DIR_TARGET=gs://encode-pipeline-genome-data
for G in "hg19" "mm9" "hg38" "mm10" "hg38_chr19_chrM" "mm10_chr19_chrM";
do
  TSV1=${G}_klab.tsv
  TSV2=${G}_$PLATFORM.tsv
  cp $TSV1 $TSV2
  sed -i "s~$DIR_SRC~$DIR_TARGET~g" $TSV2
  gsutil cp $DIR_SRC/$G/ataqc/* $DIR_TARGET/$G/ataqc/
  gsutil cp $DIR_SRC/$G/*.chrM.fa.gz $DIR_TARGET/$G/
  gsutil cp $DIR_SRC/$G/bowtie2_index/*.chrM.* $DIR_TARGET/$G/bowtie2_index/
  gsutil cp $DIR_SRC/$G/bwa_index/*.chrM.* $DIR_TARGET/$G/bwa_index/
done
gsutil cp *_gcp.tsv $DIR_TARGET/
gsutil cp *_caper.tsv $DIR_TARGET/

PLATFORM=aws
DIR_TARGET=s3://encode-pipeline-genome-data
for G in "hg19" "mm9" "hg38" "mm10" "hg38_chr19_chrM" "mm10_chr19_chrM";
do
  TSV1=${G}_klab.tsv
  TSV2=${G}_$PLATFORM.tsv
  cp $TSV1 $TSV2
  sed -i "s~$DIR_SRC~$DIR_TARGET~g" $TSV2
  gsutil cp $DIR_SRC/$G/ataqc/* $DIR_TARGET/$G/ataqc/
  gsutil cp $DIR_SRC/$G/*.chrM.fa.gz $DIR_TARGET/$G/
  gsutil cp $DIR_SRC/$G/bowtie2_index/*.chrM.* $DIR_TARGET/$G/bowtie2_index/
  gsutil cp $DIR_SRC/$G/bwa_index/*.chrM.* $DIR_TARGET/$G/bwa_index/
done
gsutil cp *_aws.tsv $DIR_TARGET/

PLATFORM=sherlock
PREFIX=leepc12@login.sherlock.stanford.edu
DIR_TARGET=/home/groups/cherry/encode/pipeline_genome_data
for G in "hg19" "mm9" "hg38" "mm10" "hg38_chr19_chrM" "mm10_chr19_chrM";
do
  TSV1=${G}_klab.tsv
  TSV2=${G}_$PLATFORM.tsv
  cp $TSV1 $TSV2
  sed -i "s~$DIR_SRC~$DIR_TARGET~g" $TSV2
  scp $DIR_SRC/$G/ataqc/* $PREFIX:$DIR_TARGET/$G/ataqc/
  scp $DIR_SRC/$G/*.chrM.fa.gz $PREFIX:$DIR_TARGET/$G/
  scp $DIR_SRC/$G/bowtie2_index/*.chrM.* $PREFIX:$DIR_TARGET/$G/bowtie2_index/
  scp $DIR_SRC/$G/bwa_index/*.chrM.* $PREFIX:$DIR_TARGET/$G/bwa_index/
done
scp *_sherlock.tsv $PREFIX:$DIR_TARGET/

PLATFORM=scg
PREFIX=leepc12@login.scg.stanford.edu
DIR_TARGET=/reference/ENCODE/pipeline_genome_data
for G in "hg19" "mm9" "hg38" "mm10" "hg38_chr19_chrM" "mm10_chr19_chrM";
do
  TSV1=${G}_klab.tsv
  TSV2=${G}_$PLATFORM.tsv
  cp $TSV1 $TSV2
  sed -i "s~$DIR_SRC~$DIR_TARGET~g" $TSV2
  scp $DIR_SRC/$G/ataqc/* $PREFIX:$DIR_TARGET/$G/ataqc/
  scp $DIR_SRC/$G/*.chrM.fa.gz $PREFIX:$DIR_TARGET/$G/
  scp $DIR_SRC/$G/bowtie2_index/*.chrM.* $PREFIX:$DIR_TARGET/$G/bowtie2_index/
  scp $DIR_SRC/$G/bwa_index/*.chrM.* $PREFIX:$DIR_TARGET/$G/bwa_index/
done
scp *_scg.tsv $PREFIX:$DIR_TARGET/

# LOG-IN ON ENCODE OFFICIAL DX (AWS) project 
PLATFORM=dx
DIR_TARGET=/pipeline-genome-data
for G in "hg19" "mm9" "hg38" "mm10" "hg38_chr19_chrM" "mm10_chr19_chrM";
do
  #dx upload $DIR_SRC/$G/ataqc/* --path $DIR_TARGET/$G/ataqc/
  dx upload $DIR_SRC/$G/*.chrM.fa.gz --path $DIR_TARGET/$G/
  dx upload $DIR_SRC/$G/bowtie2_index/*.chrM.fa.tar --path $DIR_TARGET/$G/bowtie2_index/
  dx upload $DIR_SRC/$G/bwa_index/*.chrM.fa.tar --path $DIR_TARGET/$G/bwa_index/
done

# REMOVE *.tsv files on DX (AWS) project dir /pipeline-genome-data
# COPY NEW TSVs to it
dx upload *_dx.tsv --path pipeline-genome-data/

# LOG-IN ON ENCODE OFFICIAL DX (AZURE) project 
PLATFORM=dx_azure
DIR_TARGET=/pipeline-genome-data
for G in "hg19" "mm9" "hg38" "mm10" "hg38_chr19_chrM" "mm10_chr19_chrM";
do
  #dx upload $DIR_SRC/$G/ataqc/* --path $DIR_TARGET/$G/ataqc/
  dx upload $DIR_SRC/$G/*.chrM.fa.gz --path $DIR_TARGET/$G/
  dx upload $DIR_SRC/$G/bowtie2_index/*.chrM.fa.tar --path $DIR_TARGET/$G/bowtie2_index/
  dx upload $DIR_SRC/$G/bwa_index/*.chrM.fa.tar --path $DIR_TARGET/$G/bwa_index/
done

# REMOVE *.tsv files on DX (AZURE) project dir /pipeline-genome-data
# COPY NEW TSVs to it
dx upload *_dx_azure.tsv --path pipeline-genome-data/

```

## Replacing strings in a deprecated method

```


cp ./scg/ENCSR356KRQ_subsampled_scg.json ../../example_input_json/scg/ENCSR356KRQ_subsampled.scg.json
cp ./dx/ENCSR356KRQ_subsampled_dx.json ../../example_input_json/dx/ENCSR356KRQ_subsampled.dx.json
cp ./google/ENCSR356KRQ_subsampled.json ../../example_input_json/gcp/ENCSR356KRQ_subsampled.gcp.json
cp ./sherlock/ENCSR356KRQ_subsampled_sherlock.json ../../example_input_json/ENCSR356KRQ_subsampled.sherlock.json
cp ./dx_azure/ENCSR356KRQ_subsampled_dx_azure.json ../../example_input_json/ENCSR356KRQ_subsampled.dx_azure.json
cp ./caper/ENCSR356KRQ_subsampled.json ../../example_input_json/

cp ./scg/ENCSR936XTK_subsampled_chr19_only_scg.json ../../example_input_json/scg/
cp ./dx/ENCSR936XTK_subsampled_chr19_only_dx.json ../../example_input_json/dx/
cp ./dx_azure/ENCSR000DYI_subsampled_chr19_only_dx_azure.json ../../example_input_json/dx_azure/
cp ./google/ENCSR936XTK_subsampled_chr19_only.json ../../example_input_json/gcp/ENCSR936XTK_subsampled_chr19_only_gcp.json
cp ./caper/ENCSR936XTK_subsampled_chr19_only.json ../../example_input_json/caper/ENCSR936XTK_subsampled_chr19_only_caper.json
cp ./sherlock/ENCSR936XTK_subsampled_chr19_only_sherlock.json ../../example_input_json/sherlock/
cp ./scg/ENCSR936XTK_subsampled_chr19_only_scg.json ../../example_input_json/scg/

for f in $(grep -rl "examples/" --include=*.md)
do
  sed -i "s/examples\//dev\/examples\//g" ${f}
done

for f in $(grep -rl "backends/backend" --include=*.md)
do
  sed -i "s/backends\/backend/dev\/backends\/backend/g" ${f}
done


for f in $(grep -rl "workflow_opts/" --include=*.md)
do
  sed -i "s/workflow_opts\//dev\/workflow_opts\//g" ${f}
done

```


