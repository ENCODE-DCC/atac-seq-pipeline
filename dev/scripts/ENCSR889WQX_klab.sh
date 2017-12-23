WORKDIR=/srv/scratch/shared/kadru/leepc12/run/atac-seq-pipeline/cromwell/ENCSR889WQX/from_fastq
mkdir -p $WORKDIR && cd $WORKDIR
source activate atac-seq-pipeline
java -jar -Dconfig.file=$CODE/atac-seq-pipeline/backends/default.conf $CODE/atac-seq-pipeline/cromwell-30-x.jar run $CODE/atac-seq-pipeline/atac.wdl -i $CODE/atac-seq-pipeline/examples/ENCSR889WQX_klab.json -o $CODE/atac-seq-pipeline/workflow_opts/non_docker.json -m output_ENCSR889WQX.json
source deactivate
