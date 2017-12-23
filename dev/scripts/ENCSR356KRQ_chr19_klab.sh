WORKDIR=/srv/scratch/shared/kadru/leepc12/run/atac-seq-pipeline/cromwell/ENCSR356KRQ_chr19/from_fastq
mkdir -p $WORKDIR && cd $WORKDIR
source activate atac-seq-pipeline
java -jar -Dconfig.file=$CODE/atac-seq-pipeline/backends/default.conf $CODE/atac-seq-pipeline/cromwell-30-x.jar run $CODE/atac-seq-pipeline/atac.wdl -i $CODE/atac-seq-pipeline/examples/klab/ENCSR356KRQ_chr19_klab.json -o $CODE/atac-seq-pipeline/workflow_opts/non_docker.json -m output_ENCSR356KRQ.json
source deactivate