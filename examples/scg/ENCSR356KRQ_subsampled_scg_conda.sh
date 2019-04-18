#!/bin/bash

# do not touch these settings
#  number of tasks and nodes are fixed at 1
#SBATCH -n 1
#SBATCH --ntasks-per-node=1

# job name for pipeline
#  this name will appear when you monitor jobs with "squeue -u $USER"
#SBATCH --job-name=ENCSR356KRQ_subsampled

# walltime for your job
#  give long time enough to finish your pipeline
#  <12 hr: small/test samples
#  >24 hr: large samples
#SBATCH --time=12:00:00

# total amount of memory
#  depends on the size of your FASTQs
#  but should be <= NUM_CONCURRENT_TASK x 20GB for big samples
#  or <= NUM_CONCURRENT_TASK x 10GB for small samples
#  do not request too much memory
#  cluster will not accept your job
#SBATCH --mem=20G

# max number of cpus for each pipeline
#  should be <= NUM_CONCURRENT_TASK x "atac.bowtie2_cpu" in input JSON file
#  since bowtie2 is a bottlenecking task in the pipeline
#  "atac.bowtie2_cpu" is a number of cpus per replicate
#SBATCH --cpus-per-task=2

# email notification for job status
#SBATCH --mail-type=END,FAIL

# load java module if it exists
module load java
module load miniconda/3

# activate pipeline's Conda environment if Conda env exists
source activate encode-atac-seq-pipeline

# use input JSON for a small test sample
#  you make an input JSON for your own sample
#  start from any of two templates for single-ended and paired-ended samples
#  (examples/template_se.json, examples/template_pe.json)
#  do not use an input JSON file for a test sample (ENCSR356KRQ)
#  it's a sample with multimapping reads
INPUT=examples/scg/ENCSR356KRQ_subsampled_scg.json

# If this pipeline fails, then use this metadata JSON file to resume a failed pipeline from where it left 
# See details in /utils/resumer/README.md
PIPELINE_METADATA=metadata.json

# limit number of concurrent tasks
#  we recommend to use a number of replicates here
#  so that all replicates are processed in parellel at the same time.
#  make sure that resource settings in your input JSON file
#  are consistent with SBATCH resource settings (--mem, --cpus-per-task) 
#  in this script
NUM_CONCURRENT_TASK=2

# run pipeline
#  you can monitor your jobs with "squeue -u $USER"
java -jar -Dconfig.file=backends/backend.conf \
-Dbackend.providers.Local.config.concurrent-job-limit=${NUM_CONCURRENT_TASK} \
$HOME/cromwell-34.jar run atac.wdl -i ${INPUT} -m ${PIPELINE_METADATA}