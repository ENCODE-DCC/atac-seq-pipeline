Tutorial for Stanford SCG4 cluster
==========================================

This tutorial shows how to run pipelines on SCG4. You may need to have a paid account on it because SCG4 does not offer any free of charge SLURM partition. We recommend that free users use [Sherlock](tutorial_sherlock.md) instead.

All test samples and genome data are shared on Stanford SCG4 cluster based on SLURM. You don't have to download any data for testing our pipeline on it.

1. SSH to SCG's login node.
    ```bash
      $ ssh login.scg.stanford.edu
    ```

2. Git clone this pipeline and move into it.
    ```bash
      $ git clone https://github.com/ENCODE-DCC/atac-seq-pipeline
      $ cd atac-seq-pipeline
    ```

3. Download [cromwell](https://github.com/broadinstitute/cromwell).
    ```bash
      $ wget https://github.com/broadinstitute/cromwell/releases/download/34/cromwell-34.jar
      $ chmod +rx cromwell-34.jar
    ```

Our pipeline supports both [Conda](https://conda.io/docs/) and [Singularity](https://singularity.lbl.gov/).

## For Conda users

4. [Install Conda](https://conda.io/miniconda.html)

5. Install Conda dependencies.
    ```bash
      $ bash conda/uninstall_dependencies.sh  # to remove any existing pipeline env
      $ bash conda/install_dependencies.sh
    ```

6. Make a shell script with name `run_test_pipeline.sh` on pipeline's git directory. This script will run without any modification, but read through it carefully if you want to run pipelines for your own sample which will be much bigger than the test sample.
    ```bash
      #!/bin/bash

      # job name for pipeline
      #  the shorter the better
      #  this name will appear when you monitor jobs with "squeue -u $USER"
      #SBATCH --job-name=ENCSR356KRQ_subsampled

      # do not touch these settings
      #  number of tasks and nodes are fixed at 1
      #SBATCH -n 1
      #SBATCH --ntasks-per-node=1

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

      # your paid account on SCG4
      #  free account doesn't seem to work
      #  this is not a partition specified with -p
      #  but you can also specify a partition if you need      
      #SBATCH --account [YOUR_PAID_ACCOUNT_ON_SCG4]

      # email notification about job status
      #SBATCH --mail-type=END,FAIL
      #SBATCH --mail-user=YOUR_EMAIL_ADDRESS

      # load java module
      module load java

      # activate pipeline's Conda environment
      source activate encode-atac-seq-pipeline

      # use input JSON for a small test sample
      #  you make an input JSON for your own sample
      #  start from any of two templates for single-ended and paired-ended samples
      #  (examples/template_se.json, examples/template_pe.json)
      #  do not use an input JSON file for a test sample (ENCSR356KRQ)
      #  it's a sample with multimapping reads
      INPUT=examples/scg/ENCSR356KRQ_subsampled_scg.json

      # limit number of concurrent tasks
      #  we recommend to use a number of replicates for it
      #  so that all replicates are processed in parellel at the same time.
      #  make sure that resource settings in your input JSON file
      #  are consistent with SBATCH resource settings (--mem, --cpus-per-task) 
      #  in this script
      NUM_CONCURRENT_TASK=2

      # run pipeline
      #  you can monitor your jobs with "squeue -u $USER"
      java -jar -Dconfig.file=backends/backend.conf \
      -Dbackend.providers.Local.config.concurrent-job-limit=${NUM_CONCURRENT_TASK} \
      ./cromwell-34.jar run atac.wdl -i ${INPUT}
    ```

## For singularity users

5. Pull a singularity container for the pipeline. This will pull pipeline's docker container first and build a singularity one on `~/.singularity`.
    ```bash
      $ sdev    # SCG cluster does not allow building a container on login node
      $ mkdir -p ~/.singularity && cd ~/.singularity && SINGULARITY_CACHEDIR=~/.singularity SINGULARITY_PULLFOLDER=~/.singularity singularity pull --name atac-seq-pipeline-v1.1.5.simg -F docker://quay.io/encode-dcc/atac-seq-pipeline:v1.1.5
      $ exit
    ```

6. Make a shell script with name `run_test_pipeline.sh` on pipeline's git directory. This script will run without any modification, but read through it carefully if you want to run pipelines for your own sample which will be much bigger than the test sample.
    ```bash
      #!/bin/bash

      # job name for pipeline
      #  the shorter the better
      #  this name will appear when you monitor jobs with "squeue -u $USER"
      #SBATCH --job-name=ENCSR356KRQ_subsampled

      # do not touch these settings
      #  number of tasks and nodes are fixed at 1
      #SBATCH -n 1
      #SBATCH --ntasks-per-node=1

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

      # your paid account on SCG4
      #  free account doesn't seem to work
      #  this is not a partition specified with -p
      #  but you can also specify a partition if you need
      #SBATCH --account [YOUR_PAID_ACCOUNT_ON_SCG4]

      # email notification about job status
      #SBATCH --mail-type=END,FAIL
      #SBATCH --mail-user=YOUR_EMAIL_ADDRESS

      # load java module
      module load java

      # use input JSON for a small test sample
      #  you make an input JSON for your own sample
      #  start from any of two templates for single-ended and paired-ended samples
      #  (examples/template_se.json, examples/template_pe.json)
      #  do not use an input JSON file for a test sample (ENCSR356KRQ)
      #  it's a sample with multimapping reads
      INPUT=examples/scg/ENCSR356KRQ_subsampled_scg.json

      # limit number of concurrent jobs
      #  we recommend to use a number of replicates for it
      #  so that all replicates are processed in parellel at the same time.
      #  make sure that resource settings in your input JSON file
      #  are consistent with SBATCH resource settings in this script
      #    * total memory: --mem=
      #    * max number of cpus: -cpus-per-task
      NUM_CONCURRENT_TASK=2

      # run pipeline
      #  you can monitor your jobs with "squeue -u $USER"
      java -jar -Dconfig.file=backends/backend.conf -Dbackend.default=singularity \
      -Dbackend.providers.singularity.config.concurrent-job-limit=${NUM_CONCURRENT_TASK} \
      ./cromwell-34.jar run atac.wdl -i ${INPUT} -o workflow_opts/singularity.json
    ```

## For all users

7. Run a pipeline for a SUBSAMPLED (1/400) paired-end sample of [ENCSR356KRQ](https://www.encodeproject.org/experiments/ENCSR356KRQ/).
    ```bash
      $ sbatch run_test_pipeline.sh
    ```

8. It will take about an hour. You will be able to find all outputs on `cromwell-executions/atac/[RANDOM_HASH_STRING]/`. See [output directory structure](output.md) for details. You can monitor your jobs with the following command:
    ```bash
      $ squeue -u $USER
    ```

9. See full specification for [input JSON file](input.md).

## For singularity users

10. IF YOU WANT TO RUN PIPELINES WITH YOUR OWN INPUT DATA/GENOME DATABASE, PLEASE ADD THEIR DIRECTORIES TO `workflow_opts/singularity.json`. For example, you have input FASTQs on `/your/input/fastqs/` and genome database installed on `/your/genome/database/` then add `/your/` to `singularity_bindpath`. You can also define multiple directories there. It's comma-separated.
    ```javascript
      {
          "default_runtime_attributes" : {
              "singularity_container" : "~/.singularity/atac-seq-pipeline-v1.1.5.simg",
              "singularity_bindpath" : "/scratch/users,/srv/gsfs0,/your/,YOUR_OWN_DATA_DIR1,YOUR_OWN_DATA_DIR1,..."
          }
      }
    ```
