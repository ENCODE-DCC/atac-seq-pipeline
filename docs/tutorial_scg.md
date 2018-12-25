# Tutorial for Stanford SCG4 cluster


This tutorial shows how to run pipelines on SCG4. You may need to have a paid account on it because SCG4 does not offer any free of charge SLURM partition. We recommend that free users use [Sherlock](tutorial_sherlock.md) instead.

All test samples and genome data are shared on Stanford SCG4 cluster based on SLURM. You don't have to download any data for testing our pipeline on it.

1. Download [cromwell](https://github.com/broadinstitute/cromwell) on your `$HOME` directory.
    ```bash
    $ cd 
    $ wget https://github.com/broadinstitute/cromwell/releases/download/34/cromwell-34.jar
    $ chmod +rx cromwell-34.jar
    ```

2. SSH to SCG's login node.
    ```bash
    $ ssh login.scg.stanford.edu
    ```

3. Git clone this pipeline and move into its directory.
    ```bash
    $ cd
    $ git clone https://github.com/ENCODE-DCC/atac-seq-pipeline
    $ cd atac-seq-pipeline
    ```

Our pipeline supports both [Conda](https://conda.io/docs/) and [Singularity](https://singularity.lbl.gov/).

## For Conda users

4. [Install Conda](https://conda.io/miniconda.html)

5. Install Conda dependencies.
    ```bash
    $ bash conda/uninstall_dependencies.sh  # to remove any existing pipeline env
    $ bash conda/install_dependencies.sh
    ```

6. Run a pipeline for the test sample. You must have a paid account on SCG4.
    ```bash
    $ sbatch --account [YOUR_PAID_ACCOUNT_ON_SCG4] examples/scg/ENCSR356KRQ_subsampled_scg_conda.sh
    ```

## For singularity users

5. Pull a singularity container for the pipeline. This will pull pipeline's docker container first and build a singularity one on `~/.singularity`.
    ```bash
    $ sdev    # SCG cluster does not allow building a container on login node
    $ rm -rf ~/.singularity && mkdir -p ~/.singularity && cd ~/.singularity && SINGULARITY_CACHEDIR=~/.singularity SINGULARITY_PULLFOLDER=~/.singularity singularity pull --name atac-seq-pipeline-v1.1.5.simg -F docker://quay.io/encode-dcc/atac-seq-pipeline:v1.1.5
    $ exit
    ```

6. Run a pipeline for the test sample. You must have a paid account on SCG4.
    ```bash
    $ sbatch --account [YOUR_PAID_ACCOUNT_ON_SCG4] examples/scg/ENCSR356KRQ_subsampled_scg_singularity.sh
    ```

## For all users

7. It will take about an hour. You will be able to find all outputs on `cromwell-executions/atac/[RANDOM_HASH_STRING]/`. See [output directory structure](output.md) for details. You can monitor your jobs with the following command:
    ```bash
    $ squeue -u $USER
    ```

8. See full specification for [input JSON file](input.md).

## For singularity users

9. IF YOU WANT TO RUN PIPELINES WITH YOUR OWN INPUT DATA/GENOME DATABASE, PLEASE ADD THEIR DIRECTORIES TO `workflow_opts/scg.json`. For example, you have input FASTQs on `/your/input/fastqs/` and genome database installed on `/your/genome/database/` then add `/your/` to `singularity_bindpath`. You can also define multiple directories there. It's comma-separated.
    ```javascript
    {
        "default_runtime_attributes" : {
            "singularity_container" : "~/.singularity/atac-seq-pipeline-v1.1.5.simg",
            "singularity_bindpath" : "/reference/ENCODE,/scratch,/srv/gsfs0,YOUR_OWN_DATA_DIR1,YOUR_OWN_DATA_DIR1,..."
        }
    }
    ```
