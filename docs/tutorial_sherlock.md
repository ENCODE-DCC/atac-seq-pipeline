# Tutorial for Stanford Sherlock cluster

This tutorial shows how to run pipelines on Sherlock.

All test samples and genome data are shared on Stanford Sherlock cluster based on SLURM. You don't have to download any data for testing our pipeline on it.

1. SSH to Sherlock's login node.
    ```bash
    $ ssh login.sherlock.stanford.edu
    ```

2. Download [cromwell](https://github.com/broadinstitute/cromwell) on your `$HOME` directory.
    ```bash
    $ cd 
    $ wget https://github.com/broadinstitute/cromwell/releases/download/34/cromwell-34.jar
    $ chmod +rx cromwell-34.jar
    ```

3. Git clone this pipeline and move into its directory.
    ```bash
    $ cd
    $ git clone https://github.com/ENCODE-DCC/atac-seq-pipeline
    $ cd atac-seq-pipeline
    ```

Our pipeline supports both [Conda](https://conda.io/docs/) and [Singularity](https://singularity.lbl.gov/).

## For Conda users

4. [Install Conda](https://conda.io/miniconda.html). Skip this if you already have equivalent Conda alternatives (Anaconda Python). Download and run the [installer](https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh). Agree to the license term by typing `yes`. It will ask you about the installation location. On Stanford clusters (Sherlock and SCG4), we recommend to install it outside of your `$HOME` directory since its filesystem is slow and has very limited space. At the end of the installation, choose `yes` to add Miniconda's binary to `$PATH` in your BASH startup script.
    ```bash
    $ wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
    $ bash Miniconda3-latest-Linux-x86_64.sh
    ```

5. Install Conda dependencies.
    ```bash
    $ bash conda/uninstall_dependencies.sh  # to remove any existing pipeline env
    $ bash conda/install_dependencies.sh
    ```

6. Run a pipeline for the test sample.
    ```bash
    $ sbatch --partition normal examples/sherlock/ENCSR356KRQ_subsampled_sherlock_conda.sh
    ```

## For singularity users

6. Run a pipeline for the test sample.
    ```bash
    $ sbatch --partition normal examples/sherlock/ENCSR356KRQ_subsampled_sherlock_singularity.sh
    ```

## For all users

7. It will take about an hour. You will be able to find all outputs on `cromwell-executions/atac/[RANDOM_HASH_STRING]/`. See [output directory structure](output.md) for details. You can monitor your jobs with the following command:
    ```bash
    $ squeue -u $USER
    ```

8. See full specification for [input JSON file](input.md).

9. You can resume a failed pipeline from where it left off by using `PIPELINE_METADATA`(`metadata.json`) file. This file is created for each pipeline run. See [here](../utils/resumer/README.md) for details. Once you get a new input JSON file from the resumer, then edit your shell script (`examples/sherlock/ENCSR356KRQ_subsampled_sherlock_*.sh`) to use it `INPUT=resume.[FAILED_WORKFLOW_ID].json` instead of `INPUT=examples/...`.

## For singularity users

10. IF YOU WANT TO RUN PIPELINES WITH YOUR OWN INPUT DATA/GENOME DATABASE, PLEASE ADD THEIR DIRECTORIES TO `workflow_opts/sherlock.json`. For example, you have input FASTQs on `/your/input/fastqs/` and genome database installed on `/your/genome/database/` then add `/your/` to `singularity_bindpath`. You can also define multiple directories there. It's comma-separated.
    ```javascript
    {
        "default_runtime_attributes" : {
            "singularity_container" : "~/.singularity/atac-seq-pipeline-v1.3.0.simg",
            "singularity_bindpath" : "/scratch,/lscratch,/oak/stanford,/home/groups/cherry/encode,/your/,YOUR_OWN_DATA_DIR1,YOUR_OWN_DATA_DIR1,..."
        }
    }
    ```
