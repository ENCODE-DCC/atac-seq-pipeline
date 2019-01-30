# Tutorial for Stanford SCG4 cluster

All test samples and genome data are shared on Stanford SCG cluster based on SLURM. You don't have to download any data for testing our pipeline on it.

1. SSH to SCG's login node.
    ```bash
    $ ssh login.scg.stanford.edu
    ```

2. Download [cromwell](https://github.com/broadinstitute/cromwell).
    ```bash
    $ cd
    $ wget https://github.com/broadinstitute/cromwell/releases/download/34/cromwell-34.jar
    $ chmod +rx cromwell-34.jar
    ```

3. Git clone this pipeline and move into it.
    ```bash
    $ cd
    $ git clone https://github.com/ENCODE-DCC/atac-seq-pipeline
    $ cd atac-seq-pipeline
    ```

4. Set your account in `workflow_opts/scg.json`. [PIPELINE WILL NOT WORK WITHOUT A SLURM ACCOUNT](https://web.stanford.edu/group/scgpm/cgi-bin/informatics/wiki/index.php/Getting_A_Cluster_Account). Ignore other runtime attributes for singularity.
    ```javascript
    {
        "default_runtime_attributes" : {
            "slurm_account" : "YOUR_SLURM_ACCOUNT"
        }
    }
    ```

Our pipeline supports both [Conda](https://conda.io/docs/) and [Singularity](https://singularity.lbl.gov/).

## For Conda users

5. [Install Conda](https://conda.io/miniconda.html)

6. Install Conda dependencies.
    ```bash
    $ bash conda/uninstall_dependencies.sh  # to remove any existing pipeline env
    $ bash conda/install_dependencies.sh
    ```

7. Run a pipeline for a SUBSAMPLED (1/400) paired-end sample of [ENCSR356KRQ](https://www.encodeproject.org/experiments/ENCSR356KRQ/).
    ```bash
    $ source activate encode-atac-seq-pipeline # IMPORTANT!
    $ INPUT=examples/scg/ENCSR356KRQ_subsampled_scg.json
    $ java -jar -Xmx1G -Dconfig.file=backends/backend.conf -Dbackend.default=slurm cromwell-34.jar run atac.wdl -i ${INPUT} -o workflow_opts/scg.json
    ```

8. It will take about an hour. You will be able to find all outputs on `cromwell-executions/atac/[RANDOM_HASH_STRING]/`. See [output directory structure](output.md) for details.

9. See full specification for [input JSON file](input.md).

## For singularity users

5. Pull a singularity container for the pipeline. This will pull pipeline's docker container first and build a singularity one on `~/.singularity`.
    ```bash
    $ sdev    # SCG cluster does not allow building a container on login node
    $ mkdir -p ~/.singularity && cd ~/.singularity && SINGULARITY_CACHEDIR=~/.singularity SINGULARITY_PULLFOLDER=~/.singularity singularity pull --name atac-seq-pipeline-v1.1.7.simg -F docker://quay.io/encode-dcc/atac-seq-pipeline:v1.1.7
    $ exit
    ```

6. Run a pipeline for a SUBSAMPLED (1/400) paired-end sample of [ENCSR356KRQ](https://www.encodeproject.org/experiments/ENCSR356KRQ/).
    ```bash
    $ INPUT=examples/scg/ENCSR356KRQ_subsampled_scg.json
    $ java -jar -Xmx1G -Dconfig.file=backends/backend.conf -Dbackend.default=slurm_singularity cromwell-34.jar run atac.wdl -i ${INPUT} -o workflow_opts/scg.json
    ```

7. It will take about an hour. You will be able to find all outputs on `cromwell-executions/atac/[RANDOM_HASH_STRING]/`. See [output directory structure](output.md) for details.

8. See full specification for [input JSON file](input.md).

9. IF YOU WANT TO RUN PIPELINES WITH YOUR OWN INPUT DATA/GENOME DATABASE, PLEASE ADD THEIR DIRECTORIES TO `workflow_opts/scg.json`. For example, you have input FASTQs on `/your/input/fastqs/` and genome database installed on `/your/genome/database/` then add `/your/` to `singularity_bindpath`. You can also define multiple directories there. It's comma-separated.
    ```javascript
    {
        "default_runtime_attributes" : {
            "singularity_container" : "~/.singularity/chip-seq-pipeline-v1.1.7.simg",
            "singularity_bindpath" : "/scratch/users,/srv/gsfs0,/your/,YOUR_OWN_DATA_DIR1,YOUR_OWN_DATA_DIR1,..."
        }
    }
    ```

## Running multiple pipelines with cromwell server mode

1. If you want to run multiple (>10) pipelines, then run a cromwell server on an interactive node. We recommend to use `screen` or `tmux` to keep your session alive and note that all running pipelines will be killed after walltime. Run a Cromwell server with the following commands.

    ```bash
    $ srun -n 2 --mem 5G -t 3-0 --qos normal --account [YOUR_SCG_ACCOUNT] --pty /bin/bash -i -l    # 2 CPU, 5 GB RAM and 3 day walltime
    $ hostname -f    # to get [CROMWELL_SVR_IP]
    ```

    For Conda users,
    ```bash
    $ source activate encode-atac-seq-pipeline 
    $ _JAVA_OPTIONS="-Xmx5G" java -jar -Dconfig.file=backends/backend.conf -Dbackend.default=slurm cromwell-34.jar server
    ```

    For singularity users,
    ```bash
    $ _JAVA_OPTIONS="-Xmx5G" java -jar -Dconfig.file=backends/backend.conf -Dbackend.default=slurm_singularity cromwell-34.jar server
    ```


2. You can modify `backend.providers.slurm.concurrent-job-limit` or `backend.providers.slurm_singularity.concurrent-job-limit` in `backends/backend.conf` to increase maximum concurrent jobs. This limit is **NOT PER SAMPLE**. It's for all sub-tasks of all submitted samples.

3. On a login node, submit jobs to the cromwell server. You will get `[WORKFLOW_ID]` as a return value. Keep these workflow IDs for monitoring pipelines and finding outputs for a specific sample later.  
    ```bash  
    $ INPUT=YOUR_INPUT.json
    $ curl -X POST --header "Accept: application/json" -v "[CROMWELL_SVR_IP]:8000/api/workflows/v1" \
      -F workflowSource=@atac.wdl \
      -F workflowInputs=@${INPUT} \
      -F workflowOptions=@workflow_opts/scg.json
    ```

  To monitor pipelines, see [cromwell server REST API description](http://cromwell.readthedocs.io/en/develop/api/RESTAPI/#cromwell-server-rest-api>) for more details. `squeue` will not give you enough information for monitoring jobs per sample.
    ```bash
    $ curl -X GET --header "Accept: application/json" -v "[CROMWELL_SVR_IP]:8000/api/workflows/v1/[WORKFLOW_ID]/status"
    ```
