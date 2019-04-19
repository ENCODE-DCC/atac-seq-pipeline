# Tutorial for Sun GridEngine (SGE) clusters

1. Download [cromwell](https://github.com/broadinstitute/cromwell).
    ```bash
    $ cd
    $ wget https://github.com/broadinstitute/cromwell/releases/download/34/cromwell-34.jar
    $ chmod +rx cromwell-34.jar
    ```

2. Git clone this pipeline and move into it.
    ```bash
    $ git clone https://github.com/ENCODE-DCC/atac-seq-pipeline
    $ cd atac-seq-pipeline
    ```

3. Download a SUBSAMPLED (1/400) paired-end sample of [ENCSR356KRQ](https://www.encodeproject.org/experiments/ENCSR356KRQ/).
    ```bash
    $ wget https://storage.googleapis.com/encode-pipeline-test-samples/encode-atac-seq-pipeline/ENCSR356KRQ/ENCSR356KRQ_fastq_subsampled.tar
    $ tar xvf ENCSR356KRQ_fastq_subsampled.tar
    ```

4. Download pre-built genome database for hg38.
    ```bash
    $ wget https://storage.googleapis.com/encode-pipeline-genome-data/test_genome_database_hg38_atac.tar
    $ tar xvf test_genome_database_hg38_atac.tar
    ```

5. Set your parallel environment (PE) and queue in `workflow_opts/sge.json`. If your SGE cluster does not have any PE, ask your admin to add one for our pipeline. If you don't want to specify any queue then remove `, "sge_queue" : "YOUR_SGE_QUEUE"` from the file. See [here](how_to_config_sge.md) to find details about how to configure SGE for the pipeline.
    ```javascript
    {
        "default_runtime_attributes" : {
            "sge_pe" : "YOUR_SGE_PE",
            "sge_queue" : "YOUR_SGE_QUEUE"
        }
    }
    ```

Our pipeline supports both [Conda](https://conda.io/docs/) and [Singularity](https://singularity.lbl.gov/).

## For Conda users,

6. [Install Conda](https://conda.io/miniconda.html)

7. Install Conda dependencies.
    ```bash
    $ bash conda/uninstall_dependencies.sh  # to remove any existing pipeline env
    $ bash conda/install_dependencies.sh
    ```

8. Run a pipeline for the test sample.
    ```bash
    $ source activate encode-atac-seq-pipeline # IMPORTANT!
    $ INPUT=examples/local/ENCSR356KRQ_subsampled.json
    $ java -jar -Xmx1G -Dconfig.file=backends/backend.conf -Dbackend.default=sge cromwell-34.jar run atac.wdl -i ${INPUT} -o workflow_opts/sge.json
    ```

9. It will take about an hour. You will be able to find all outputs on `cromwell-executions/atac/[RANDOM_HASH_STRING]/`. See [output directory structure](output.md) for details.

10. See full specification for [input JSON file](input.md).

## For singularity users,

6. CHECK YOUR SINGULARITY VERSION FIRST AND UPGRADE IT TO A VERSION `>=2.5.2` OR PIPELINE WILL NOT WORK CORRECTLY.
    ```bash
    $ singularity --version
    ```

7. Pull a singularity container for the pipeline. This will pull pipeline's docker container first and build a singularity one on `~/.singularity`.
    ```bash
    $ mkdir -p ~/.singularity && cd ~/.singularity && SINGULARITY_CACHEDIR=~/.singularity SINGULARITY_PULLFOLDER=~/.singularity singularity pull --name atac-seq-pipeline-v1.2.0.simg -F docker://quay.io/encode-dcc/atac-seq-pipeline:v1.1
    ```

8. Run a pipeline for the test sample.
    ```bash
    $ INPUT=examples/local/ENCSR356KRQ_subsampled.json
    $ java -jar -Xmx1G -Dconfig.file=backends/backend.conf -Dbackend.default=sge_singularity cromwell-34.jar run atac.wdl -i ${INPUT} -o workflow_opts/sge.json
    ```

9. It will take about an hour. You will be able to find all outputs on `cromwell-executions/atac/[RANDOM_HASH_STRING]/`. See [output directory structure](output.md) for details.

10. See full specification for [input JSON file](input.md).

11. IF YOU WANT TO RUN PIPELINES WITH YOUR OWN INPUT DATA/GENOME DATABASE, PLEASE ADD THEIR DIRECTORIES TO `workflow_opts/sge.json`. For example, you have input FASTQs on `/your/input/fastqs/` and genome database installed on `/your/genome/database/` then add `/your/` to `singularity_bindpath`. You can also define multiple directories there. It's comma-separated.
    ```javascript
    {
        "default_runtime_attributes" : {
            "singularity_container" : "~/.singularity/chip-seq-pipeline-v1.1.simg",
            "singularity_bindpath" : "/your/,YOUR_OWN_DATA_DIR1,YOUR_OWN_DATA_DIR2,..."
        }
    }
    ```

## Running multiple pipelines with cromwell server mode

1. If you want to run multiple (>10) pipelines, then run a cromwell server on an interactive node. We recommend to use `screen` or `tmux` to keep your session alive and note that all running pipelines will be killed after walltime. Run a Cromwell server with the following commands.
    ```bash
    $ qlogin -h_vmem=5G -h_rt=72:00:00 # long walltime      
    $ hostname -f    # to get [CROMWELL_SVR_IP]
    ```

    For Conda users,
    ```bash
    $ source activate encode-atac-seq-pipeline
    $ _JAVA_OPTIONS="-Xmx5G" java -jar -Dconfig.file=backends/backend.conf -Dbackend.default=sge cromwell-34.jar server
    ```
    For singularity users,
    ```bash
    $ _JAVA_OPTIONS="-Xmx5G" java -jar -Dconfig.file=backends/backend.conf -Dbackend.default=sge_singularity cromwell-34.jar server
    ```

2. You can modify `backend.providers.sge.concurrent-job-limit` or `backend.providers.sge_singularity.concurrent-job-limit` in `backends/backend.conf` to increase maximum concurrent jobs. This limit is **not per sample**. It's for all sub-tasks of all submitted samples.

3. On a login node, submit jobs to the cromwell server. You will get `[WORKFLOW_ID]` as a return value. Keep these workflow IDs for monitoring pipelines and finding outputs for a specific sample later.  
    ```bash  
    $ INPUT=YOUR_INPUT.json
    $ curl -X POST --header "Accept: application/json" -v "[CROMWELL_SVR_IP]:8000/api/workflows/v1" \
      -F workflowSource=@atac.wdl \
      -F workflowInputs=@${INPUT} \
      -F workflowOptions=@workflow_opts/sge.json
    ```

  To monitor pipelines, see [cromwell server REST API description](http://cromwell.readthedocs.io/en/develop/api/RESTAPI/#cromwell-server-rest-api>) for more details. `squeue` will not give you enough information for monitoring jobs per sample.
    ```bash
    $ curl -X GET --header "Accept: application/json" -v "[CROMWELL_SVR_IP]:8000/api/workflows/v1/[WORKFLOW_ID]/status"
    ```
