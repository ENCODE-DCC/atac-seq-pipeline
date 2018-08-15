Tutorial for Stanford Sherlock 2.0 cluster
==========================================

All test samples and genome data are shared on Stanford Sherlock cluster. You don't have to download any data for testing our pipeline on it.

1. Git clone this pipeline.
    ```
      $ git clone https://github.com/ENCODE-DCC/atac-seq-pipeline
    ```

2. Move to pipeline's directory.
    ```
      $ cd atac-seq-pipeline
    ```

3. Set your partition in `workflow_opts/sherlock.json`.

    ```
      {
        "default_runtime_attributes" : {
          "slurm_partition": "YOUR_SLURM_PARTITON"
        }
      }
    ```

4. [Install Conda](https://conda.io/miniconda.html)

5. Install Conda dependencies.
    ```
      $ bash installers/install_dependencies.sh
    ```

6. Download cromwell.
    ```
      $ wget https://github.com/broadinstitute/cromwell/releases/download/34/cromwell-34.jar
      $ chmod +rx cromwell-34.jar
    ```

7. Run a pipeline for a SUBSAMPLED (1/400) paired-end sample of [ENCSR356KRQ](https://www.encodeproject.org/experiments/ENCSR356KRQ/).
    ```
      $ source activate encode-atac-seq-pipeline # IMPORTANT!
      $ INPUT=examples/sherlock/ENCSR356KRQ_subsampled.json
      $ java -jar -Dconfig.file=backends/backend.conf -Dbackend.default=slurm cromwell-34.jar run atac.wdl -i ${INPUT} -o workflow_opts/sherlock.json
    ```

8. It will take about an hour. You will be able to find all outputs on `cromwell-executions/atac/[RANDOM_HASH_STRING]/`. See [output directory structure](output.md) for details.

9. See full specification for [input JSON file](input.md).

## Extras for advanced users

1. If you want to run multiple (>10) pipelines, then run a cromwell server on an interactive node. We recommend to use `screen` or `tmux` to keep your session alive and note that all running pipelines will be killed after walltime. Run a Cromwell server with the following commands.

    ```
      $ srun -n 2 --mem 5G -t 3-0 --qos normal -p [YOUR_SLURM_PARTITION] --pty /bin/bash -i -l    # 2 CPU, 5 GB RAM and 3 day walltime
      $ hostname -f    # to get [CROMWELL_SVR_IP]

      $ source activate encode-atac-seq-pipeline
      $ _JAVA_OPTIONS="-Xmx5G" java -jar -Dconfig.file=backends/backend/conf -Dbackend.default=slurm cromwell-34.jar server
    ```

2. You can modify `backend.providers.slurm.concurrent-job-limit` in `backends/backend.conf` to increase maximum concurrent jobs. This limit is **not per sample**. It's for all sub-tasks of all submitted samples.

3. On a login node, submit jobs to the cromwell server. You will get `[WORKFLOW_ID]` as a return value. Keep these workflow IDs for monitoring pipelines and finding outputs for a specific sample later.  
    ```  
      $ INPUT=YOUR_INPUT.json
      $ curl -X POST --header "Accept: application/json" -v "[CROMWELL_SVR_IP]:8000/api/workflows/v1" \
        -F workflowSource=@atac.wdl \
        -F workflowInputs=@${INPUT} \
        -F workflowOptions=@workflow_opts/sherlock.json
    ```

  To monitor pipelines, see [cromwell server REST API description](http://cromwell.readthedocs.io/en/develop/api/RESTAPI/#cromwell-server-rest-api>) for more details. `squeue` will not give you enough information for monitoring jobs per sample.
    ```
      $ curl -X GET --header "Accept: application/json" -v "[CROMWELL_SVR_IP]:8000/api/workflows/v1/[WORKFLOW_ID]/status"
    ```
