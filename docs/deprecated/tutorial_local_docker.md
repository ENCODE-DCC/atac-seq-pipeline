# Tutorial for general UNIX computers with docker

1. Download [cromwell](https://github.com/broadinstitute/cromwell).
    ```bash
    $ cd
    $ wget https://github.com/broadinstitute/cromwell/releases/download/38/cromwell-38.jar
    $ chmod +rx cromwell-38.jar
    ```

2. Git clone this pipeline and move into it.
    ```bash
    $ cd
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
    
5. Run a pipeline for the test sample.
    ```bash
    $ INPUT=dev/examples/local/ENCSR356KRQ_subsampled.json
    $ PIPELINE_METADATA=metadata.json
    $ java -jar -Dconfig.file=dev/backends/backend.conf cromwell-38.jar run atac.wdl -i ${INPUT} -o dev/workflow_opts/docker.json -m ${PIPELINE_METADATA}
    ```

6. It will take about an hour. You will be able to find all outputs on `cromwell-executions/atac/[RANDOM_HASH_STRING]/`. See [output directory structure](output.md) for details.

7. See full specification for [input JSON file](input.md).

8. You can resume a failed pipeline from where it left off by using `PIPELINE_METADATA`(`metadata.json`) file. This file is created for each pipeline run. See [here](../utils/resumer/README.md) for details. Once you get a new input JSON file from the resumer, use it `INPUT=resume.[FAILED_WORKFLOW_ID].json` instead of `INPUT=dev/examples/local/ENCSR356KRQ_subsampled.json`.