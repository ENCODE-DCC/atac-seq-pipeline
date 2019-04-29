# Tutorial for general UNIX computers with singularity

1. Download [cromwell](https://github.com/broadinstitute/cromwell).
    ```bash
    $ cd
    $ wget https://github.com/broadinstitute/cromwell/releases/download/34/cromwell-34.jar
    $ chmod +rx cromwell-34.jar
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

5. CHECK YOUR SINGULARITY VERSION FIRST AND UPGRADE IT TO A VERSION `>=2.5.2` OR PIPELINE WILL NOT WORK CORRECTLY.
    ```bash
    $ singularity --version
    ```

6. Pull a singularity container for the pipeline. This will pull pipeline's docker container first and build a singularity one on `~/.singularity`.
    ```bash
    $ mkdir -p ~/.singularity && cd ~/.singularity && SINGULARITY_CACHEDIR=~/.singularity SINGULARITY_PULLFOLDER=~/.singularity singularity pull --name atac-seq-pipeline-v1.3.0.simg -F docker://quay.io/encode-dcc/atac-seq-pipeline:v1.3.0
    ```

7. Run a pipeline for the test sample.
    ```bash
    $ INPUT=examples/local/ENCSR356KRQ_subsampled.json
    $ PIPELINE_METADATA=metadata.json
    $ java -jar -Xmx1G -Dconfig.file=backends/backend.conf -Dbackend.default=singularity cromwell-34.jar run atac.wdl -i ${INPUT} -o workflow_opts/singularity.json -m ${PIPELINE_METADATA}
    ```

8. It will take about an hour. You will be able to find all outputs on `cromwell-executions/atac/[RANDOM_HASH_STRING]/`. See [output directory structure](output.md) for details.

9. See full specification for [input JSON file](input.md).

10. You can resume a failed pipeline from where it left off by using `PIPELINE_METADATA`(`metadata.json`) file. This file is created for each pipeline run. See [here](../utils/resumer/README.md) for details. Once you get a new input JSON file from the resumer, use it `INPUT=resume.[FAILED_WORKFLOW_ID].json` instead of `INPUT=examples/local/ENCSR356KRQ_subsampled.json`.

11. IF YOU WANT TO RUN PIPELINES WITH YOUR OWN INPUT DATA/GENOME DATABASE, PLEASE ADD THEIR DIRECTORIES TO `workflow_opts/singularity.json`. For example, you have input FASTQs on `/your/input/fastqs/` and genome database installed on `/your/genome/database/` then add `/your/` to `singularity_bindpath`. You can also define multiple directories there. It's comma-separated.
    ```javascript
    {
        "default_runtime_attributes" : {
            "singularity_container" : "~/.singularity/chip-seq-pipeline-v1.3.0.simg",
            "singularity_bindpath" : "/your/,YOUR_OWN_DATA_DIR1,YOUR_OWN_DATA_DIR1,..."
        }
    }
    ```