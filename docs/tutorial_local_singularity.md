Tutorial for general UNIX computers with docker
===============================================

1. Git clone this pipeline.
    ```
      $ git clone https://github.com/ENCODE-DCC/atac-seq-pipeline
    ```

2. Move to pipeline's directory.
    ```
      $ cd atac-seq-pipeline
    ```

3. Download cromwell.
    ```
      $ wget https://github.com/broadinstitute/cromwell/releases/download/34/cromwell-34.jar
      $ chmod +rx cromwell-34.jar
    ```

4. Download a SUBSAMPLED (1/400) paired-end sample of [ENCSR356KRQ](https://www.encodeproject.org/experiments/ENCSR356KRQ/).
    ```
      $ wget https://storage.googleapis.com/encode-pipeline-test-samples/encode-atac-seq-pipeline/ENCSR356KRQ/ENCSR356KRQ_fastq_subsampled.tar
      $ tar xvf ENCSR356KRQ_fastq_subsampled.tar
    ```

5. Download pre-built genome database for hg38.
    ```
      $ wget https://storage.googleapis.com/encode-pipeline-genome-data/test_genome_database_hg38_atac.tar
      $ tar xvf test_genome_database_hg38_atac.tar
    ```

6. Pull a singularity container for the pipeline. This will pull pipeline's docker container first and build a singularity one on `~/.singularity`.
    ```
      $ SINGULARITY_PULLFOLDER=~/.singularity singularity pull docker://quay.io/encode-dcc/atac-seq-pipeline:v1.1
    ```

7. Run a pipeline for the test sample.
    ```
      $ INPUT=examples/local/ENCSR356KRQ_subsampled.json
      $ java -jar -Dconfig.file=backends/backend.conf -Dbackend.default=singularity cromwell-34.jar run atac.wdl -i ${INPUT} -o workflow_opts/singularity.json
    ```

8. It will take about an hour. You will be able to find all outputs on `cromwell-executions/atac/[RANDOM_HASH_STRING]/`. See [output directory structure](output.md) for details.

9. See full specification for [input JSON file](input.md).

## Extras for advanced users

1. To make singularity have access to input files on a different file system (e.g. NFS mounted), you need to specify `"singularity_command_options"` as "--bind ROOT_DIR_FOR_YOUR_DATA". All subdirectories/files under`ROOT_DIR_FOR_YOUR_DATA` will be bound recursively. You can also specify multiple directories (comma-separated) to be bound inside a container.

    ```
      {
          "default_runtime_attributes" : {
              "singularity_container" : "~/.singularity/atac-seq-pipeline-v1.1.img",
              "singularity_command_options" : "--bind ROOT_DIR_FOR_YOUR_DATA,DIR2,DIR3,..."
          }
      }
    ```
