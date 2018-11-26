Tutorial for general UNIX computers without docker
==================================================

1. Git clone this pipeline and move into it.
    ```
      $ git clone https://github.com/ENCODE-DCC/atac-seq-pipeline
      $ cd atac-seq-pipeline
    ```

2. Download [cromwell](https://github.com/broadinstitute/cromwell).
    ```
      $ wget https://github.com/broadinstitute/cromwell/releases/download/34/cromwell-34.jar
      $ chmod +rx cromwell-34.jar
    ```

3. Download a SUBSAMPLED (1/400) paired-end sample of [ENCSR356KRQ](https://www.encodeproject.org/experiments/ENCSR356KRQ/).
    ```
      $ wget https://storage.googleapis.com/encode-pipeline-test-samples/encode-atac-seq-pipeline/ENCSR356KRQ/ENCSR356KRQ_fastq_subsampled.tar
      $ tar xvf ENCSR356KRQ_fastq_subsampled.tar
    ```

4. Download pre-built genome database for hg38.
    ```
      $ wget https://storage.googleapis.com/encode-pipeline-genome-data/test_genome_database_hg38_atac.tar
      $ tar xvf test_genome_database_hg38_atac.tar
    ```

5. [Install Conda](https://conda.io/miniconda.html)

6. Install Conda dependencies.
    ```
      $ bash conda/uninstall_dependencies.sh  # to remove any existing pipeline env
      $ bash conda/install_dependencies.sh
    ```
    
7. Run a pipeline for the test sample.
    ```
      $ source activate encode-atac-seq-pipeline # IMPORTANT!
      $ INPUT=examples/local/ENCSR356KRQ_subsampled.json
      $ java -jar -Dconfig.file=backends/backend.conf cromwell-34.jar run atac.wdl -i ${INPUT}
    ```

8. It will take about an hour. You will be able to find all outputs on `cromwell-executions/atac/[RANDOM_HASH_STRING]/`. See [output directory structure](output.md) for details.

9. See full specification for [input JSON file](input.md).
