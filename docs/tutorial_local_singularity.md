Tutorial for general UNIX computers with singularity
====================================================

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

5. Pull a singularity container for the pipeline. This will pull pipeline's docker container first and build a singularity one on `~/.singularity`.
    ```
      $ SINGULARITY_PULLFOLDER=~/.singularity singularity pull docker://quay.io/encode-dcc/atac-seq-pipeline:v1.1
    ```

6. **DO NOT SKIP THIS STEP OR PIPEPILE WILL FAIL.** Look for `--bind` in `singularity_command_options`. If you want to run pipelines with your own input data and genome database then you may need to add their directdories to `--bind`. It's comma-separated and all sub-directories under these directories will be bound recursively to singularity. Therefore, pick a root directory (e.g. your scratch folder `/scratch/users/$USER`) for your data.
    ```
      {
          "default_runtime_attributes" : {
              "singularity_container" : "~/.singularity/atac-seq-pipeline-v1.1.simg",
              "singularity_command_options" : "--bind YOUR_WORK_DIR"
          }
      }
    ```

7. Run a pipeline for the test sample.
    ```
      $ INPUT=examples/local/ENCSR356KRQ_subsampled.json
      $ java -jar -Dconfig.file=backends/backend.conf -Dbackend.default=singularity cromwell-34.jar run atac.wdl -i ${INPUT} -o workflow_opts/singularity.json
    ```

8. It will take about an hour. You will be able to find all outputs on `cromwell-executions/atac/[RANDOM_HASH_STRING]/`. See [output directory structure](output.md) for details.

9. See full specification for [input JSON file](input.md).

10. Please read through the next section. It's VERY IMPORTANT.

## Binding directories for singularity

1. Singularity automatically binds two directories (your `$HOME` and your current working directory `$PWD`) RECURSIVELY. So it has access to any physical/hard-linked data under these directories. See item 2) for soft-linked (linked with `ln -s`) data.

2. However, if your input data live on a different location then you need to manually bind them for singularity. For example, you git cloned our pipeline on your `$HOME/code/atac-seq-pipeline` and run pipelines. You want to use `/scratch/fastq.gz` as an input (which is outside of your `$HOME` and `$PWD`) then singularity will not be able find it. Even if your input is `$HOME/fastq.gz` but it's a soft link to `/scratch2/any/sub-dir/fastq.gz` then singularity cannot find it either.

3. Therefore, you need to modify `workflow_opts/singularity.json` to specify such bindings (`"singularity_command_options"` as `"--bind ROOT_DIR_FOR_YOUR_DATA"`). All subdirectories/files under `ROOT_DIR_FOR_YOUR_DATA` will be bound RECURSIVELY. You can also specify multiple directories (comma-separated). For the example of the item 2) see the following JSON:
    ```
      {
          "default_runtime_attributes" : {
              "singularity_container" : "~/.singularity/atac-seq-pipeline-v1.1.img",
              "singularity_command_options" : "--bind /scratch,/scratch2"
          }
      }
    ```