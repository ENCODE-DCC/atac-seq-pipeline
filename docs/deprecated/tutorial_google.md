# Tutorial for Google Cloud Platform

All test samples and genome data are shared on our public Google Cloud buckets. You don't have to download any data for testing our pipeline on Google Cloud.

1. Sign up for a Google account.
2. Go to [Google Project](https://console.developers.google.com/project) page and click "SIGN UP FOR FREE TRIAL" on the top left and agree to terms.
3. Set up a payment method and click "START MY FREE TRIAL".
4. Create a [Google Project](https://console.developers.google.com/project) `[YOUR_PROJECT_NAME]` and choose it on the top of the page.
5. Create a [Google Cloud Storage bucket](https://console.cloud.google.com/storage/browser) `gs://[YOUR_BUCKET_NAME]` by clicking on a button "CREATE BUCKET" and create it to store pipeline outputs.
6. Find and enable following APIs in your [API Manager](https://console.developers.google.com/apis/library). Click a back button on your web brower after enabling each.
    * Compute Engine API
    * Google Cloud Storage (DO NOT click on "Create credentials")
    * Google Cloud Storage JSON API
    * Genomics API

7. Install [Google Cloud Platform SDK](https://cloud.google.com/sdk/downloads) and authenticate through it. You will be asked to enter verification keys. Get keys from the URLs they provide.
    ```bash
    $ gcloud auth login --no-launch-browser
    $ gcloud auth application-default login --no-launch-browser
    ```

8. If you see permission errors at runtime, then unset environment variable `GOOGLE_APPLICATION_CREDENTIALS` or add it to your BASH startup scripts (`$HOME/.bashrc` or `$HOME/.bash_profile`).
    ```bash
      unset GOOGLE_APPLICATION_CREDENTIALS
    ```

7. Set your default Google Cloud Project. Pipeline will provision instances on this project.
    ```bash
    $ gcloud config set project [YOUR_PROJECT_NAME]
    ```

8. Download [cromwell](https://github.com/broadinstitute/cromwell).
    ```bash
    $ cd
    $ wget https://github.com/broadinstitute/cromwell/releases/download/38/cromwell-38.jar
    $ chmod +rx cromwell-38.jar
    ```

9. Git clone this pipeline and move into it.
    ```bash
    $ cd
    $ git clone https://github.com/ENCODE-DCC/atac-seq-pipeline
    $ cd atac-seq-pipeline
    ```

10. Run a pipeline for a SUBSAMPLED (1/400) paired-end sample of [ENCSR356KRQ](https://www.encodeproject.org/experiments/ENCSR356KRQ/).
    ```bash
    $ PROJECT=[YOUR_PROJECT_NAME]
    $ BUCKET=gs://[YOUR_BUCKET_NAME]/ENCSR356KRQ_subsampled
    $ INPUT=dev/examples/google/ENCSR356KRQ_subsampled.json
    $ PIPELINE_METADATA=metadata.json

    $ java -jar -Dconfig.file=dev/backends/backend.conf -Dbackend.default=google -Dbackend.providers.google.config.project=${PROJECT} -Dbackend.providers.google.config.root=${BUCKET} cromwell-38.jar run atac.wdl -i ${INPUT} -o dev/workflow_opts/docker.json -m ${PIPELINE_METADATA}
    ```

11. It will take about an hour. You will be able to find all outputs on your Google Cloud bucket. Final QC report/JSON will be written on `gs://[YOUR_BUCKET_NAME]/ENCSR356KRQ_subsampled/atac/[SOME_HASH_STRING]/call-qc_report/execution/glob*/qc.html` or `qc.json`. See [output directory structure](output.md) for details.

12. See full specification for [input JSON file](input.md).

13. You can resume a failed pipeline from where it left off by using `PIPELINE_METADATA`(`metadata.json`) file. This file is created for each pipeline run. See [here](../utils/resumer/README.md) for details. Once you get a new input JSON file from the resumer, use it `INPUT=resume.[FAILED_WORKFLOW_ID].json` instead of `INPUT=dev/examples/google/ENCSR356KRQ_subsampled.json`.

## Extras for advanced users

1. Set quota for [Google Compute Engine API](https://console.cloud.google.com/iam-admin/quotas) per region. Increase quota for SSD/HDD storage, number of vCPUs to process more samples faster simulateneouly.
    * CPUs
    * Persistent Disk Standard (GB)
    * Persistent Disk SSD (GB)
    * In-use IP addresses
    * Networks

2. Set `default_runtime_attributes.zones` in `dev/workflow_opts/docker.json` as your preferred Google Cloud zone.
    ```javascript
    {
      "default_runtime_attributes" : {
        ...
        "zones": "us-west1-a us-west1-b us-west1-c",
        ...
    }
    ```

3. Set `default_runtime_attributes.preemptible` as `"0"` to disable preemptible instances. This value means a number of retrial for failures in a preemtible instance. Pipeline defaults not to use [preemptible instances](https://cloud.google.com/compute/docs/instances/preemptible). If all retrial fails then the instance will be upgraded to a regular one. **Disabling preemtible instances will cost you significantly more** but you can get your samples processed much faster and stabler. Preemptible instance is disabled by default. Some hard tasks like `bowtie2`, `bwa` and `spp` will not be executed on preemtible instances since they can take longer than the limit (24 hours) of preemptible instances.
    ```javascript
    {
      "default_runtime_attributes" : {
        ...
        "preemptible": "0",
        ...
    }
    ```