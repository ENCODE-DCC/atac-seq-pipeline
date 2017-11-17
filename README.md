ENCODE ATAC-seq pipeline
===================================================

# Directories

* `backends/` : Backend configuration files (`.conf`)
* `workflow_opts/` : Workflow option files (`.json`)
* `examples/` : input JSON examples (SE and PE)
* `genome/` : genome data TSV files
* `src/` : Python script for each task in WDL
* `installers/` : dependency/genome data installers for Local, SGE and SLURM
* `docker_image/` : Dockerfile

# Usage

Choose `[BACKEND_CONF]` and `[WORKFLOW_OPT]` according to your platform and presence of `docker`.
  ```
  $ java -jar -Dconfig.file=[BACKEND_CONF] cromwell-*.jar run atac.wdl -i input.json -o [WORKFLOW_OPT]
  ```

### Local computer with `Docker`

1) Install [genome data](#genome-data-installation).
2) Run a pipeline.
  ```
  $ java -jar -Dconfig.file=backends/default.conf cromwell-*.jar run atac.wdl -i input.json -o workflow_opts/docker.json
  ```

### Local computer without `Docker`

1) Install [dependencies](#dependency-installation).
2) Install [genome data](#genome-data-installation).
3) Run a pipeline.
  ```
  $ source activate atac-seq-pipeline
  $ java -jar -Dconfig.file=backends/default.conf cromwell-*.jar run atac.wdl -i input.json -o workflow_opts/non_docker.json
  $ source deactivate
  ```

### Kundaje Lab cluster with `Docker`

Jobs will run locally without being submitted to Sun GridEngine (SGE). Genome data have already been installed and shared.
1) Run a pipeline. 
  ```
  $ java -jar -Dconfig.file=backends/default.conf cromwell-*.jar run atac.wdl -i input.json -o workflow_opts/docker.json
  ```

### Kundaje Lab cluster with Sun GridEngine (SGE)

Jobs will be submitted to Sun GridEngine (SGE) and distributed to all server nodes. Genome data have already been installed and shared.
1) Install [dependencies](#dependency-installation).
2) Run a pipeline.
  ```
  $ source activate atac-seq-pipeline
  $ java -jar -Dconfig.file=backends/sge.conf cromwell-*.jar run atac.wdl -i input.json -o workflow_opts/non_docker.json
  $ source deactivate
  ```

### Sun GridEngine (SGE)

1) Check if your SGE has a parallel environment named `shm` by `$ qconf -spl`. If it does not have a PE `shm` then ask your SGE admin to create one or change the name of PE in `backends/sge.conf`.
2) Install [dependencies](#dependency-installation).
3) Install [genome data](#genome-data-installation).
4) Run a pipeline.
  ```
  $ source activate atac-seq-pipeline
  $ java -jar -Dconfig.file=backends/sge.conf cromwell-*.jar run atac.wdl -i input.json -o workflow_opts/non_docker.json
  $ source deactivate
  ```

### SLURM

1) Install [dependencies](#dependency-installation).
2) Install [genome data](#genome-data-installation).
3) Run a pipeline.
  ```
  $ source activate atac-seq-pipeline
  $ java -jar -Dconfig.file=backends/slurm.conf cromwell-*.jar run atac.wdl -i input.json -o workflow_opts/non_docker.json
  $ source deactivate
  ```

### Sun GridEngine (SGE) on Stanford SCG4 cluster

Genome data have already been installed and shared.
1) Install [dependencies](#dependency-installation).
2) Run a pipeline.
  ```
  $ source activate atac-seq-pipeline
  $ java -jar -Dconfig.file=backends/sge.conf cromwell-*.jar run atac.wdl -i input.json -o workflow_opts/non_docker.json
  $ source deactivate
  ```

### SLURM on Stanford Sherlock-2 cluster

Genome data have already been installed and shared.
1) Install [dependencies](#dependency-installation).
2) Run a pipeline.
  ```
  $ source activate atac-seq-pipeline
  $ java -jar -Dconfig.file=backends/slurm.conf cromwell-*.jar run atac.wdl -i input.json -o workflow_opts/non_docker.json
  $ source deactivate
  ```

### Google Cloud Platform

1) Create a [Google Project](https://console.developers.google.com/project).
2) Set up a [Google Cloud Storage bucket](https://console.cloud.google.com/storage/browser) to store outputs.
3) Enable the following API's in your [API Manager](https://console.developers.google.com/apis/library).
  * Google Compute Engine
  * Google Cloud Storage
  * Genomics API
4) If you are already on a VM instance on your Google Project. Skip step 5 and 6.
5) Install [Google Cloud Platform SDK](https://cloud.google.com/sdk/downloads) and authenticate through it. You will be asked to enter verification keys. Get keys from the URLs they provide.
  ```
  $ gcloud auth login --no-launch-browser
  $ gcloud auth application-default login --no-launch-browser
  ```
6) Get on the Google Project.
  ```
  $ gcloud config set project [PROJ_NAME]
  ```
7) You don't have to repeat step 1-6 for next pipeline run. Credential information will be stored in `$HOME/.config/gcloud`. Go directly to step 8.
8) Run a pipeline. Use any string for `[SAMPLE_NAME]` to distinguish between multiple samples.
  ```
  $ java -jar -Dconfig.file=backends/google.conf -Dbackend.providers.JES.config.project=[PROJ_NAME] -Dbackend.providers.JES.config.root=[OUT_BUCKET]/[SAMPLE_NAME] cromwell-*.jar run atac.wdl -i input.json -o workflow_opts/docker_google.json
  ```

### AWS

Not supported yet.

### GA4GH TES

Not supported yet.

### DNANexus

Not supported yet.

# Input JSON

* For DNase-Seq, set `"bam2ta.disable_tn5_shift"=true`
* Pipeline can start from any type of genome data (`fastq`, `bam`, `nodup_bam`, `ta` and `peak`). WDL currently does not allow optional arrays in a workflow level. Therefore, DO NOT remove input file arrays (`fastqs`, `adapters`, `bams`, `nodup_bams`, `tas`, `peaks`, `peaks_pr1`, `peaks_pr2`) from `input.json`. Also DO NOT remove `adapters` from `input.json` even if you are not starting from fastqs.
choose one of (`fastqs`, `bams`, `nodup_bams`, `tas`, `peaks`) to start with but set others as `[]`.
* `fastqs` is an 3-dimensional array to allow merging of fastqs per replicate/endedness. Length of 3rd dimension must be 1 ([R1]) for SE and 2 ([R1, R2]) for PE.
  - 1st dimension: replicate id
  - 2nd dimension: merge id (will reduce after merging)
  - 3rd dimension: R1, R2 (single ended or paired end)
* Other input types are just 1-dimensional arrays
  - 1st dimension: replicate id
* Structure/dimension of `adapters` must match with that of `fastqs`. If no adapters are given then set `"atac.adapters" = []` in `input.json`. If only some adapters are known then specify them in `adapters` and leave other entries empty (`""`) while keeping the same structure/dimension as in `fastqs`. All specified/non-empty adapters will be trimmed without auto detection.
* Set `"trim_adapter.auto_detect_adapter"=true` to automatically detect/trim adapters for empty entries in `adapters`. There will be no auto detection for non-empty entries in `adapters`. if `adapters`==`[]`, adapters will be detected/trimmed for all fastqs.
* If starting from peaks then always specify `peaks`. Specify `peaks_pr1`, `peaks_pr2`, `peak_pooled`, `peak_ppr1` and `peak_ppr2`according to the following rules:
  ```
  if num_rep>1:
      if true_rep_only: peak_pooled, 
      else: peaks_pr1[], peaks_pr2[], peak_pooled, peak_ppr1, peak_ppr2
  else:
      if true_rep_only: "not the case!"
      else: peaks_pr1[], peaks_pr2[]
  ```

# Dependency installation

**WE DO NOT RECOMMEND RUNNING OUR PIPELINE WITHOUT `DOCKER`!** Use it with caution.
1) **Our pipeline is for BASH only. Set your default shell as BASH**.
2) For Mac OSX users, do not install dependencies and just install `Docker` and use our pipeline with it.
3) Remove any Conda (Anaconda Python and Miniconda) from your `PATH`.**Pipeline will not work if you have other version of Conda binaries in `PATH`**.
4) Install Miniconda3 for 64-bit Linux on your system. Miniconda2 will not work. If your system is 32-bit Linux then try with `x86_32`.
   ```
   $ wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
   $ bash Miniconda3-latest-Linux-x86_64.sh -b -p [MINICONDA3_INSTALL_DIR]
   ```
5) Add `PATH` for our pipeline Python scripts and Miniconda3 to one of your bash startup scripts (`$HOME/.bashrc`, `$HOME/.bash_profile`...). 
   ```
   unset PYTHONPATH
   export PATH=[ATAC_SEQ_PIPELINE_DIR]/src:$PATH
   export PATH=[MINICONDA3_INSTALL_DIR]/bin:$PATH
   ```
6) Re-login.
7) Make sure that conda correctly points to `[MINICONDA3_INSTALL_DIR]/bin/conda`.
   ```
   $ which conda
   ```
8) Install dependencies on Minconda3 environment. Java 8 JDK and Cromwell-29 are included in the installation.
   ```
   $ cd installers/
   $ bash install_dependencies.sh
   ```
9) **ACTIVATE MINICONDA3 ENVIRONMENT** and run a pipeline.
   ```
   $ source activate atac-seq-pipeline
   $ java -jar -Dconfig.file=[BACKEND_OPT] $(which cromwell-29.jar) run atac.wdl -i input.json -o [WORKFLOW_OPT]
   $ source deactivate
   ```

# Genome data installation

**WE DO NOT RECOMMEND RUNNING OUR PIPELINE WITH LOCALLY INSTALLED/BUILT GENOME DATA!** Use it with caution. **We will provide an official downloader for all genome data later**. Cromwell is planning to support AWS buckets (`s3://`). Until then, use this installer.
**On Google Cloud TSV** files are already installed and shared on a bucket `gs://atac-seq-pipeline-genome-data`.

Supported genomes: hg38 (from ENCODE), mm10 (from ENCODE), hg19 and mm9. A TSV file will be generated under `[DEST_DIR]`. Use it for `atac.genomv_tsv` value in pipeline's input JSON file

1) Do not install genome data on Stanford clusters (Sherlock-2 and SCG4). They already have all genome data installed. Use `genome/[GENOME]_sherlock.tsv` or `genome/[GENOME]_scg4.tsv` as your TSV file.
2) For Mac OSX users, if [dependency installation](#dependency-installation) does not work then post an issue on the repo.
3) Install [dependencies](#dependency-installation) first .
4) Install genome data.
   ```
   $ cd installers/
   $ source activate atac-seq-pipeline
   $ bash install_genome_data.sh [GENOME] [DEST_DIR]
   $ source deactivate
   ```

# Examples

```
java -jar -Dconfig.file=backends/google.conf -Dbackend.providers.JES.config.project=atac-seq-pipeline -Dbackend.providers.JES.config.root="gs://atac-seq-pipeline-workflows/ENCSR889WQX" cromwell-29.jar run atac.wdl -i examples/ENCSR889WQX_google.json -o workflow_opts/docker_google.json
```