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

Choose `[BACKEND_CONF]` and `[WORKFLOW_OPT]` according to your platform and presence of `Docker`.
    ```
    $ java -jar -Dconfig.file=[BACKEND_CONF] cromwell-*.jar run atac.wdl -i input.json -o [WORKFLOW_OPT]
    ```

### Google Cloud Platform

1) Create a [Google Project](https://console.developers.google.com/project).
2) Set up a [Google Cloud Storage bucket](https://console.cloud.google.com/storage/browser) to store outputs.
3) Enable the following API's in your [API Manager](https://console.developers.google.com/apis/library).
    * Google Compute Engine
    * Google Cloud Storage
    * Genomics API
4) Set quota for `Google Compute Engine API` on https://console.cloud.google.com/iam-admin/quotas per region. Increase quota for SSD/HDD storage, number of vCPUs to process more samples faster simulateneouly.
    * CPUs
    * Persistent Disk Standard (GB)
    * Persistent Disk SSD (GB)
    * In-use IP addresses
    * Networks
5) Set `default_runtime_attributes.zones` in `workflow_opts/docker_google.json` as your preferred Google Cloud zone.
    ```
    {
      "default_runtime_attributes" : {
        ...
        "zones": "us-west1-a us-west1-b us-west1-c",
        ...
    }
    ```
6) Set `default_runtime_attributes.preemptible` as `"0"` to disable preemptible instances. Pipeline defaults to use [preemptible instances](https://cloud.google.com/compute/docs/instances/preemptible) with 10 retries. **Disabling it will cost you significantly more** but you can get your samples processed much faster and stabler. Preemptible instance is disabled in `atac.wdl` for `bowtie2` task since it can take longer than the limit (24 hours) of preemptible instances.
    ```
    {
      "default_runtime_attributes" : {
        ...
        "preemptible": "10",
        ...
    }
    ```
7) If you are already on a VM instance on your Google Project. Skip step 8 and 9.
8) Install [Google Cloud Platform SDK](https://cloud.google.com/sdk/downloads) and authenticate through it. You will be asked to enter verification keys. Get keys from the URLs they provide.
    ```
    $ gcloud auth login --no-launch-browser
    $ gcloud auth application-default login --no-launch-browser
    ```
9) Get on the Google Project.
    ```
    $ gcloud config set project [PROJ_NAME]
    ```
10) You don't have to repeat step 1-9 for next pipeline run. Credential information will be stored in `$HOME/.config/gcloud`. Go directly to step 11.
11) Run a pipeline. Use any string for `[SAMPLE_NAME]` to distinguish between multiple samples.
    ```
    $ java -jar -Dconfig.file=backends/google.conf -Dbackend.providers.JES.config.project=[PROJ_NAME] -Dbackend.providers.JES.config.root=[OUT_BUCKET]/[SAMPLE_NAME] cromwell-*.jar run atac.wdl -i input.json -o workflow_opts/docker_google.json
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

### Sun GridEngine (SGE)

Genome data have already been installed and shared on Stanford SCG4. You can skip step 3 on SCG4.
1) Set your parallel environment (`default_runtime_attributes.sge_pe`) and queue (`default_runtime_attributes.sge_queue`) in `workflow_opts/non_docker.json`. If there is no parallel environment on your SGE then ask your SGE admin to create one.
    ```
    $ qconf -spl
    ```
2) Install [dependencies](#dependency-installation).
3) Install [genome data](#genome-data-installation).
4) Run a pipeline.
    ```
    $ source activate atac-seq-pipeline
    $ java -jar -Dconfig.file=backends/sge.conf cromwell-*.jar run atac.wdl -i input.json -o workflow_opts/non_docker.json
    $ source deactivate
    ```

### SLURM

Genome data have already been installed and shared on Stanford Sherlock-2. You can skip step 3 on Sherlock-2.
1) Set your partition (`default_runtime_attributes.slurm_partition`) in `workflow_opts/non_docker.json`.
2) Install [dependencies](#dependency-installation).
3) Install [genome data](#genome-data-installation).
4) Run a pipeline.
    ```
    $ source activate atac-seq-pipeline
    $ java -jar -Dconfig.file=backends/slurm.conf cromwell-*.jar run atac.wdl -i input.json -o workflow_opts/non_docker.json
    $ source deactivate
    ```

### Kundaje lab cluster with `Docker`

Jobs will run locally without being submitted to Sun GridEngine (SGE). Genome data have already been installed and shared.
1) Run a pipeline. 
    ```
    $ java -jar -Dconfig.file=backends/default.conf cromwell-*.jar run atac.wdl -i input.json -o workflow_opts/docker.json
    ```

### Kundaje lab cluster with Sun GridEngine (SGE)

Jobs will be submitted to Sun GridEngine (SGE) and distributed to all server nodes. Genome data have already been installed and shared.
1) Install [dependencies](#dependency-installation).
2) Run a pipeline.
    ```
    $ source activate atac-seq-pipeline
    $ java -jar -Dconfig.file=backends/sge.conf cromwell-*.jar run atac.wdl -i input.json -o workflow_opts/non_docker.json
    $ source deactivate
    ```

### AWS

Not supported yet.

### GA4GH TES

Not supported yet.

### DNANexus

Not supported yet.

# Input JSON

Optional parameters and flags are marked with `?`.

1) Reference genome

    Currently supported genomes:

    * hg38: ENCODE [GRCh38_no_alt_analysis_set_GCA_000001405](https://www.encodeproject.org/files/GRCh38_no_alt_analysis_set_GCA_000001405.15/@@download/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta.gz)
    * mm10: ENCODE [mm10_no_alt_analysis_set_ENCODE](https://www.encodeproject.org/files/mm10_no_alt_analysis_set_ENCODE/@@download/mm10_no_alt_analysis_set_ENCODE.fasta.gz)
    * hg19: ENCODE [GRCh37/hg19](http://hgdownload.cse.ucsc.edu/goldenpath/hg19/encodeDCC/referenceSequences/male.hg19.fa.gz)
    * mm9: [mm9, NCBI Build 37](http://hgdownload.cse.ucsc.edu/goldenPath/mm9/bigZips/mm9.2bit)

    This TSV file has all genome specific data parameters and file path/URIs. Choose one of TSVs in `genome` directory.

    * `"atac.genome_tsv"` : TSV file path/URI.

2) Input genome data files
    Choose any genome data type you want to start with and set all others as `[]`.

    * `"atac.fastqs"` : 3-dimensional array with FASTQ file path/URI.
        - 1st dimension: replicate ID
        - 2nd dimension: merge ID (this dimension will be reduced after merging FASTQs)
        - 3rd dimension: endedness ID (0 for SE and 0,1 for PE)
    * `"atac.bams"` : Array of raw (unfiltered) BAM file path/URI.
        - 1st dimension: replicate ID
    * `"atac.nodup_bams"` : Array of filtered (deduped) BAM file path/URI.
        - 1st dimension: replicate ID
    * `"atac.tas"` : Array of TAG-ALIGN file path/URI.
        - 1st dimension: replicate ID
    * `"atac.peaks"` : Array of NARROWPEAK file path/URI.
        - 1st dimension: replicate ID
    * `"atac.peaks_pr1"` : Array of NARROWPEAK file path/URI for 1st self pseudo replicate of replicate ID.
        - 1st dimension: replicate ID
    * `"atac.peaks_pr2"` : Array of NARROWPEAK file path/URI for 2nd self pseudo replicate of replicate ID.
        - 1st dimension: replicate ID
    * `"atac.peak_ppr1"`? : NARROWPEAK file path/URI for pooled 1st pseudo replicates.
    * `"atac.peak_ppr2"`? : NARROWPEAK file path/URI for pooled 2nd pseudo replicates.
    * `"atac.peak_pooled"`? : NARROWPEAK file path/URI for pooled replicate.

    If starting from peaks then always define `"atac.peaks"`. Define `"atac.peaks_pr1"`, `"atac.peaks_pr2"`, `"atac.peak_pooled"`, `"atac.peak_ppr1"` and `"atac.peak_ppr2"` according to the following rules:    

    ```
    if num_rep>1:
        if true_rep_only: peak_pooled, 
        else: peaks_pr1[], peaks_pr2[], peak_pooled, peak_ppr1, peak_ppr2
    else:
        if true_rep_only: "not the case!"
        else: peaks_pr1[], peaks_pr2[]
    ```

3) Pipeline settings

    Input data endedness.

    * `"atac.paired_end"` : Set it as `true` if input data are paired end, otherwise `false`.

    Pipeline type (ATAC-Seq or DNase-Seq) : The only difference between two types is TN5 shifting.

    * `"atac.bam2ta.disable_tn5_shift"`? : Set it as `true` for DNase-Seq or `false` for ATAC-Seq (default).

    Other important settings.

    * `"atac.multimapping"`? : Multimapping reads.
    * `"atac.true_rep_only"`? : Set it as `true` to disable all analyses (including IDR, naive-overlap and reproducibility QC) related to pseudo replicates. This flag suppresses `"atac.enable_idr"`.

4) Adapter trimmer settings

    Structure/dimension of `"atac.adapters` must match with that of `"atac.fastqs"`. If no adapters are given then set `"atac.adapters"` as `[]` in `input.json`. If some adapters are known then define them in `"atac.adapters"` and leave other entries empty (`""`) while keeping the same structure/dimension as in `"atac.fastqs"`. All undefined/non-empty adapters will be trimmed without auto detection.
    
    * `"atac.trim_adapter.auto_detect_adapter"` : Set it as `true` to automatically detect/trim adapters for empty entries in `"atac.adapters"`. There will be no auto detection for non-empty entries it. if it is set as `[]`, adapters will be detected/trimmed for all fastqs.
    * `"atac.trim_adapter.min_trim_len"`? : Minimum trim length for `cutadapt -m`.
    * `"atac.trim_adapter.err_rate"`? : Maximum allowed adapter error rate for `cutadapt -e`.

5) Bowtie2 settings

    * `"atac.bowtie2.score_min"`? : Min. acceptable alignment score function w.r.t read length.

6) Filter/dedup (post-alignment) settings

    * `"atac.filter.dup_marker"`? : Dup marker. Choose between `picard` (default) and `sambamba`.
    * `"atac.filter.mapq_thresh"`? : Threshold for low MAPQ reads removal.
    * `"atac.filter.no_dup_removal"`? : No dup reads removal when filtering BAM.

7) BAM-2-TAGALIGN settings

    * `"atac.bam2ta.disable_tn5_shift"`? : No TN5 shifting. This is for DNase-Seq.
    * `"atac.bam2ta.regex_grep_v_ta"`? : Perl-style regular expression pattern to remove matching reads from TAGALIGN.
    * `"atac.bam2ta.subsample"`? : Number of reads to subsample TAGALIGN. Subsampled TAGALIGN will be used for all downstream analysis (MACS2, IDR, naive-overlap).

8) Cross correlation analysis settings

    * `"atac.xcor.subsample"`? : Number of reads to subsample TAGALIGN. Only one end (R1) will be used for cross correlation analysis. This will not affect downstream analysis.

9) MACS2 settings

    **DO NOT DEFINE MACS2 PARAMETERS IN `"atac.macs2"` SCOPE**. All MACS2 parameters must be defined in `"atac"` scope.

    * `"atac.cap_num_peak"`? : Cap number of raw peaks called from MACS2.
    * `"atac.pval_thresh"`? : P-value threshold.
    * `"atac.smooth_win"`? : Size of smoothing window.

10) IDR settings

    **DO NOT DEFINE IDR PARAMETERS IN `"atac.idr"` SCOPE**. All IDR parameters must be defined in `"atac"` scope.

    * `"atac.enable_idr"`? : Set it as `true` to enable IDR on raw peaks.
    * `"atac.idr_thresh"`? : IDR threshold.

11) Resources

    CPU (`cpu`), memory (`mem_mb`) settings are used for submitting jobs to cluster engines (SGE and SLURM) and Cloud platforms (Google Cloud Platform, AWS, ...). VM instance type on cloud platforms will be automatically chosen according to each task's `cpu` and `mem_mb`. Number of cores for tasks without `cpu` parameter is fixed at 1.

    * `"atac.trim_adapter.cpu"`? : Number of cores for `trim_adapter` (default: 4).
    * `"atac.bowtie2.cpu"`? : Number of cores for `bowtie2` (default: 4).
    * `"atac.filter.cpu"`? : Number of cores for `filter` (default: 2).
    * `"atac.bam2ta.cpu"`? : Number of cores for `bam2ta` (default: 2).
    * `"atac.xcor.cpu"`? : Number of cores for `xcor` (default: 2).
    * `"atac.trim_adapter.mem_mb"`? : Max. memory limit in MB for `trim_adapter` (default: 10000).
    * `"atac.bowtie2.mem_mb"`? : Max. memory limit in MB for `bowtie2` (default: 20000).
    * `"atac.filter.mem_mb"`? : Max. memory limit in MB for `filter` (default: 20000).
    * `"atac.bam2ta.mem_mb"`? : Max. memory limit in MB for `bam2ta` (default: 10000).
    * `"atac.xcor.mem_mb"`? : Max. memory limit in MB for `xcor` (default: 10000).
    * `"atac.macs2_mem_mb"`? : Max. memory limit in MB for `macs2` (default: 16000).

    Disks (`disks`) is used for Cloud platforms (Google Cloud Platforms, AWS, ...).

    * `"atac.trim_adapter.disks"`? : Disks for `trim_adapter` (default: "local-disk 100 HDD").
    * `"atac.bowtie2.disks"`? : Disks for `bowtie2` (default: "local-disk 100 HDD").
    * `"atac.filter.disks"`? : Disks for `filter` (default: "local-disk 100 HDD").
    * `"atac.bam2ta.disks"`? : Disks for `bam2ta` (default: "local-disk 100 HDD").
    * `"atac.xcor.disks"`? : Disks for `xcor` (default: "local-disk 100 HDD").
    * `"atac.macs2_disks"`? : Disks for `macs2` (default: "local-disk 100 HDD").

    Walltime (`time`) settings (for SGE and SLURM only).

    * `"atac.trim_adapter.time_hr"`? : Walltime for `trim_adapter` (default: 24).
    * `"atac.bowtie2.time_hr"`? : Walltime for `bowtie2` (default: 48).
    * `"atac.filter.time_hr"`? : Walltime for `filter` (default: 24).
    * `"atac.bam2ta.time_hr"`? : Walltime for `bam2ta` (default: 6).
    * `"atac.xcor.time_hr"`? : Walltime for `xcor` (default: 6).
    * `"atac.macs2_time_hr"`? : Walltime for `macs2` (default: 24).

# Dependency installation

**WE DO NOT RECOMMEND RUNNING OUR PIPELINE WITHOUT `DOCKER`!** Use it with caution.
1) **Our pipeline is for BASH only. Set your default shell as BASH**.
2) For Mac OSX users, do not install dependencies and just install `Docker` and use our pipeline with it.
3) Remove any Conda (Anaconda Python and Miniconda) from your `PATH`. **PIPELINE WILL NOT WORK IF YOU HAVE OTHER VERSION OF CONDA BINARIES IN `PATH`**.
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

Supported genomes:

  * hg38: ENCODE [GRCh38_no_alt_analysis_set_GCA_000001405](https://www.encodeproject.org/files/GRCh38_no_alt_analysis_set_GCA_000001405.15/@@download/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta.gz)
  * mm10: ENCODE [mm10_no_alt_analysis_set_ENCODE](https://www.encodeproject.org/files/mm10_no_alt_analysis_set_ENCODE/@@download/mm10_no_alt_analysis_set_ENCODE.fasta.gz)
  * hg19: ENCODE [GRCh37/hg19](http://hgdownload.cse.ucsc.edu/goldenpath/hg19/encodeDCC/referenceSequences/male.hg19.fa.gz)
  * mm9: [mm9, NCBI Build 37](http://hgdownload.cse.ucsc.edu/goldenPath/mm9/bigZips/mm9.2bit)

A TSV file will be generated under `[DEST_DIR]`. Use it for `atac.genomv_tsv` value in pipeline's input JSON file.

1) Do not install genome data on Stanford clusters (Sherlock-2 and SCG4). They already have all genome data installed. Use `genome/[GENOME]_sherlock.tsv` or `genome/[GENOME]_scg4.tsv` as your TSV file.
2) For Mac OSX users, if [dependency installation](#dependency-installation) does not work then post an issue on the repo.
3) Install [dependencies](#dependency-installation) first.
4) Install genome data.
   ```
   $ cd installers/
   $ source activate atac-seq-pipeline
   $ bash install_genome_data.sh [GENOME] [DEST_DIR]
   $ source deactivate
   ```

# Dev

### Example command lines

```
# SE
java -jar -Dconfig.file=backends/google.conf -Dbackend.providers.JES.config.project=atac-seq-pipeline -Dbackend.providers.JES.config.root="gs://atac-seq-pipeline-workflows/ENCSR889WQX" cromwell-29.jar run atac.wdl -i examples/ENCSR889WQX_google.json -o workflow_opts/docker_google.json -m output_ENCSR889WQX.json

# PE
java -jar -Dconfig.file=backends/google.conf -Dbackend.providers.JES.config.project=atac-seq-pipeline -Dbackend.providers.JES.config.root="gs://atac-seq-pipeline-workflows/ENCSR356KRQ" cromwell-29.jar run atac.wdl -i examples/ENCSR356KRQ_google.json -o workflow_opts/docker_google.json -m output_ENCSR356KRQ.json
```

### Docker build

```
docker build -f docker_image/Dockerfile .
```