ENCODE ATAC-seq pipeline
===================================================

# Directories

* `docker/` : Dockerfile
* `examples/` : input JSON examples (SE and PE)
* `genome/` : genome data TSV files
* `non-docker/` : configuration file for cluster engines (SGE and SLURM)
* `src/` : Python scripts for each task in WDL

# Usage

General usage.

```
$ java -jar [BACKEND_OPTS] cromwell-29.jar run atac.wdl -i input.json -o [WORKFLOW_OPTS]
```

* For platforms supporting docker:

     ```
     $ java -jar cromwell-29.jar run atac.wdl -i input.json -o docker.json
     ```

* For platforms without docker support:

     ```
     $ java -jar cromwell-29.jar run atac.wdl -i input.json
     ```

* For SLURM without docker:
 
     ```
     $ java -jar -Dconfig.file=backends/slurm.conf cromwell-29.jar run atac.wdl -i input.json
     ```

* For SGE without docker: Modify `backends/sge.conf` if your parallel environment (PE) is not `shm`.

     ```
     $ java -jar -Dconfig.file=backends/sge.conf cromwell-29.jar run atac.wdl -i input.json
     ```

# Input JSON

IMPORTANT NOTE on specifying input files in `input.json`:

1) For DNase-Seq, set `"bam2ta.disable_tn5_shift"=true`

2) Pipeline can start from any type of genome data (`fastq`, `bam`, `nodup_bam`, `ta` and `peak`). WDL currently does not allow optional arrays in a workflow level. Therefore, DO NOT remove input file arrays (`fastqs`, `adapters`, `bams`, `nodup_bams`, `tas`, `peaks`, `peaks_pr1`, `peaks_pr2`) from `input.json`. Also DO NOT remove `adapters` from `input.json` even if you are not starting from fastqs.
choose one of (`fastqs`, `bams`, `nodup_bams`, `tas`, `peaks`) to start with but set others as `[]`.

3) `fastqs` is an 3-dimensional array to allow merging of fastqs per replicate/endedness. Length of 3rd dimension must be 1 ([R1]) for SE and 2 ([R1, R2]) for PE.
  - 1st dimension: replicate id
  - 2nd dimension: merge id (will reduce after merging)
  - 3rd dimension: R1, R2 (single ended or paired end)

4) Other input types are just 1-dimensional arrays
  - 1st dimension: replicate id

5) Structure/dimension of `adapters` must match with that of `fastqs`. If no adapters are given then set `"atac.adapters" = []` in `input.json`. If only some adapters are known then specify them in `adapters` and leave other entries empty (`""`) while keeping the same structure/dimension as in `fastqs`. All specified/non-empty adapters will be trimmed without auto detection.

6) Set `"trim_adapter.auto_detect_adapter"=true` to automatically detect/trim adapters for empty entries in `adapters`. There will be no auto detection for non-empty entries in `adapters`. if `adapters`==`[]`, adapters will be detected/trimmed for all fastqs.

7) If starting from peaks then always specify `peaks`. Specify `peaks_pr1`, `peaks_pr2`, `peak_pooled`, `peak_ppr1` and `peak_ppr2`according to the following rules:
     ```
     if num_rep>1:
       if true_rep_only: peak_pooled, 
     else: peaks_pr1[], peaks_pr2[], peak_pooled, peak_ppr1, peak_ppr2
     else:
       if true_rep_only: not the case
     else: peaks_pr1[], peaks_pr2[]
     ```



# Dependency installation for systems without docker support

**WE DO NOT RECOMMEND RUNNIG OUR PIPELINE WITHOUT DOCKER!** Use it with caution.

1) **Our pipeline is for BASH only. Set your default shell as BASH**.

2) For Mac OS users, do not install dependencies and just install docker and use our pipeline with docker.

3) Remove any Conda (Anaconda Python and Miniconda) from your `PATH`. Pipeline will not work if other version of Conda binaries is in `PATH`.

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
   $ cd non-docker/
	 $ bash install_dependencies.sh
	 ```

9) Activate Minconda3 environment and run a pipeline
	 ```
	 $ source activate atac-seq-pipeline
	 $ java -jar [BACKEND_OPTS] $(which cromwell-29.jar) run atac.wdl -i input.json
   $ source deactivate
	 ```

# Genome data installation

**WE DO NOT RECOMMEND RUNNIG OUR PIPELINE WITH LOCALLY INSTALLED/BUILT GENOME DATA!** Use it with caution.

**We will provide an official downloader for all genome data later**. On Google Cloud, TSV files for buckets (`gs://`) will be provided. Cromwell is planning to support AWS buckets (`s3://`) too. Until then, use this installer.

Supported genomes: hg38 (ENCODE), mm10 (ENCODE), hg19 and mm9. A TSV file will be generated under `[DEST_DIR]`. Use it for `atac.genomv_tsv` value in pipeline's input JSON file

1) Do not install genome data on Stanford clusters (Sherlock and SCG4). They already have all genome data installed. Use `genome/[GENOME]_sherlock.tsv` or `genome/[GENOME]_scg4.tsv` as your TSV file.

2) For Mac OS users, if [dependency installation](#dependency-installation-for-systems-without-docker-support) does not work then post an issue on the repo.

3) Install dependencies as described in the [previous section](#dependency-installation-for-systems-without-docker-support).

4) Install genome data.
   ```
   $ cd non-docker/
   $ source activate atac-seq-pipeline
   $ bash install_genome_data.sh [GENOME] [DEST_DIR]
   $ source deactivate
   ```
