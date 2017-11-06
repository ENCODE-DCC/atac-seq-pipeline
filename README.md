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

Take a careful look at input definition in [`atac.wdl`](atac.wdl). Read through all comments at the top of the code.

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
