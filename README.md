# ENCODE ATAC-seq pipeline

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.156534.svg)](https://doi.org/10.5281/zenodo.156534)[![CircleCI](https://circleci.com/gh/ENCODE-DCC/atac-seq-pipeline/tree/master.svg?style=svg)](https://circleci.com/gh/ENCODE-DCC/atac-seq-pipeline/tree/master)


## Updated genome TSV files (v3 -> v4)



## Download new Caper>=2.1

New Caper is out. You need to update your Caper to work with the latest ENCODE ATAC-seq pipeline.
```bash
$ pip install caper --upgrade
```

## Local/HPC users and new Caper>=2.1

There are tons of changes for local/HPC backends: `local`, `slurm`, `sge`, `pbs` and `lsf`(added). Make a backup of your current Caper configuration file `~/.caper/default.conf` and run `caper init`. Local/HPC users need to reset/initialize Caper's configuration file according to your chosen backend. Edit the configuration file and follow instructions in there.
```bash
$ cd ~/.caper
$ cp default.conf default.conf.bak
$ caper init [YOUR_BACKEND]
```

In order to run a pipeline, you need to add one of the following flags to specify the environment to run each task within. i.e. `--conda`, `--singularity` and `--docker`. These flags are not required for cloud backend users (`aws` and `gcp`).
```bash
# for example
$ caper run ... --singularity
```

For Conda users, **RE-INSTALL PIPELINE'S CONDA ENVIRONMENT AND DO NOT ACTIVATE CONDA ENVIRONMENT BEFORE RUNNING PIPELINES**. Caper will internally call `conda run -n ENV_NAME CROMWELL_JOB_SCRIPT`. Just make sure that pipeline's new Conda environments are correctly installed.
```bash
$ scripts/uninstall_conda_env.sh
$ scripts/install_conda_env.sh
```

## Introduction

This pipeline is designed for automated end-to-end quality control and processing of ATAC-seq and DNase-seq data. The pipeline can be run on compute clusters with job submission engines as well as on stand alone machines. It inherently makes uses of parallelized/distributed computing. Pipeline installation is also easy as most dependencies are automatically installed. The pipeline can be run end-to-end, starting from raw FASTQ files all the way to peak calling and signal track generation using a single caper submit command. One can also start the pipeline from intermediate stages (for example, using alignment files as input). The pipeline supports both single-end and paired-end data as well as replicated or non-replicated datasets. The outputs produced by the pipeline include 1) formatted HTML reports that include quality control measures specifically designed for ATAC-seq and DNase-seq data, 2) analysis of reproducibility, 3) stringent and relaxed thresholding of peaks, 4) fold-enrichment and pvalue signal tracks. The pipeline also supports detailed error reporting and allows for easy resumption of interrupted runs. It has been tested on some human, mouse and yeast ATAC-seq datasets as well as on human and mouse DNase-seq datasets.

The ATAC-seq pipeline protocol specification is [here](https://docs.google.com/document/d/1f0Cm4vRyDQDu0bMehHD7P7KOMxTOP-HiNoIvL1VcBt8/edit?usp=sharing). Some parts of the ATAC-seq pipeline were developed in collaboration with Jason Buenrostro, Alicia Schep and Will Greenleaf at Stanford.

### Features

* **Portability**: The pipeline run can be performed across different cloud platforms such as Google, AWS and DNAnexus, as well as on cluster engines such as SLURM, SGE and PBS.
* **User-friendly HTML report**: In addition to the standard outputs, the pipeline generates an HTML report that consists of a tabular representation of quality metrics including alignment/peak statistics and FRiP along with many useful plots (IDR/TSS enrichment). An example of the [HTML report](https://storage.googleapis.com/encode-pipeline-test-samples/encode-atac-seq-pipeline/ENCSR889WQX/example_output/qc.html). The [json file](https://storage.googleapis.com/encode-pipeline-test-samples/encode-atac-seq-pipeline/ENCSR889WQX/example_output/qc.json) used in generating this report.
* **Supported genomes**: Pipeline needs genome specific data such as aligner indices, chromosome sizes file and blacklist. We provide a genome database downloader/builder for hg38, hg19, mm10, mm9. You can also use this [builder](docs/build_genome_database.md) to build genome database from FASTA for your custom genome.

## Installation

1) Make sure that you have Python>=3.6. Caper does not work with Python2. Install Caper and check its version >=2.0.
	```bash
	$ python --version
	$ pip install caper
	```
2) Make a backup of your Caper configuration file `~/.caper/default.conf` if you are upgrading from old Caper(<2.0.0). Reset/initialize Caper's configuration file. Read Caper's [README](https://github.com/ENCODE-DCC/caper/blob/master/README.md) carefully to choose a backend for your system. Follow the instruction in the configuration file.
	```bash
	# make a backup of ~/.caper/default.conf if you already have it
	$ caper init [YOUR_BACKEND]

	# then edit ~/.caper/default.conf
	$ vi ~/.caper/default.conf
	```

3) Git clone this pipeline.
	> **IMPORTANT**: use `~/atac-seq-pipeline/atac.wdl` as `[WDL]` in Caper's documentation.

	```bash
	$ cd
	$ git clone https://github.com/ENCODE-DCC/atac-seq-pipeline
	```

4) (Optional for Conda users) Install pipeline's Conda environments if you don't have Singularity or Docker installed on your system. We recommend to use Singularity instead of Conda. If you don't have Conda on your system, install [Miniconda3](https://docs.conda.io/en/latest/miniconda.html).
	```bash
	$ cd atac-seq-pipeline
	# uninstall old environments (<2.0.0)
	$ bash scripts/uninstall_conda_env.sh
	$ bash scripts/install_conda_env.sh
	```


## Input JSON file specification

> **IMPORTANT**: DO NOT BLINDLY USE A TEMPLATE/EXAMPLE INPUT JSON. READ THROUGH THE FOLLOWING GUIDE TO MAKE A CORRECT INPUT JSON FILE. ESPECIALLY FOR AUTODETECTING/DEFINING ADAPTERS.

An input JSON file specifies all the input parameters and files that are necessary for successfully running this pipeline. This includes a specification of the path to the genome reference files and the raw data fastq file. Please make sure to specify absolute paths rather than relative paths in your input JSON files.

1) [Input JSON file specification (short)](docs/input_short.md)
2) [Input JSON file specification (long)](docs/input.md)


## Running on local computer/HPCs

You can use URIs(`s3://`, `gs://` and `http(s)://`) in Caper's command lines and input JSON file then Caper will automatically download/localize such files. Input JSON file example: https://storage.googleapis.com/encode-pipeline-test-samples/encode-atac-seq-pipeline/ENCSR356KRQ_subsampled.json

According to your chosen platform of Caper, run Caper or submit Caper command line to the cluster. You can choose other environments like `--singularity` or `--docker` instead of `--conda`. But you must define one of the environments.

The followings are just examples. Please read [Caper's README](https://github.com/ENCODE-DCC/caper) very carefully to find an actual working command line for your chosen platform.
    ```bash
    # Run it locally with Conda (You don't need to activate it, make sure to install Conda envs first)
    $ caper run atac.wdl -i https://storage.googleapis.com/encode-pipeline-test-samples/encode-atac-seq-pipeline/ENCSR356KRQ_subsampled.json --conda

    # Or submit it as a leader job (with long/enough resources) to SLURM (Stanford Sherlock) with Singularity
    # It will fail if you directly run the leader job on login nodes
    $ sbatch -p [SLURM_PARTITON] -J [WORKFLOW_NAME] --export=ALL --mem 4G -t 4-0 --wrap "caper run atac.wdl -i https://storage.googleapis.com/encode-pipeline-test-samples/encode-atac-seq-pipeline/ENCSR356KRQ_subsampled.json --singularity"

    # Check status of your leader job
    $ squeue -u $USER | grep [WORKFLOW_NAME]

    # Cancel the leader node to close all of its children jobs
    $ scancel -j [JOB_ID]
    ```


## Running and sharing on Truwl

You can run this pipeline on [truwl.com](https://truwl.com/). This provides a web interface that allows you to define inputs and parameters, run the job on GCP, and monitor progress. To run it you will need to create an account on the platform then request early access by emailing [info@truwl.com](mailto:info@truwl.com) to get the right permissions. You can see the example case from this repo at [https://truwl.com/workflows/instance/WF_e85df4.f10.8880/command](https://truwl.com/workflows/instance/WF_e85df4.f10.8880/command). The example job (or other jobs) can be forked to pre-populate the inputs for your own job.

If you do not run the pipeline on Truwl, you can still share your use-case/job on the platform by getting in touch at [info@truwl.com](mailto:info@truwl.com) and providing your inputs.json file.


## Running on Terra/Anvil (using Dockstore)

Visit our pipeline repo on [Dockstore](https://dockstore.org/workflows/github.com/ENCODE-DCC/atac-seq-pipeline). Click on `Terra` or `Anvil`. Follow Terra's instruction to create a workspace on Terra and add Terra's billing bot to your Google Cloud account.

Download this [test input JSON for Terra](https://storage.googleapis.com/encode-pipeline-test-samples/encode-atac-seq-pipeline/ENCSR356KRQ_subsampled.terra.json) and upload it to Terra's UI and then run analysis.

If you want to use your own input JSON file, then make sure that all files in the input JSON are on a Google Cloud Storage bucket (`gs://`). URLs will not work.


## Running on DNAnexus (using Dockstore)

Sign up for a new account on [DNAnexus](https://platform.dnanexus.com/) and create a new project on either AWS or Azure. Visit our pipeline repo on [Dockstore](https://dockstore.org/workflows/github.com/ENCODE-DCC/atac-seq-pipeline). Click on `DNAnexus`. Choose a destination directory on your DNAnexus project. Click on `Submit` and visit DNAnexus. This will submit a conversion job so that you can check status of it on `Monitor` on DNAnexus UI.

Once conversion is done download one of the following input JSON files according to your chosen platform (AWS or Azure) for your DNAnexus project:
- AWS: https://storage.googleapis.com/encode-pipeline-test-samples/encode-atac-seq-pipeline/ENCSR356KRQ_subsampled_dx.json
- Azure: https://storage.googleapis.com/encode-pipeline-test-samples/encode-atac-seq-pipeline/ENCSR356KRQ_subsampled_dx_azure.json

You cannot use these input JSON files directly. Go to the destination directory on DNAnexus and click on the converted workflow `atac`. You will see input file boxes in the left-hand side of the task graph. Expand it and define FASTQs (`fastq_repX_R1` and also `fastq_repX_R2` if it's paired-ended) and `genome_tsv` as in the downloaded input JSON file. Click on the `common` task box and define other non-file pipeline parameters. e.g. `auto_detect_adapters` and `paired_end`.

We have a separate project on DNANexus to provide example FASTQs and `genome_tsv` for `hg38` and `mm10`. We recommend to make copies of these directories on your own project.

`genome_tsv`
- AWS: https://platform.dnanexus.com/projects/BKpvFg00VBPV975PgJ6Q03v6/data/pipeline-genome-data/genome_tsv/v4
- Azure: https://platform.dnanexus.com/projects/F6K911Q9xyfgJ36JFzv03Z5J/data/pipeline-genome-data/genome_tsv/v4

Example FASTQs
- AWS: https://platform.dnanexus.com/projects/BKpvFg00VBPV975PgJ6Q03v6/data/pipeline-test-samples/encode-atac-seq-pipeline/ENCSR356KRQ/fastq_subsampled
- Azure: https://platform.dnanexus.com/projects/F6K911Q9xyfgJ36JFzv03Z5J/data/pipeline-test-samples/encode-atac-seq-pipeline/ENCSR356KRQ/fastq_subsampled


## Running on DNAnexus (using our pre-built workflows)

See [this](docs/tutorial_dx_web.md) for details.


## How to organize outputs

Install [Croo](https://github.com/ENCODE-DCC/croo#installation). **You can skip this installation if you have installed pipeline's Conda environment and activated it**. Make sure that you have python3(> 3.4.1) installed on your system. Find a `metadata.json` on Caper's output directory.

```bash
$ pip install croo
$ croo [METADATA_JSON_FILE]
```

## How to make a spreadsheet of QC metrics

Install [qc2tsv](https://github.com/ENCODE-DCC/qc2tsv#installation). Make sure that you have python3(> 3.4.1) installed on your system. 

Once you have [organized output with Croo](#how-to-organize-outputs), you will be able to find pipeline's final output file `qc/qc.json` which has all QC metrics in it. Simply feed `qc2tsv` with multiple `qc.json` files. It can take various URIs like local path, `gs://` and `s3://`.

```bash
$ pip install qc2tsv
$ qc2tsv /sample1/qc.json gs://sample2/qc.json s3://sample3/qc.json ... > spreadsheet.tsv
```

QC metrics for each experiment (`qc.json`) will be split into multiple rows (1 for overall experiment + 1 for each bio replicate) in a spreadsheet.
