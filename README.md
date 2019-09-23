# ENCODE ATAC-seq pipeline

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.156534.svg)](https://doi.org/10.5281/zenodo.156534)[![CircleCI](https://circleci.com/gh/ENCODE-DCC/atac-seq-pipeline/tree/master.svg?style=svg)](https://circleci.com/gh/ENCODE-DCC/atac-seq-pipeline/tree/master)

## Introduction

This pipeline is designed for automated end-to-end quality control and processing of ATAC-seq or DNase-seq data. The pipeline can be run on compute clusters with job submission engines or stand alone machines. It inherently makes uses of parallelized/distributed computing. Pipeline installation is also easy as most dependencies are automatically installed. The pipeline can be run end-to-end i.e. starting from raw FASTQ files all the way to peak calling and signal track generation; or can be started from intermediate stages as well (e.g. alignment files). The pipeline supports single-end or paired-end ATAC-seq or DNase-seq data (with or without replicates). The pipeline produces formatted HTML reports that include quality control measures specifically designed for ATAC-seq and DNase-seq data, analysis of reproducibility, stringent and relaxed thresholding of peaks, fold-enrichment and pvalue signal tracks.  The pipeline also supports detailed error reporting and easy resuming of runs. The pipeline has been tested on human, mouse and yeast ATAC-seq data and human and mouse DNase-seq data.

The ATAC-seq pipeline protocol definition is [here](https://docs.google.com/document/d/1f0Cm4vRyDQDu0bMehHD7P7KOMxTOP-HiNoIvL1VcBt8/edit?usp=sharing). Some parts of the ATAC-seq pipeline were developed in collaboration with Jason Buenrostro, Alicia Schep and Will Greenleaf at Stanford.

### Features

* **Portability**: The pipeline run can be performed across different cloud platforms such as Google, AWS and DNAnexus, as well as on cluster engines such as SLURM, SGE and PBS.
* **User-friendly HTML report**: In addition to the standard outputs, the pipeline generates an HTML report that consists of a tabular representation of quality metrics including alignment/peak statistics and FRiP along with many useful plots (IDR/TSS enrichment). An example of the [HTML report](https://storage.googleapis.com/encode-pipeline-test-samples/encode-atac-seq-pipeline/ENCSR889WQX/example_output/qc.html). The [json file](docs/example_output/v1.1.5/qc.json) used in generating this report.
* **Supported genomes**: Pipeline needs genome specific data such as aligner indices, chromosome sizes file and blacklist. We provide a genome database downloader/builder for hg38, hg19, mm10, mm9. You can also use this [builder](docs/build_genome_database.md) to build genome database from FASTA for your custom genome.

## Installation

1) Git clone this repo.

	```bash
	$ cd
	$ git clone https://github.com/ENCODE-DCC/atac-seq-pipeline
	```

2) Install [Caper](https://github.com/ENCODE-DCC/caper#installation). Caper is a python wrapper for [Cromwell](https://github.com/broadinstitute/cromwell).

	> **IMPORTANT**: Make sure that you have python3(> 3.4.1) installed on your system.

	```bash
	$ pip install caper  # use pip3 if it doesn't work
	```

3) Read through [Caper's README](https://github.com/ENCODE-DCC/caper) carefully. Find an instruction for your platform. 
	> **IMPORTANT**: Configure your Caper configuration file `~/.caper/default.conf` correctly for your platform.

## Running a pipeline locally with Caper

1) Prepare an input JSON file. We will use a subsampled example input JSON based on URLs. Caper will automatically download all fastqs and reference human genome data recursively.
	```bash
	$ INPUT_JSON=https://storage.googleapis.com/encode-pipeline-test-samples/encode-atac-seq-pipeline/ENCSR356KRQ_subsampled_caper.json
	```

2-1) **Conda**: Run a workflow with Conda. Make sure that you have followed [this instruction](docs/install_conda.md) to install Conda and its environments.
	> **WARNING**: We no longer recommend Conda for resolving dependencies and plan to phase out Conda support. Instead we recommend using Docker or Singularity. You can install Singularity and use it for our pipeline with Caper (by adding `--singularity` to command line arguments).

	```bash
	$ conda activate encode-atac-seq-pipeline
	$ caper run atac.wdl -i ${INPUT_JSON}
	```

2-2) **Singularity (RECOMMENDED)**: Run a workflow with Singularity.

	```bash
	$ caper run atac.wdl -i ${INPUT_JSON} --singularity
	```

	> **HPCs**: To run multiple workflows on HPCs (e.g. Stanford Sherlock and SCG) see details at [Caper's README](https://github.com/ENCODE-DCC/caper/blob/master/README.md#how-to-run-it-on-slurm-cluster). Do not run Caper on login nodes. Your workflows will get killed. There is a learning curve to understand server/client structure of Caper/Cromwell.


2-3) **Docker**: Run a workflow with Docker.

  ```bash
	$ caper run atac.wdl -i ${INPUT_JSON} --docker
  ```

3) You can also run a workflow on cloud platforms such as AWS (`aws`) and Google Cloud Platform (`gcp`) if Caper's configuration file is correctly configured for them. See details at [Caper's README](https://github.com/ENCODE-DCC/caper).


## Running a pipeline without Caper

Caper uses the cromwell workflow execution engine to run the workflow on the platform you specify. While we recommend you use caper, if you want to run cromwell directly without caper you can learn about that [here](docs/deprecated/OLD_METHOD.md).

## Running a pipeline on DNAnexus

You can also run our pipeline on DNAnexus without using Caper or Cromwell. There are two ways to build a workflow on DNAnexus based on our WDL.

1) [dxWDL CLI](docs/tutorial_dx_cli.md)
2) [DNAnexus Web UI](docs/tutorial_dx_web.md)

## Input JSON file

An input JSON file includes all genomic data files, input parameters and metadata for running pipelines. Always use absolute paths in an input JSON.

[Input JSON file specification](docs/input.md)

## How to organize outputs

Install [Croo](https://github.com/ENCODE-DCC/croo#installation). Make sure that you have python3(> 3.4.1) installed on your system. Find a `metadata.json` on Caper's output directory.

```bash
$ pip install croo
$ croo [METADATA_JSON_FILE]
```

## Useful tools

There are some useful tools to post-process outputs of the pipeline.

### qc_jsons_to_tsv

[This tool](utils/qc_jsons_to_tsv/README.md) recursively finds and parses all `qc.json` (pipeline's [final output](docs/example_output/v1.1.5/qc.json)) found from a specified root directory. It generates a TSV file that has all quality metrics tabulated in rows for each experiment and replicate. This tool also estimates overall quality of a sample by [a criteria definition JSON file](utils/qc_jsons_to_tsv/criteria.default.json) which can be a good guideline for QC'ing experiments.
