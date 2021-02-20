# ENCODE ATAC-seq pipeline

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.156534.svg)](https://doi.org/10.5281/zenodo.156534)[![CircleCI](https://circleci.com/gh/ENCODE-DCC/atac-seq-pipeline/tree/master.svg?style=svg)](https://circleci.com/gh/ENCODE-DCC/atac-seq-pipeline/tree/master)

## Important notice for Conda users

If it takes too long to resolve Conda package conflicts while installing pipeline's Conda environment, then try with `mamba` instead. Add `mamba` to the install command line.
```bash
$ scripts/install_conda_env.sh mamba
```

For every new pipeline release, Conda users always need to update pipeline's Conda environment (`encode-atac-seq-pipeline`) even though they don't use new added features.
```bash
$ cd atac-seq-pipeline
$ scripts/update_conda_env.sh
```

For pipelines >= v1.4.0, Conda users also need to manually install Caper and Croo **INSIDE** the environment. These two tools have been removed from pipeline's Conda environment since v1.9.2.
```bash
$ source activate encode-atac-seq-pipeline
$ pip install caper croo
```

## Introduction

This pipeline is designed for automated end-to-end quality control and processing of ATAC-seq and DNase-seq data. The pipeline can be run on compute clusters with job submission engines as well as on stand alone machines. It inherently makes uses of parallelized/distributed computing. Pipeline installation is also easy as most dependencies are automatically installed. The pipeline can be run end-to-end, starting from raw FASTQ files all the way to peak calling and signal track generation using a single caper submit command. One can also start the pipeline from intermediate stages (for example, using alignment files as input). The pipeline supports both single-end and paired-end data as well as replicated or non-replicated datasets. The outputs produced by the pipeline include 1) formatted HTML reports that include quality control measures specifically designed for ATAC-seq and DNase-seq data, 2) analysis of reproducibility, 3) stringent and relaxed thresholding of peaks, 4) fold-enrichment and pvalue signal tracks. The pipeline also supports detailed error reporting and allows for easy resumption of interrupted runs. It has been tested on some human, mouse and yeast ATAC-seq datasets as well as on human and mouse DNase-seq datasets.

The ATAC-seq pipeline protocol specification is [here](https://docs.google.com/document/d/1f0Cm4vRyDQDu0bMehHD7P7KOMxTOP-HiNoIvL1VcBt8/edit?usp=sharing). Some parts of the ATAC-seq pipeline were developed in collaboration with Jason Buenrostro, Alicia Schep and Will Greenleaf at Stanford.

### Features

* **Portability**: The pipeline run can be performed across different cloud platforms such as Google, AWS and DNAnexus, as well as on cluster engines such as SLURM, SGE and PBS.
* **User-friendly HTML report**: In addition to the standard outputs, the pipeline generates an HTML report that consists of a tabular representation of quality metrics including alignment/peak statistics and FRiP along with many useful plots (IDR/TSS enrichment). An example of the [HTML report](https://storage.googleapis.com/encode-pipeline-test-samples/encode-atac-seq-pipeline/ENCSR889WQX/example_output/qc.html). The [json file](https://storage.googleapis.com/encode-pipeline-test-samples/encode-atac-seq-pipeline/ENCSR889WQX/example_output/qc.json) used in generating this report.
* **Supported genomes**: Pipeline needs genome specific data such as aligner indices, chromosome sizes file and blacklist. We provide a genome database downloader/builder for hg38, hg19, mm10, mm9. You can also use this [builder](docs/build_genome_database.md) to build genome database from FASTA for your custom genome.

## Installation

1) Git clone this pipeline.
	> **IMPORTANT**: use `~/atac-seq-pipeline/atac.wdl` as `[WDL]` in Caper's documentation.

	```bash
	$ cd
	$ git clone https://github.com/ENCODE-DCC/atac-seq-pipeline
	```

2) Install pipeline's [Conda environment](docs/install_conda.md) if you want to use Conda instead of Docker/Singularity. Conda is recommneded on local computer and HPCs (e.g. Stanford Sherlock/SCG).
	> **IMPORTANT**: use `encode-atac-seq-pipeline` as `[PIPELINE_CONDA_ENV]` in Caper's documentation.

3) **Skip this step if you have installed pipeline's Conda environment**. Caper is already included in the Conda environment. [Install Caper](https://github.com/ENCODE-DCC/caper#installation). Caper is a python wrapper for [Cromwell](https://github.com/broadinstitute/cromwell).
	> **IMPORTANT**: Make sure that you have python3(>= 3.6.0) installed on your system.

	```bash
	$ pip install caper  # use pip3 if it doesn't work
	```

4) Follow [Caper's README](https://github.com/ENCODE-DCC/caper) carefully. Find an instruction for your platform.
	> **IMPORTANT**: Configure your Caper configuration file `~/.caper/default.conf` correctly for your platform.

## Test input JSON file

Use `https://storage.googleapis.com/encode-pipeline-test-samples/encode-atac-seq-pipeline/ENCSR356KRQ_subsampled.json` as `[INPUT_JSON]` in Caper's documentation.

## Input JSON file

> **IMPORTANT**: DO NOT BLINDLY USE A TEMPLATE/EXAMPLE INPUT JSON. READ THROUGH THE FOLLOWING GUIDE TO MAKE A CORRECT INPUT JSON FILE.

An input JSON file specifies all the input parameters and files that are necessary for successfully running this pipeline. This includes a specification of the path to the genome reference files and the raw data fastq file. Please make sure to specify absolute paths rather than relative paths in your input JSON files.

1) [Input JSON file specification (short)](docs/input_short.md)
2) [Input JSON file specification (long)](docs/input.md)


## Running and sharing on Truwl
You can run this pipeline on [truwl.com](https://truwl.com/). This provides a web interface that allows you to define inputs and parameters, run the job on GCP, and monitor progress. To run it you will need to create an account on the platform then request early access by emailing [info@truwl.com](mailto:info@truwl.com) to get the right permissions. You can see the example case from this repo at [https://truwl.com/workflows/instance/WF_e85df4.f10.8880/command](https://truwl.com/workflows/instance/WF_e85df4.f10.8880/command). The example job (or other jobs) can be forked to pre-populate the inputs for your own job.

If you do not run the pipeline on Truwl, you can still share your use-case/job on the platform by getting in touch at [info@truwl.com](mailto:info@truwl.com) and providing your inputs.json file.


## Running a pipeline on DNAnexus

You can also run this pipeline on DNAnexus without using Caper or Cromwell. There are two ways to build a workflow on DNAnexus based on our WDL.

1) [dxWDL CLI](docs/tutorial_dx_cli.md)
2) [DNAnexus Web UI](docs/tutorial_dx_web.md)

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
