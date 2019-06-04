# ENCODE ATAC-seq pipeline

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.156534.svg)](https://doi.org/10.5281/zenodo.156534)[![CircleCI](https://circleci.com/gh/ENCODE-DCC/atac-seq-pipeline/tree/master.svg?style=svg)](https://circleci.com/gh/ENCODE-DCC/atac-seq-pipeline/tree/master)

## Introduction

This pipeline is designed for automated end-to-end quality control and processing of ATAC-seq or DNase-seq data. The pipeline can be run on compute clusters with job submission engines or stand alone machines. It inherently makes uses of parallelized/distributed computing. Pipeline installation is also easy as most dependencies are automatically installed. The pipeline can be run end-to-end i.e. starting from raw FASTQ files all the way to peak calling and signal track generation; or can be started from intermediate stages as well (e.g. alignment files). The pipeline supports single-end or paired-end ATAC-seq or DNase-seq data (with or without replicates). The pipeline produces formatted HTML reports that include quality control measures specifically designed for ATAC-seq and DNase-seq data, analysis of reproducibility, stringent and relaxed thresholding of peaks, fold-enrichment and pvalue signal tracks.  The pipeline also supports detailed error reporting and easy [resuming](utils/qc_jsons_to_tsv/README.md) of runs. The pipeline has been tested on human, mouse and yeast ATAC-seq data and human and mouse DNase-seq data.

The ATAC-seq pipeline specification is also the official pipeline specification of the Encyclopedia of DNA Elements (ENCODE) consortium. The ATAC-seq pipeline protocol definition is [here](https://docs.google.com/document/d/1f0Cm4vRyDQDu0bMehHD7P7KOMxTOP-HiNoIvL1VcBt8/edit?usp=sharing). Some parts of the ATAC-seq pipeline were developed in collaboration with Jason Buenrostro, Alicia Schep and Will Greenleaf at Stanford.

### Features

* **Portability**: Support for many cloud platforms (Google/DNAnexus) and cluster engines (SLURM/SGE/PBS).
* **User-friendly HTML report**: tabulated quality metrics including alignment/peak statistics and FRiP along with many useful plots (IDR/cross-correlation measures).
  - Examples: [HTML](https://storage.googleapis.com/encode-pipeline-test-samples/encode-atac-seq-pipeline/ENCSR889WQX/example_output/qc.html), [JSON](docs/example_output/v1.1.5/qc.json)
* **ATAqC**: Annotation-based analysis including TSS enrichment and comparison to Roadmap DNase.
* **Genomes**: Pre-built database for GRCh38, hg19, mm10, mm9 and additional support for custom genomes.

## Installation

1) Install [Caper](https://github.com/ENCODE-DCC/caper#installation). Caper is a python wrapper for [Cromwell](https://github.com/broadinstitute/cromwell). Make sure that you have python3(> 3.4.1) installed on your system.

  ```bash
  $ pip install caper
  ```

2) Read through [Caper's README](https://github.com/ENCODE-DCC/caper) carefully.

3) Run a pipeline with Caper.

## Running pipelines without Caper

Caper uses the cromwell workflow execution engine to run the workflow on the platform you specify.  While we recommend you use caper, if you want to run cromwell directly without caper you can learn about that [here](docs/deprecated/OLD_METHOD.md).

## DNAnexus

You can also run our pipeline on DNAnexus without using Caper or Cromwell. There are two ways to build a workflow on DNAnexus based on our WDL.

1) [dxWDL CLI](docs/tutorial_dx_cli.md)
2) [DNAnexus Web UI](docs/tutorial_dx_web.md)

## Conda

We no longer recommend Conda for resolving dependencies and plan to phase out Conda support. Instead we recommend using Docker or Singularity. You can install Singularity and use it for our pipeline with Caper (by adding `--use-singularity` to command line arguments).

1) Install [Conda](https://docs.conda.io/en/latest/miniconda.html).

2) Install Conda environment for pipeline.

  ```bash
  $ conda/install_dependencies.sh
  ```

3) Initialize Conda and re-login.

  ```bash
  $ conda init bash
  $ exit
  ```

4) Configure pipeline's python2 and python3 environments.

  ```bash
  $ conda/config_conda_env.sh
  $ conda/config_conda_env_py3.sh
  ```

5) Update pipeline's Conda environment with pipeline's python source code. You need to run this step everytime you update (`git pull`) this pipeline.

  ```bash
  $ conda/update_conda_env.sh
  ```

## Tutorial

Make sure that you have configured Caper correctly.
> **WARNING**: DO NOT RUN THIS ON HPC LOGIN NODES. YOUR JOBS WILL BE KILLED.

Run it. Due to `--deepcopy` all files in `examples/caper/ENCSR356KRQ_subsampled.json` will be recursively copied into Caper's temporary folder (`--tmp-dir`).
```bash
$ caper run atac.wdl -i examples/caper/ENCSR356KRQ_subsampled.json --deepcopy --use-singularity
```

If you use Conda or Docker (on cloud platforms) then remove `--use-singularity` from the command line and activate it before running a pipeline.
```bash
$ conda activate encode-atac-seq-pipeline
$ caper run atac.wdl -i examples/caper/ENCSR356KRQ_subsampled.json --deepcopy
```

To run it on an HPC (e.g. Stanford Sherlock and SCG). See details at [Caper's README](https://github.com/ENCODE-DCC/caper/blob/master/README.md#how-to-run-it-on-slurm-cluster).

## Input JSON file

Always use absolute paths in an input JSON.

[Input JSON file specification](docs/input.md)

## How to organize outputs

Install [Croo](https://github.com/ENCODE-DCC/croo#installation). Make sure that you have python3(> 3.4.1) installed on your system.

```bash
$ pip install croo
```

Find a `metadata.json` on Caper's output directory.

```bash
$ croo [METADATA_JSON_FILE]
```

## How to build/download genome database

You need to specify a genome data TSV file in your input JSON. Such TSV can be generated/downloaded with actual genome database files.

Use genome database [downloader](genome/download_genome_data.sh) or [builder](docs/build_genome_database.md) for your own genome.

## Useful tools

There are some useful tools to post-process outputs of the pipeline.

### qc_jsons_to_tsv

[This tool](utils/qc_jsons_to_tsv/README.md) recursively finds and parses all `qc.json` (pipeline's [final output](docs/example_output/v1.1.5/qc.json)) found from a specified root directory. It generates a TSV file that has all quality metrics tabulated in rows for each experiment and replicate. This tool also estimates overall quality of a sample by [a criteria definition JSON file](utils/qc_jsons_to_tsv/criteria.default.json) which can be a good guideline for QC'ing experiments.

### ENCODE downloader

[This tool](https://github.com/kundajelab/ENCODE_downloader) downloads any type (FASTQ, BAM, PEAK, ...) of data from the ENCODE portal. It also generates a metadata JSON file per experiment which will be very useful to make an input JSON file for the pipeline.