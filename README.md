# ENCODE ATAC-seq pipeline

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.156534.svg)](https://doi.org/10.5281/zenodo.156534)[![CircleCI](https://circleci.com/gh/ENCODE-DCC/atac-seq-pipeline/tree/master.svg?style=svg)](https://circleci.com/gh/ENCODE-DCC/atac-seq-pipeline/tree/master)

## Introduction

This pipeline is designed for automated end-to-end quality control and processing of ATAC-seq or DNase-seq data. The pipeline can be run on compute clusters with job submission engines or stand alone machines. It inherently makes uses of parallelized/distributed computing. Pipeline installation is also easy as most dependencies are automatically installed. The pipeline can be run end-to-end i.e. starting from raw FASTQ files all the way to peak calling and signal track generation; or can be started from intermediate stages as well (e.g. alignment files). The pipeline supports single-end or paired-end ATAC-seq or DNase-seq data (with or without replicates). The pipeline produces formatted HTML reports that include quality control measures specifically designed for ATAC-seq and DNase-seq data, analysis of reproducibility, stringent and relaxed thresholding of peaks, fold-enrichment and pvalue signal tracks.  The pipeline also supports detailed error reporting and easy resuming of runs. The pipeline has been tested on human, mouse and yeast ATAC-seq data and human and mouse DNase-seq data.

The ATAC-seq pipeline specification is also the official pipeline specification of the Encyclopedia of DNA Elements (ENCODE) consortium. The ATAC-seq pipeline protocol definition is [here](https://docs.google.com/document/d/1f0Cm4vRyDQDu0bMehHD7P7KOMxTOP-HiNoIvL1VcBt8/edit?usp=sharing). Some parts of the ATAC-seq pipeline were developed in collaboration with Jason Buenrostro, Alicia Schep and Will Greenleaf at Stanford.

### Features

* **Flexibility**: Support for `docker`, `singularity` and `Conda`.
* **Portability**: Support for many cloud platforms (Google/DNAnexus) and cluster engines (SLURM/SGE/PBS).
* **User-friendly HTML report**: tabulated quality metrics including alignment/peak statistics and FRiP along with many useful plots (IDR/cross-correlation measures).
  - Examples: [HTML](https://storage.googleapis.com/encode-pipeline-test-samples/encode-atac-seq-pipeline/ENCSR889WQX/example_output/qc.html), [JSON](docs/example_output/v1.1.5/qc.json)
* **ATAqC**: Annotation-based analysis including TSS enrichment and comparison to Roadmap DNase.
* **Genomes**: Pre-built database for GRCh38, hg19, mm10, mm9 and additional support for custom genomes.

## Installation and tutorial

This pipeline supports many cloud platforms and cluster engines. It also supports `docker`, `singularity` and `Conda` to resolve complicated software dependencies for the pipeline. A tutorial-based instruction for each platform will be helpful to understand how to run pipelines. There are special instructions for two major Stanford HPC servers (SCG4 and Sherlock).

* Cloud platforms
  * Web interface
    * [DNAnexus Platform](docs/tutorial_dx_web.md)
  * CLI (command line interface)
    * [Google Cloud Platform](docs/tutorial_google.md)
    * [DNAnexus Platform](docs/tutorial_dx_cli.md)
* Stanford HPC servers (CLI)
  * [Stanford SCG4](docs/tutorial_scg.md)
  * [Stanford Sherlock 2.0](docs/tutorial_sherlock.md)
* Cluster engines (CLI)
  * [SLURM](docs/tutorial_slurm.md)
  * [Sun GridEngine (SGE/PBS)](docs/tutorial_sge.md)
* Local Linux computers (CLI)
  * [Local system with `singularity`](docs/tutorial_local_singularity.md)
  * [Local system with `docker`](docs/tutorial_local_docker.md)
  * [Local system with `Conda`](docs/tutorial_local_conda.md)
* Local Windows computers (CLI)
  * [Windows 10 Pro with `docker`](docs/tutorial_windows_docker.md)
  * [Windows 10 Pro/Home with `Conda`](docs/tutorial_windows_conda.md)

## Input JSON file

[Input JSON file specification](docs/input.md)

## Output directories

[Output directory specification](docs/output.md)

## Useful tools

There are some useful tools to post-process outputs of the pipeline.

### qc_jsons_to_tsv

[This tool](utils/qc_jsons_to_tsv/README.md) recursively finds and parses all `qc.json` (pipeline's [final output](docs/example_output/v1.1.5/qc.json)) found from a specified root directory. It generates a TSV file that has all quality metrics tabulated in rows for each experiment and replicate. This tool also estimates overall quality of a sample by [a criteria definition JSON file](utils/qc_jsons_to_tsv/criteria.default.json) which can be a good guideline for QC'ing experiments.
