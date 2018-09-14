ENCODE ATAC-seq pipeline
========================

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.156534.svg)](https://doi.org/10.5281/zenodo.156534)

This pipeline is designed for automated end-to-end quality control and processing of ATAC-seq or DNase-seq data. The pipeline can be run on compute clusters with job submission engines or stand alone machines. It inherently makes uses of parallelized/distributed computing. Pipeline installation is also easy as most dependencies are automatically installed. The pipeline can be run end-to-end i.e. starting from raw FASTQ files all the way to peak calling and signal track generation; or can be started from intermediate stages as well (e.g. alignment files). The pipeline supports single-end or paired-end ATAC-seq or DNase-seq data (with or without replicates). The pipeline produces pretty HTML reports that include quality control measures specifically designed for ATAC-seq and DNase-seq data, analysis of reproducibility, stringent and relaxed thresholding of peaks, fold-enrichment and pvalue signal tracks.  The pipeline also supports detailed error reporting and easy resuming of runs. The pipeline has been tested on human, mouse and yeast ATAC-seq data and human and mouse DNase-seq data.

The ATAC-seq pipeline specification is also the official pipeline specification of the Encyclopedia of DNA Elements (ENCODE) consortium. The ATAC-seq pipeline protocol definition is [here](https://docs.google.com/document/d/1f0Cm4vRyDQDu0bMehHD7P7KOMxTOP-HiNoIvL1VcBt8/edit?usp=sharing). Some parts of the ATAC-seq pipeline were developed in collaboration with Jason Buenrostro, Alicia Schep and Will Greenleaf at Stanford.

# Installation / Tutorial

* [Google Cloud Platform](docs/tutorial_google.md)
* [DNANexus Platform with dxWDL CLI](docs/tutorial_dx_cli.md)
* [DNANexus Platform with Web UI](docs/tutorial_dx_web.md)
* [Stanford SCG](docs/tutorial_scg.md)
* [Stanford Sherlock 2.0](docs/tutorial_sherlock.md)
* [SLURM](docs/tutorial_slurm.md)
* [Sun GridEngine (SGE)](docs/tutorial_sge.md)
* [Local system with singularity](docs/tutorial_local_singularity.md)
* [Local system with docker](docs/tutorial_local_docker.md)
* [Local system with Conda](docs/tutorial_local_conda.md)

# [Input JSON](docs/input.md)

# [Output](docs/output.md)
