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
$ java -jar [JAVA_OPTS] ../cromwell-29.jar run atac.wdl -i input.json -o [WORKFLOW_OPTS]
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
     $ java -jar -Dconfig.file=non_docker/slurm.conf cromwell-29.jar run atac.wdl -i input.json
     ```

* For SGE without docker: Modify `non_docker/sge.conf` if your parallel environment (PE) is not `shm`.

     ```
     $ java -jar -Dconfig.file=non_docker/sge.conf cromwell-29.jar run atac.wdl -i input.json
     ```

# Dependency installation for platforms without docker support

**We do not recommend running our pipeline without docker!**. Use it with caution.

1) Install Miniconda on your system: [Instruction](https://github.com/kundajelab/atac_dnase_pipelines#conda)
2) Download [ATAC-Seq/DNase-Seq pipeline](https://github.com/kundajelab/atac_dnase_pipelines):
	 ```
	 $ git clone --recursive https://github.com/kundajelab/atac_dnase_pipelines
	 ```
3) Install dependencies.
	 ```
	 $ bash install_dependencies.sh
	 ```
4) Activate Conda environment and run a pipeline
	 ```
	 $ source activate bds_atac
	 $ java -jar [JAVA_OPTS] ../cromwell-29.jar run atac.wdl -i input.json
	 ```