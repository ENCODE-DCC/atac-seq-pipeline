ENCODE ATAC-seq pipeline
===================================================

# Usage

General usage.

```
$ java -jar [JAVA_OPTS] ../cromwell-29.jar run atac.wdl -i input.json -o [CROMWELL_WF_OPTS]
```

* For platforms supporting docker 

     ```
     $ java -jar cromwell-29.jar run atac.wdl -i input.json -o docker.json
     ```

* For SLURM without docker
 
     ```
     $ java -jar -Dconfig.file=non_docker/slurm.conf cromwell-29.jar run atac.wdl -i input.json
     ```

* For SGE without docker

     ```
     $ java -jar -Dconfig.file=non_docker/sge.conf cromwell-29.jar run atac.wdl -i input.json
     ```
