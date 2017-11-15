Issues and discussion
===================================================

# 

# Issues

## Resumability (call-caching)

- One MYSQL database on a leader node manages all.
- We lose resumability
  - If docker tag (hash) changes.
  - If outputs moved.

## Output directory structure on cloud platforms

- bucket for outputs
```
  [accession_id_or_sample_name]
    cromwell-executions
      atac
        [hash_str1]
          report.html # standalone HTML which embeds all images (svg, png)
          output.json
        [hash_str2]
          report.html
          output.json
        [hash_str3]
          report.html
          output.json
      ... 
```
- bucket for metadata (output.json only)
```
  [accession_id_or_sample_name]
    output.json -> [accession_id_or_sample_name]/.../[hash_str_latest]/out.json
```

atac_accession.wdl (atac_accession.py)
  check md5sum of all files to be accessioned in output.json.
    if md5sum matches that in the portal then ignore it
    otherwise, acession it
      remove all QC objects attached to it
      naming for new alias (md5sum?) for new file
      upload QC objects


## Accession








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
