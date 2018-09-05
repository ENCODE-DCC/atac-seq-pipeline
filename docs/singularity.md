Binding directories to singularity
==================================

CHECK YOUR SINGULARITY VERSION FIRST AND UPGRADE YOUR SINGULARITY TO A VERSION `>=2.5.2` OR PIPELINE WILL NOT WORK CORRECTLY.
```
$ singularty --version
```

While [docker](https://www.docker.com/) is supported as a native backend in [cromwell](https://github.com/broadinstitute/cromwell), singularity runs as a custom backend. Therefore, users need to manually bind their input/output directories to singuarlity.

Singularity defaults to bind your `$HOME` directory and working directory (where you run a pipeline command line) automatically. This means that if you have input data/genome data located on your `$HOME` or `$PWD` then you don't need to define `singularity_command_options` in a workflow options JSON file below.

* Input data directory (your FASTQs, ...).
* Genome database directory (built from [this instruction](build_genome_database.md)).
* Output directory (where you run your pipelines).

A binding of directories is configured in a workflow options JSON file in `/workflow_opts`. See an example of workflow options file `workflow_opts/scg.json` for Stanford SCG cluster.
```
  {
      "default_runtime_attributes" : {
          "singularity_container" : "~/.singularity/atac-seq-pipeline-v1.1.simg",
          "singularity_command_options" : "--bind /ifs/scratch/leepc12/pipeline_genome_data,/ifs/scratch/leepc12/pipeline_test_samples"
      }
  }
```

Look for `--bind` in `singularity_command_options`. You will find that two ROOT directories are already defined there for binding pre-built genome database (`/ifs/scratch/leepc12/pipeline_genome_data`) and test input data (`/ifs/scratch/leepc12/pipeline_test_samples`).

If you don't even find `singularity_command_options` in your workflow options JSON then add one. If you don't find these two data directories then add your own to `singularity_command_options`. It's comma-separated.