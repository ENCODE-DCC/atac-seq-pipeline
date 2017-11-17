Issues and discussion
===================================================

# WDL workflow/task structure

```
workflow wf {
  File f1 # = "/path/to/file1"
  call t1 { input : f = f1 }
  call t2 { input : f = t1.out }
}

task t1 {
  File f
  command {
    ...
  }
  output {
    File out = glob("*.bam")[0]
  }
}

task t2 {
  File f
  command {
    ...
  }
  output {
    ...
  }
}

input.json 
{
  "wf.f1" : "/path/to/file1"
}
...
```

# Cromwell file localization/duplication strategy

```
hard link -> soft link -> copy
```

Every task has its own working directory. If task-B's input comes from task-A's output then Cromwell creates a `input/` directory for task-B's and make a hard or soft link of task-A's output on `input/`.

For hard-linking, original file and `input/` directory must be on the same file system. For example, if original file is on NFS mounted directory and `input/` is on local machine then cromwell always makes a soft link.

You can change this strategy in backend `.conf` file. See [reference.conf](https://github.com/broadinstitute/cromwell/blob/29_hotfix/core/src/main/resources/reference.conf) for details.


# Variables

Java-like variables are available in WDL. There are type coercions among these variables. For exmaple, `String` to `File`, `Int` to `String`, ...
```
File, String, Int, Float, ...
```

Entries in `Array[] a` and `Map[] a` are accessible by `a[index_or_key]`. Entries in a `Pair[A,B] a` are acccesible by `a.left` and `a.right`.
```
Array[File], Array[String], Array[Array[File]], ...
Map[String, File], Map[String, Int], ...
Pair[String, Int], ...
```

## mandatory input variable
All mandatory variables must be specified in workflow/task level or in input.json

```
workflow wf {
  String f1
  File f2 = "/path/to/file2"

  call t1 { input : f1 = f1 }
}

task t1 {
  File f1
  File f2
  command {
    ...
  }
  output {
    ...
  }
}

input.json
{
  "wf.f" : "/path/to/file1"
  "wf.t1.f2" : "/path/to/file2"
}
```

## optional input variable and `if()`
There are two important functions to deal with optional variables. Optional variables default to `null` if not specified.
- `select_first([A,B,C,...])` : returns the first valid (non-null) variable in the array.
- `defined(A)` : return A!=null

- `if() then () else ()` : only for RHS
- `if() {}` : allowed as a code block only. No `else {}` is allowed

```
workflow wf {
  File f1
  File? f2
  Array[File]? arr  # error there is no type coercion from opt. array to array
  Boolean? flag1
  Boolean? flag2

  if ( flag1 ) {
    # error when flag1 is not specified in input.json
    call t1 { input: f1 = f1 }
  }
  else {
    # error. there is no if else both in workflow and task namespace
    # else is only allowed in RHS
  }
  if ( !flag1 ) { # this works
  }

  if ( select_first([flag1,false]) ) {
    call t1 as t1_special { input : f1 = if defined(flag2) then f2 else f1 } 
  }
}

task t1 {
  File f1
  command {
    ...
  }
  output {
    File out = ...
  }
}

input.json {
  "wf.flag1" : true
}
```

# Useful functions

`read_json()`, `write_json()` are not supported yet. If there is any undefined optional variable in `${}` then whole expression `${}` is ignored.

* `stdout()` : returns stdout of a task
* `glob()` : returns list of files
* `read_*()` : `read_lines()`, `read_tsv()`
* `write_*()` : `write_lines()`, `write_tsv()`
* `sep=delimiter arr` : returns `delimiter.join(arr)`

```
task t1 {
  Array[Array[String]] arr1
  Array[String] arr2
  String? opt

  command {
    python process_2_dim_array.py \
      ${"--opt " + opt} \
      --tsv ${write_tsv(arr1)} \
      --out-tsv arr1.tsv
    python process_1_dim_array.py \
      ${sep=' ' arr2} \
      --out-txt arr2.txt
    echo "hello world"
  }
  output {
    Array[Array[String]] out1 = read_tsv("arr1.tsv")
    Array[String] out2 = read_lines("arr2.txt")
    String out = stdout()
  }
}
```

# Parallelization `scatter()`

`scatter()` works very similar to MPI `scatter()`, but there is no `gather()` function because output of a task is automtically gathered in the namespace of a task.

Parallelization is **ALWAYS MAXIMIZED** in WDL. Cromwell first makes a task input/output dependency tree and start **ALL** tasks with **given inputs**. If task-B's input comes from another task-A's output, then task-B waits for task-A.


# ATAC-Seq pipeline.WDL

Goals
1) can start from any type of input files
2) can customize/fine-tune pipeline (need to have most flags and parameters in original BDS pipeline)
3) can run on any platforms (local, SGE, SLURM, Google Cloud, AWS, ...) with or w/o docker
4) separate genome specific files and parameters from workflow input JSON

# Issues

## Resumability (call-caching)

* One MYSQL database on a leader node manages all.
* We lose resumability
  - if docker tag (hash) changes.
  - if there is any change in `command {}` of a task.
  - if contents of workflow output directory move.

So workaround for this is to make the workflow start with any type of inputs (fastq, bam, tag-align, peak, ...) so that it can resume from any stage of the pipeline.

## Output directory structure on cloud platforms

* bucket for outputs
```
  [accession_id_or_sample_name]
    cromwell-executions
      atac
        [hash_str1]
          report.html # standalone HTML with embedded text and images (svg, png)
          output.json
          input_to_resume_from_[STAGE].json # to resume workflow from each stage of the workflow
        [hash_str2] # resumed from a certain [STAGE], including file and QC info of previous pipeline runs
          report.html
          output.json
          ...
        [hash_str3]
          report.html
          output.json
          ...
        ... 
```

* bucket for metadata (`output.json` only)
```
  [accession_id_or_sample_name]
    output.json -> [accession_id_or_sample_name]/.../[hash_str_latest]/output.json
```

## dxWDL

- Three dimensional array for `fastqs` and `adapters` are supported.
- `if {}` in `if {}` not supported
- `if {}` in `scatter() {}` not supported
- `scatter() {}` in `if {}` not supported

## Accession

```
atac_accession.wdl (atac_accession.py)
  check md5sum of all files to be accessioned in `output.json`.
    if md5sum matches that in the portal then ignore it
    otherwise, acession it
      remove all QC objects attached to it
      naming for new alias (md5sum?) for new file
      upload QC objects
      ...
```
