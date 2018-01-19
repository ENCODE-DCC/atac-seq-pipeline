ENCODE ATAC-seq pipeline test
===================================================

# How to generate reference outputs

1) Generate base reference outputs by running the following `.sh`'s.
```
$ bash ENCSR356KRQ_chr19.sh # PE test sample
$ bash ENCSR889WQX_chr19.sh # SE test sample
```

2) Wait until 1) is done, run the others.

3) Link generated input/reference outputs (starting with `gs://`) on GC (Google Cloud) into `test_atac.json`. `se_` prefix is for SE samples and `pe_` for PE ones.

# How to run cromwell server on GC

1) Create/restart an instance.

2) SSH to the server and run:
```
$ cd atac-seq-pipeline/test
$ bash run_cromwell_server.sh
```

# How `test_atac.sh` works

`test_atac.wdl` is not stand-alone. Cromwell/WDL supports importing a sub-workflow by `import` but the sub-workflow WDL to be imported should be in the remote workig directory (or URL) where cromwell server runs on GC. There is no way for Jenkins to update `atac.wdl` in this working directory on GC. So we are not using this `import "../atac.wdl"` feature.

`test_atac.sh` appends all task blocks (`task {}`) in `../atac.wdl` to `test_atac.wdl` and sends POST message with the appended `test_atac_standalone.wdl` to the cromwell server running on GC.

Task blocks in `../atac.wdl` should be wrapped in decorator-like comments:
```
#@TASK_DEF_BEGIN

task trim_adapter {
	...
}

task bowtie2 {
	...
}
...

#@TASK_DEF_END
```
