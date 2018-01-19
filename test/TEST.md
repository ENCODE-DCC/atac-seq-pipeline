ENCODE ATAC-seq pipeline test
===================================================

# How to generate reference outputs

1) Specify correct file paths in `subsample_fastq.sh` and run to subsample test samples (to 1/200 reads).
```
$ cd test_sample
$ bash subsample_fastq.sh
```

2) Generate base reference outputs by running the following shell scripts. These test samples are subsampled down to 1/200 reads.
```
$ cd test_sample
$ bash ENCSR356KRQ.sh
$ bash ENCSR889WQX.sh
```

3) Wait until 2) is done. Link outputs of 2) to JSON files in `test_sample/*.sh`, run other shell scripts.
```
$ cd test_sample
$ bash ENCSR356KRQ_disable_tn5_shift.sh
$ bash ENCSR356KRQ_no_dup_removal.sh
$ bash ENCSR356KRQ_no_multimapping.sh
$ bash ENCSR356KRQ_subsample.sh
$ bash ENCSR356KRQ_subsample_xcor.sh
$ bash ENCSR889WQX_disable_tn5_shift.sh
$ bash ENCSR889WQX_no_dup_removal.sh
$ bash ENCSR889WQX_no_multimapping.sh
$ bash ENCSR889WQX_subsample.sh
$ bash ENCSR889WQX_subsample_xcor.sh
```

4) Link generated input/reference outputs (starting with `gs://`) on GC (Google Cloud) into `test_atac.json`. `se_` prefix is for SE samples and `pe_` for PE ones.

# How to run cromwell server on GC

1) Create/restart an instance with name `encode-cromwell-test-server`. Choose 1vCPU and 4GB memory. Choose zone `us-west1-a`. Choose image `Ubuntu 16.04 (xenial)` with `Standard persistent disk 20GB`. Check the followings in Firewall section.
```
Allow HTTP traffic
Allow HTTPS traffic
```

2) SSH to the instance and run the followings to install cromwell, Docker and Java 8:
```
$ sudo su
$ cd
$ wget https://github.com/broadinstitute/cromwell/releases/download/30.1/cromwell-30.1.jar
$ chmod +x cromwell*.jar
$ echo 'export PATH=$PATH:/root'>> ~/.bashrc
$ source ~/.bashrc
$ apt-get update
$ apt-get install docker.io default-jre
```

3) Clone pipeline and run `MySQL` container.
```
$ sudo su
$ cd
$ mkdir /cromwell_db
$ git clone https://github.com/ENCODE-DCC/atac-seq-pipeline
$ docker run -d --name mysql-cromwell -v /cromwell_db:/var/lib/mysql -v /root/atac-seq-pipeline/docker_image/mysql:/docker-entrypoint-initdb.d -e MYSQL_ROOT_PASSWORD=cromwell -e MYSQL_DATABASE=cromwell_db --publish 3306:3306 mysql
$ docker ps
```

4) Run Cromwell server
```
$ sudo su
$ cd /root/atac-seq-pipeline
$ git checkout develop_test_jenkins
$ cd test
$ bash run_cromwell_server.sh
```

# How `test_atac.sh` works

`test_atac.wdl` is not stand-alone. Cromwell/WDL supports importing a sub-workflow by `import` but the sub-workflow WDL to be imported should be in the remote workig directory (or URL) where cromwell server runs on GC. There is no way for Jenkins to update `atac.wdl` in this working directory on GC. So we are not using this `import "../atac.wdl"` feature.

`test_atac.sh` appends all task blocks (`task {}`) in `../atac.wdl` to `test_atac.wdl` and sends POST message with the appended `test_atac_standalone.wdl` to the cromwell server running on GC.

Task blocks in `../atac.wdl` should be wrapped in decorator-like comments:
```
#@WORKFLOW_DEF_BEGIN

workflow atac {
...
}

#@WORKFLOW_DEF_END

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
