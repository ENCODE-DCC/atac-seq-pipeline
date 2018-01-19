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

1) Create/restart an instance with the following settings.
* name : `encode-cromwell-test-server`. 
* resource: 1vCPU and 4GB memory
* zone: `us-west1-a`.
* image: `Ubuntu 16.04 (xenial)`
* disk: `Standard persistent disk 20GB`
* Network tags: add a tag `cromwell-server`.
* Cloud API access scopes: `Allow full access to all Cloud APIs`.
* External IP: any static IP address.

2) SSH to the instance and run the followings to install Docker and Java 8:
```
$ sudo apt-get update
$ sudo apt-get install docker.io default-jre
$ sudo usermod -aG docker $USER
```

3) Log out and log back in.

4) Install cromwell.
```
$ cd
$ wget https://github.com/broadinstitute/cromwell/releases/download/30.1/cromwell-30.1.jar
$ chmod +x cromwell*.jar
$ echo "export PATH=\$PATH:\$HOME">> ~/.bashrc
$ source ~/.bashrc
```

5) Clone pipeline, make DB directory (where metadata of all pipelines are stored) and run `MySQL` container.
```
$ cd
$ git clone https://github.com/ENCODE-DCC/atac-seq-pipeline
$ mkdir cromwell_db
$ docker run -d --name mysql-cromwell -v $HOME/cromwell_db:/var/lib/mysql -v $HOME/atac-seq-pipeline/docker_image/mysql:/docker-entrypoint-initdb.d -e MYSQL_ROOT_PASSWORD=cromwell -e MYSQL_DATABASE=cromwell_db --publish 3306:3306 mysql
$ docker ps
```

4) Run Cromwell server
```
$ cd $HOME/atac-seq-pipeline
$ git checkout develop_test_jenkins
$ cd test
$ screen -RD cromwell # make screen for cromwell server
$ bash run_cromwell_server_on_gc.sh
```

5) Firewall settings to open port 8000
* Go to cloud.google.com
* Go to my Console
* Choose you Project.
* Choose Networking > VPC network
* Choose "Firewalls rules"
* Choose Create Firewall Rule `encode-cromwell-test-server-open-port-8000`.
* Targets: `Specified target rags`.
* Target tags: cromwell-server
* Source IP ranges: `0.0.0.0/0`.
* Protocols and Ports: `Specified protocols and ports` with `tcp:8000`.

# How `test_atac.sh` works (for Jenkins)

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
