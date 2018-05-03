ENCODE ATAC-seq pipeline test
===================================================

# Task level test (local)

This test requires `atac-seq-pipeline-test-data` directory in `test_task/`. Git glone [a data repo](https://github.com/leepc12/atac-seq-pipeline-test-data) on `test_task/`. This repo has 1/400 subsampled test samples and chr19-chrM only bowtie2 indices and other genome data for hg38 and mm10. Make sure that you have `cromwell-31.jar` in your `$PATH` as an executable (`chmod +x`) and `Docker` installed on your system.
```
$ cd test_task/
$ git clone https://github.com/encode-dcc/atac-seq-pipeline-test-data
```

Each task in `../atac.wdl` has a corresponding pair of tester WDL/JSON (`[TASK_NAME].WDL` and [TASK_NAME].json`). You can also specify your own docker image to test each task.
```
$ cd test_task/
$ ./test.sh [WDL] [INPUT_JSON] [DOCKER_IMAGE](optional)
```

# Workflow level test (on GC)

Make sure that you have a Cromwell server running on GC. This shell script will submit `../atac.wdl` to the server and wait for a response (`result.json`). There are two input JSON files (original and subsampled) for each endedness (SE and PE). You can also check all outputs on GC bucket `gs://encode-pipeline-test-runs`.
```
$ cd test_workflow/
$ ./test_atac.sh [INPUT_JSON] [QC_JSON_TO_COMPARE] [DOCKER_IMAGE](optional)
```

Jenkins must do the following:
```
$ cd test_workflow/
# For master branch (full test sample, ~24hr)
$ ./test_atac.sh ENCSR356KRQ.json ref_output/ENCSR356KRQ_qc.json [NEW_DOCKER_IMAGE]
$ ./test_atac.sh ENCSR889WQX.json ref_output/ENCSR889WQX_qc.json [NEW_DOCKER_IMAGE]
# For develop branch (1/400 subsampled and chr19 only test sample ~30mins)
$ ./test_atac.sh ENCSR356KRQ_subsampled.json ref_output/ENCSR356KRQ_subsampled_chr19_only_qc.json [NEW_DOCKER_IMAGE]
$ ./test_atac.sh ENCSR889WQX_subsampled.json ref_output/ENCSR889WQX_subsampled_chr19_only_qc.json [NEW_DOCKER_IMAGE]
```

`test_atac.sh` will generate the following files to validate pipeline outputs. Jenkins must check if `PREFIX.qc_json_diff.txt` is empty or not.
* `PREFIX.result.json`: all outputs of `atac.wdl`.
* `PREFIX.result.qc.json`: qc summary JSON file `qc.json` of `atac.wdl`.
* `PREFIX.qc_json_diff.txt`: diff between `PREFIX.result.qc.json` and reference in `ref_output/`.

# How to run a Cromwell server on GC

1) Create/restart an instance with the following settings.
* name : `encode-cromwell-test-server`. 
* resource: 1vCPU and 4GB memory
* zone: `us-west1-a`.
* image: `Ubuntu 16.04 (xenial)`
* disk: `Standard persistent disk 20GB`
* Network tags: add a tag `cromwell-server`.
* Cloud API access scopes: `Allow full access to all Cloud APIs`.
* External IP (optional): any static IP address.

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
$ wget https://github.com/broadinstitute/cromwell/releases/download/31/cromwell-31.jar
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
* Go to Google Cloud Console
* Choose your Project.
* Choose Networking > VPC network
* Choose "Firewalls rules"
* Choose Create Firewall Rule `encode-cromwell-test-server-open-port-8000`.
* Targets: `Specified target rags`.
* Target tags: cromwell-server
* Source IP ranges: `0.0.0.0/0` (CIDR notation for allowed IP range)
* Protocols and Ports: `Specified protocols and ports` with `tcp:8000`.
