pipeline {
        agent none
        environment {
                QUAY_USER = credentials('quay-robot')
                QUAY_PASS = credentials('quay-robot-token')
        }
        stages {
        		stage('Tag Non-master') {
                        agent {label 'master-builder'}
                        when { not {branch 'master'}}
            			steps { 
                                // the tag gets built here, and can be referenced in the other stages 
                                script {
                                        TAG = sh([script: "echo quay.io/encode-dcc/atac-seq-pipeline:${env.BRANCH_NAME}_${env.BUILD_NUMBER}", returnStdout: true]).trim()
                                       }
        				echo "On non-master"
            			}
                        
        		}

                stage('Tag Master') {
                        agent {label 'master-builder'}
                        when { branch 'master'}
                        steps { 
                                // the tag gets built here, and can be referenced in the other stages 
                                script {
                                        TAG = sh([script: "echo quay.io/encode-dcc/atac-seq-pipeline:latest", returnStdout: true]).trim()
                                       }
                        echo "On non-master"
                        }
                }

                stage('Build-nonmaster') {
                        agent {label 'slave-w-docker-cromwell-60GB-ebs'} 
                        when { not { branch 'master' } }
                        steps { 
                                echo "the tag is $TAG"
                                echo "going to build a docker image now.."
                                slackSend (color: '#7CFC00', message: "started job: ${env.JOB_NAME}, build number ${env.BUILD_NUMBER} on branch: ${env.BRANCH_NAME}.")
				                slackSend "The images will be tagged as $TAG"

                                // pull the cache template image (the image is going to stay pretty much the same so it is no need to be dynamic)
                                sh "docker pull quay.io/encode-dcc/atac-seq-pipeline:v1.1"
                                sh "docker login -u=${QUAY_USER} -p=${QUAY_PASS} quay.io"
                                sh "docker build --cache-from quay.io/encode-dcc/atac-seq-pipeline:v1.1 -f docker_image/Dockerfile -t atac-seq-pipeline ."
                                sh "docker tag atac-seq-pipeline $TAG"
                                sh "docker push $TAG"
                                sh "docker logout"
                        }
                }

                stage('Build-master') {
                        agent {label 'slave-w-docker-cromwell-60GB-ebs'} 
                        when { branch 'master'}
                        steps {
                                echo "going to build a docker image now.."
                                slackSend (color: '#7CFC00', message: "started job: ${env.JOB_NAME}, build number ${env.BUILD_NUMBER} on branch: ${env.BRANCH_NAME}.")
                                slackSend "The images will be tagged as quay.io/encode-dcc/atac-seq-pipeline:latest"
                                // pull the cache template image (the image is going to stay pretty much the same so it is no need to be dynamic)
                                sh "docker pull quay.io/encode-dcc/atac-seq-pipeline:v1.1"
                                sh "docker login -u=${QUAY_USER} -p=${QUAY_PASS} quay.io"
                                sh "docker build --cache-from quay.io/encode-dcc/atac-seq-pipeline:v1.1 -f docker_image/Dockerfile -t atac-seq-pipeline ."
                                sh "docker tag atac-seq-pipeline quay.io/encode-dcc/atac-seq-pipeline:latest"
                                sh "docker push quay.io/encode-dcc/atac-seq-pipeline:latest"
                                sh "docker logout"
                        }
                }

                stage('Run-Task-Level-Tests-Non-Master'){
                        agent {label 'slave-w-docker-cromwell-60GB-ebs'} 
                        steps{
                                sh "cd test/test_task && git clone https://github.com/ENCODE-DCC/atac-seq-pipeline-test-data"
                                sh """cd test/test_task
                                      ./download_hg38_fasta_for_test_ataqc.sh
                                      ./test.sh test_bam2ta.wdl test_bam2ta.json $TAG
                                      python -c "import sys; import json; data=json.loads(sys.stdin.read()); sys.exit(int(not data[u'match_overall']))" < test_bam2ta.result.json
                                      ./test.sh test_bowtie2.wdl test_bowtie2.json $TAG
                                      python -c "import sys; import json; data=json.loads(sys.stdin.read()); sys.exit(int(not data[u'match_overall']))" < test_bowtie2.result.json
                                      ./test.sh test_filter.wdl test_filter.json $TAG
                                      python -c "import sys; import json; data=json.loads(sys.stdin.read()); sys.exit(int(not data[u'match_overall']))" < test_filter.result.json
                                      ./test.sh test_idr.wdl test_idr.json $TAG
                                      python -c "import sys; import json; data=json.loads(sys.stdin.read()); sys.exit(int(not data[u'match_overall']))" < test_idr.result.json
                                      ./test.sh test_macs2.wdl test_macs2.json $TAG
                                      python -c "import sys; import json; data=json.loads(sys.stdin.read()); sys.exit(int(not data[u'match_overall']))" < test_macs2.result.json
                                      ./test.sh test_overlap.wdl test_overlap.json $TAG
                                      python -c "import sys; import json; data=json.loads(sys.stdin.read()); sys.exit(int(not data[u'match_overall']))" < test_overlap.result.json
                                      ./test.sh test_pool_ta.wdl test_pool_ta.json $TAG
                                      python -c "import sys; import json; data=json.loads(sys.stdin.read()); sys.exit(int(not data[u'match_overall']))" < test_pool_ta.result.json
                                      ./test.sh test_reproducibility.wdl test_reproducibility.json $TAG
                                      python -c "import sys; import json; data=json.loads(sys.stdin.read()); sys.exit(int(not data[u'match_overall']))" < test_reproducibility.result.json
                                      ./test.sh test_spr.wdl test_spr.json $TAG
                                      python -c "import sys; import json; data=json.loads(sys.stdin.read()); sys.exit(int(not data[u'match_overall']))" < test_spr.result.json
                                      ./test.sh test_trim_adapter.wdl test_trim_adapter.json $TAG
                                      python -c "import sys; import json; data=json.loads(sys.stdin.read()); sys.exit(int(not data[u'match_overall']))" < test_trim_adapter.result.json
                                      ./test.sh test_xcor.wdl test_xcor.json $TAG
                                      python -c "import sys; import json; data=json.loads(sys.stdin.read()); sys.exit(int(not data[u'match_overall']))" < test_xcor.result.json
                                      ./test.sh test_ataqc.wdl test_ataqc.json $TAG
                                      python -c "import sys; import json; data=json.loads(sys.stdin.read()); sys.exit(int(not data[u'match_overall']))" < test_ataqc.result.json
                                   """
                        }
                }

                stage('Run-Workflow-Level-Quick-Tests'){
                    agent {label 'master-builder'}
                    when {not {branch 'master'}}
                    steps {
                        echo "running subsampled/chr19_only/paired_end workflow level tests when there is an event on master branch"
                        //sh """cd test/test_workflow
                        //      ./test_atac.sh ENCSR356KRQ_subsampled_chr19_only.json $TAG
                        //      python -c "import sys; import json; data=json.loads(sys.stdin.read()); sys.exit(int(not data[u'outputs'][u'atac.qc_report.qc_json_match']))" < ENCSR356KRQ_subsampled_chr19_only.result.json
                        //   """
                    }
                }

                stage('Run-Workflow-Level-Full-Tests'){
                    agent {label 'slave-w-docker-cromwell-60GB-ebs'}
                    when { branch 'master'}
                    steps {
                        echo "running subsampled/paired_end workflow level tests when there is an event on master branch"
                        sh """
                              ./test_atac.sh ENCSR356KRQ_subsampled.json $TAG
                              python -c "import sys; import json; data=json.loads(sys.stdin.read()); sys.exit(int(not data[u'outputs'][u'atac.qc_report.qc_json_match']))" < ENCSR356KRQ_subsampled.result.json
                           """
                    }
                }
        }
                
	post {
                success {
                        echo "Post build actions that run on success"
                        slackSend "Job ${env.JOB_NAME}, build number ${env.BUILD_NUMBER} on branch ${env.BRANCH_NAME} finished with"
                        slackSend (color: '#7cfc00', message: "SUCCESS")
                        slackSend "For details, visit ${env.BUILD_URL}"
                }
                failure {
                        echo "Post build actions that run on failure"
                        slackSend "Job ${env.JOB_NAME}, build number ${env.BUILD_NUMBER} on branch ${env.BRANCH_NAME} finished with"
                        slackSend (color: '#FF0000', message: "FAILURE")
                        slackSend "For details, visit ${env.BUILD_URL}"

                }

	}
}
