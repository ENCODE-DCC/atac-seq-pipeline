pipeline {
        agent none
        environment {
                QUAY_USER = credentials('quay-robot')
                QUAY_PASS = credentials('quay-robot-token')
        }
        stages {
		stage('Unit-tests') {
                        agent {label 'master-builder'}
			steps { 
                                // the tag gets built here, and can be referenced in the other stages 
                                script {
                                        TAG = sh([script: "echo quay.io/encode-dcc/atac-seq-pipeline:${env.BRANCH_NAME}_${env.BUILD_NUMBER}", returnStdout: true]).trim()
                                       }
				echo "Running unit tests.."
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
                                sh "docker pull quay.io/encode-dcc/atac-seq-pipeline:develop_test_jenkins_30"
                                sh "docker login -u=${QUAY_USER} -p=${QUAY_PASS} quay.io"
                                sh "docker build --cache-from quay.io/encode-dcc/atac-seq-pipeline:develop_test_jenkins_30 -f docker_image/Dockerfile -t atac-seq-pipeline ."
                                sh "docker tag atac-seq-pipeline $TAG"
                                sh "docker push $TAG"
                                sh "docker logout"
                                
                        }
                }
                stage('Build-master') {
                        agent {label 'master-builder'} //this will be the slave-w-docker-cromwell-60GB-ebs
                        when { branch 'master'}
                        steps {
                                echo "the tag is $TAG"
				                slackSend "started job: ${env.JOB_NAME}, build number ${env.BUILD_NUMBER} on branch: ${env.BRANCH_NAME}."
				                slackSend "The images will be tagged as $TAG"
                                echo "Running master build steps."
                        }
                }
                stage('Run-Task-Level-Tests'){
                        agent {label 'slave-w-docker-cromwell-60GB-ebs'} 
                        steps{
                                echo "running task level tests on every push on every branch"
                                sh "cd test/test_task && git clone https://github.com/leepc12/atac-seq-pipeline-test-data"
                                sh """cd test/test_task
                                      ./test.sh test_bam2ta.wdl test_bam2ta.json $TAG
                                   """
                        }
                }
                stage('Run-Workflow-Level-Tests'){
                    agent {label 'master-builder'}
                    when {branch 'master'}
                    steps{
                        echo "running long workflow level tests when there is an event on master branch"
                        //add the test run here
                    }
                }
        }
	post {
                success {
                        echo "Post build actions that run on success"
                        slackSend "Job ${env.JOB_NAME}, build number ${env.BUILD_NUMBER} on branch ${env.BRANCH_NAME} finished with"
                        slackSend (color: '#7cfc00', message: "SUCCESS")
                }
                failure {
                        echo "Post build actions that run on failure"
                        slackSend "Job ${env.JOB_NAME}, build number ${env.BUILD_NUMBER} on branch ${env.BRANCH_NAME} finished with"
                        slackSend (color: '#FF0000', message: "FAILURE")
                }
                always {
                        echo "Post build actions that run always"
                }
	}
}
