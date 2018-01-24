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
                        agent {label 'master-builder'} //this will be the slave-w-docker-cromwell-60GB-ebs 
                        when { not { branch 'master' } }
                        steps { 
                                echo "the tag is $TAG"
                                slackSend "started job: ${env.JOB_NAME}, build number ${env.BUILD_NUMBER} on branch: ${env.BRANCH_NAME}."
				slackSend "The images will be tagged as $TAG"
                                sh 'docker login -u=${QUAY_USER} -p=${QUAY_PASS} quay.io'
                                sh 'docker logout'
                                
                        }
                }
                stage('Build-master') {
                        agent {label 'slave-w-docker-cromwell-60GB-ebs'} //this will be the slave-w-docker-cromwell-60GB-ebs
                        when { branch 'master'}
                        steps {
                                echo "the tag is $TAG"
				slackSend "started job: ${env.JOB_NAME}, build number ${env.BUILD_NUMBER} on branch: ${env.BRANCH_NAME}."
				slackSend "The images will be tagged as $TAG"
                                echo "Running master build steps."
                        }
                }
                stage('Run-Cromwell-Tests'){
                        agent {label 'master-builder'} //this will actually run on master
                        steps{
                                echo "run cromwell tests"
                                //sh "cd test && ./test_atac.sh"
                                
                                
                        }
                }
        }
	post {
                success {
                        echo "Post build actions that run on success"
                        slackSend "Job ${env.JOB_NAME}, build number ${env.BUILD_NUMBER} on branch ${env.BRANCH_NAME} finished successfully."
                }
                failure {
                        echo "Post build actions that run on failure"
                        slackSend "Job ${env.JOB_NAME}, build number ${env.BUILD_NUMBER} on branch ${env.BRANCH_NAME} failed."
                }
                always {
                        echo "Post build actions that run always"
                }
	}
}
