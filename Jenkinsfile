pipeline {
        agent none
        environment {
                QUAY-USER = credentials('quay-robot')
                QUAY-PASS = credentials('quay-robot-token')
        }
        stages {
		stage('Unit-tests') {
                        agent {label 'master-builder'}
			steps { 
				echo "Running unit tests.."
			}
		}
                stage('Build-nonmaster') {
                        agent {label 'slave-w-docker-cromwell-60GB-ebs'}
                        when { not { branch 'master' } }
                        steps { 
                                slackSend "started job: ${env.JOB_NAME}, build number ${env.BUILD_NUMBER} on branch: ${env.BRANCH_NAME}."
				slackSend "The images will be tagged as ${env.BRANCH_NAME}:${env.BUILD_NUMBER}"
                                sh 'docker login -u=${QUAY-USER} -p=${QUAY-PASS} quay.io'
                                sh 'docker logout'
                                
                        }
                }
                stage('Build-master') {
                        agent {label 'slave-w-docker-cromwell-60GB-ebs'} 
                        when { branch 'master'}
                        steps {
				slackSend "started job: ${env.JOB_NAME}, build number ${env.BUILD_NUMBER} on branch: ${env.BRANCH_NAME}."
				slackSend "The images will be tagged as ${env.BRANCH_NAME}:${env.BUILD_NUMBER}"
                                echo "${env.BRANCH_NAME}"
                                echo "Running master build steps."
                        }
                }
                stage('Run-Cromwell-Tests'){
                        agent {label 'master-builder'} //this will actually run on master
                        steps{
                                echo "run cromwell tests"
                                //sh "cd test && ./test_atac.sh"
                                sh "ls"
                                
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
