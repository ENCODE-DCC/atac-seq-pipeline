pipeline {
        agent none
        stages {
		stage('Unit-tests') {
                        agent {label 'master-builder'}
			steps { 
				echo "Running unit tests.."
			}
		}
                stage('Build-nonmaster') {
                        agent {label 'master-builder'} //this will happen on slave in the real thing
                        when { not { branch 'master' } }
                        steps { 
                                slackSend "started job: ${env.JOB_NAME}, build number ${env.BUILD_NUMBER} on branch: ${env.BRANCH_NAME}."
				slackSend "The images will be tagged as ${env.BRANCH_NAME}:${env.BUILD_NUMBER}"
                                echo "${env.BRANCH_NAME}"
                                echo "Running non-master build steps."
                        }
                }
                stage('Build-master') {
                        agent {label 'master-builder'} //this will happen on slave in the real thing
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
