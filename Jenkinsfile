pipeline {
        agent {label 'slave-w-docker-cromwell-60GB-ebs'}

        stages {
		stage('Unit-tests') {
			steps { 
				echo "Running unit tests.."
			}
		}
                stage('Build-nonmaster') {
                        when { not { branch 'master' } }
                        steps { 
                                slackSend "started job: ${env.JOB_NAME}, build number ${env.BUILD_NUMBER} on branch: ${env.BRANCH_NAME}."
				slackSend "The images will be tagged as ${env.BRANCH_NAME}:${env.BUILD_NUMBER}"
                                echo "${env.BRANCH_NAME}"
                                echo "Running non-master build steps."
                        }
                }
                stage('Build-master') {
                        when { branch 'master'}
                        steps {
				slackSend "started job: ${env.JOB_NAME}, build number ${env.BUILD_NUMBER} on branch: ${env.BRANCH_NAME}."
				slackSend "The images will be tagged as ${env.BRANCH_NAME}:${env.BUILD_NUMBER}"
                                echo "${env.BRANCH_NAME}"
                                echo "Running master build steps."
                        }
                }
        }
	post {
                success {
                        echo 'Post build actions that run on success'
                        slackSend "Job ${env.JOB_NAME}, build number ${env.BUILD_NUMBER} on branch ${env.BRANCH_NAME} finished successfully."
                }
                failure {
                        echo 'Post build actions that run on failure'
                        slackSend "Job ${env.JOB_NAME}, build number ${env.BUILD_NUMBER} on branch ${env.BRANCH_NAME} failed."
                }
                always {
                        echo 'Post build actions that run always'
                }
	}
}
