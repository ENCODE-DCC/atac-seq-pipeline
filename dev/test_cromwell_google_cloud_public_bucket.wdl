# cromwell-30.jar doesn't seem to work with Google Cloud requester pays bucket

workflow test {
	String f = "gs://encode-pipeline-genome-data/test.txt"
	call t0 { input: f=f }
}

task t0 {
	File f
	command {
		zcat -f ${f}
	}
	output {
		String out = read_string(stdout())
	}
}
