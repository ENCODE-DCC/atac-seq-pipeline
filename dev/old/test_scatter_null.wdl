task t1 {
	command {		
		echo 1234 > "test.bam"
	}
	output {
		File o1 = glob('*.bam')[0]
	}
}

workflow test_wf {
	scatter( i in range(0) ) {
		call t1 {
		}
	}
}
