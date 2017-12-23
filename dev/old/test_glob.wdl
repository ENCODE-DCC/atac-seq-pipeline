workflow test {
	call t1 {}
}

task t1 {
	Boolean paired_end = true
	#Boolean paired_end = false
	command {
		echo 1234 > "R1_bout.R1.fastq.gz"
		echo 1234 > "R2_aout.R2.fastq.gz"
	}
	output {
		Array[File] a = glob("R?_*.fastq.gz")
	}
}
