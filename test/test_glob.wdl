workflow test {
	call t1 {}
}

task t1 {
	command {
		echo 1234 > "out.R2.fastq.gz"
		echo 1234 > "out.R1.fastq.gz"
	}
	output {
		Array[File] a = glob("*.R?.fastq.gz")
	}
}
