workflow test {
	call t1 {}
}

task t1 {
	Boolean paired_end = true
	#Boolean paired_end = false
	command {
		echo 1234 > "bout.R1.fastq.gz"
		echo 1234 > "aout.R2.fastq.gz"
	}
	output {
		#Array[File] a = glob("*.R?.fastq.gz")
		Array[File] a = if paired_end then 
			[glob("*.R1.fastq.gz")[0], glob("*.R2.fastq.gz")[0]]
			else glob("*.R1.fastq.gz")
	}
}
