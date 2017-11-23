workflow test {
	File? tmp
	call t0 {}
	call t1 { input: tmp = tmp, t = t0.out }
}

task t1 {
	File? tmp
	File? t
	#Boolean paired_end = true
	Boolean paired_end = false
	command {
		echo 1234 > "R1_bout.R1.fastq.gz"
		echo 1234 > "R2_aout.R2.fastq.gz"
	}
	output {
		Array[File] a = glob("R?_*.fastq.gz")
		Array[File] b = if paired_end then glob("R?_*.fastq.gz") else []
	}
}

task t0 {
	command {
		echo 1234 > "X.txt"
	}
	output {
		File out = "X.txt"
	}
}
