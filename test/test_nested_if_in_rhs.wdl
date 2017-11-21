workflow test {
	Boolean paired_end = false
	Boolean true_rep_only = false
	call t1 {
		input : 
			a = if paired_end then "p1"
				else if true_rep_only then "p2"
				else "p3"
	}
}

task t1 {
	String? a
	command {
		echo ${a} > "R1_bout.R1.fastq.gz"
	}
	output {
		Array[File] out = glob("R?_*.fastq.gz")
	}
}
