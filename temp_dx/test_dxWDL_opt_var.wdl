workflow test_dxWDL_opt_var {
	Array[File] fastqs 
	Int num_rep = length(fastqs)

	scatter(i in range(num_rep)) {
		call get_num_line { input :
			fastq = fastqs[i],
		}
	}
}

task get_num_line {
	File fastq

	command {		
		zcat ${fastq} | wc -l
	}
	output {
		String out = read_string(stdout())
	}
}
