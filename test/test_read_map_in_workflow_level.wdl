#java -jar $(which cromwell-30.jar) run *.wdl

workflow test {
	File genome_tsv
	Map[String,String] genome = read_map(genome_tsv)

	call t0 {input: i = genome['bowtie2_idx_tar'] }
}

task t0 {
	String? i
	command {
		echo ${i}
	}
	output {
		String out = stdout()
	}
}
