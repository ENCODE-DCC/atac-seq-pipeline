workflow test {
	String? null
	File genome_tsv = '/users/leepc12/code/atac-seq-pipeline/genome/hg38_klab.tsv'
	#File genome_tsv = '/users/leepc12/code/atac-seq-pipeline/genome/hg38_google.tsv'
	#File genome_tsv = 'hg38_klab.tsv'
	#File genome_tsv = 'hg38_google.tsv'
	#String genome_tsv #= 'gs://atac-seq-pipeline-genome-data/hg38_google.tsv'
	call inputs { input: genome_tsv = genome_tsv }
	String? blacklist = if inputs.genome['blacklist']=='/dev/null' then null else inputs.genome['blacklist']
	String a = select_first([inputs.genome['test'],'/dev/null'])

	call t0 { input: blacklist = blacklist, chrsz=inputs.genome['chrsz'] }
}

task t0 {
	File? blacklist
	File chrsz
	command {
		echo ${"--blacklist " + blacklist} > 'b1.txt'
		echo ${"--chrsz " + chrsz} > 'b2.txt'
		touch null
	}
	output {
		File blacklist_out = if defined(blacklist) then globe("b1.txt")[0] else glob("null")[0]
	}
}

task inputs {
	File genome_tsv
	command {
		echo
	}
	output {
		Map[String,String] genome = read_map(genome_tsv)
	}
}
