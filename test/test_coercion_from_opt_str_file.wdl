workflow test {
	#String test_file = "/dev/null"
	#String genome_tsv = "/users/leepc12/code/atac-seq-pipeline/test/mm10_klab.tsv"
	String genome_tsv = "gs://atac-seq-pipeline-genome-data/mm10_google.tsv"
	call inputs { input: genome_tsv = genome_tsv }
	#call t0 { input: test_file = test_file, make_signal = true, blacklist= if inputs.has_blacklist then [inputs.blacklist] else [] }
	call t0 { input: make_signal = true, blacklist= [] }
	call t1 { input: t0_out=t0.out, t0_out_bfilt=t0.out_bfilt }
}

task t0 {
	#File test_file
	Array[File?] blacklist
	Boolean? make_signal
	String peak_type = "narrowPeak"
	command {
		echo 123 > "tto_xyz.123.narrowPeak.gz"
		echo 123 > "tto_xyz.123.bfilt.narrowPeak.gz"
		touch null
		echo 123 > "xx.pval.signal.bigwig"
		echo 123 > "xx.fc.signal.bigwig"
	}
	output {
		File out = glob("*[!.][!f][!i][!l][!t]."+peak_type+".gz")[0]
		File out_bfilt = if length(blacklist)>0 then
				glob("*.bfilt."+peak_type+".gz")[0] else out
		Boolean make_signal_ = select_first([make_signal,false])
		File sig_pval = if make_signal_ then glob("*.pval.signal.bigwig")[0] else glob("null")[0]
		File sig_fc = if make_signal_ then glob("*.fc.signal.bigwig")[0] else glob("null")[0]
		#File sig_pval = if (select_first([make_signal,false])) then glob("*.pval.signal.bigwig")[0] else glob("null")[0]
		#File sig_fc = if (select_first([make_signal,false])) then glob("*.fc.signal.bigwig")[0] else glob("null")[0]
		#File sig_pval = glob("*.pval.signal.bigwig")[0]
		#File sig_fc = glob("*.fc.signal.bigwig")[0]
		#File sig_pval = glob("null")[0]
		#File sig_fc = glob("null")[0]
	}
}

task t1 {
	File t0_out
	File t0_out_bfilt
	command {
		echo
	}
}

task inputs {
	File genome_tsv
	command <<<
		touch null
	>>>
	output {
		Map[String,String] genome = read_map(genome_tsv)
		String ref_fa = genome['ref_fa']
		String bowtie2_idx_tar = genome['bowtie2_idx_tar']
		#String bwa_idx_tar = genome['bwa_idx_tar']
		#String blacklist = "/dev/null"
		String blacklist = genome['blacklist']
		String chrsz = genome['chrsz']
		String gensz = genome['gensz']

		Boolean has_blacklist = blacklist!="/dev/null"
	}
}
 
