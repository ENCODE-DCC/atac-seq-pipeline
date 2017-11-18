task align {
	File a
	String? b
	command {
		echo ${b} > ret.txt
		ln ${a} ${basename(a)}
	}
	output {
		File ret = basename(a)
		File? ret2 = select_first(["babo.txt",null])
	}
}

task filter {
	File a
	command {
		ln ${a} ${basename(a)}
	}
	output {
		File ret = basename(a)
	}
}

workflow test_wf {
	Array[String] fastqs
	Array[String] adapters
	String input_type
	Int num_rep

	if ( input_type in ['fastq'] ) {
		scatter( i in range(num_rep)) {
			call align {
				input: 
					a = fastqs[i],
					#b = select_first([adapters,[]])[i],
					b = adapters[i],
			}
		}
	}

	Array[String] sams
	if ( input_type=='sam' || input_type=='fastq' ) {
		scatter( i in range(num_rep)) {
			call filter {
				input: 
					a = select_first([align.ret,sams])[i]
			}
		}
	}
}
