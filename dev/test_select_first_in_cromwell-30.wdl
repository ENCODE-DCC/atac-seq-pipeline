#java -jar $(which cromwell-30.jar) run *.wdl

workflow test {
	Array[String]? fastqs
	Array[String]? adapters
	Array[String]? bams
	Array[String]? tas

	Array[String] fastqs_ = select_first([fastqs,[]])
	Array[String] adapters_ = select_first([adapters,[]])
	Int fastqs_len = length(fastqs_)
	if ( fastqs_len>0 ) {
		scatter(i in range(fastqs_len)) {
			call t0 as align { input: i = fastqs_[i], j = adapters_[i] }
		}
	}
	Array[String] bams_ = select_first([bams,align.out,[]])
	if ( length(bams_)>0 ) {
		scatter(bam in bams_) {
			call t0 as bam2ta { input: i = bam }
		}
	}
	Array[String] tas_ = select_first([tas,bam2ta.out,[]])
	if ( length(tas_)>0 ) {
		scatter(ta in tas_) {
			call t0 as macs2 { input: i = ta }
		}
	}
}

task t0 {
	String? i
	String? j
	command {
		echo ${i}${','+j}
	}
	output {
		String out = read_string(stdout())
	}
}
