import "wdl/test3_sub.wdl" as sub

task filter {
	File sam
	command {
		zcat -f ${sam} | gzip -nc > "test.bam"
	}
	output {
		File bam = glob('*.bam')[0]
	}
}

workflow test_wf {
	String input_type

	if ( input_type=='fastq' ) {
		call sub.test_sub_wf {}
	}

	Array[File] input_sams

	if ( input_type=='sam' ) {
		scatter( sam in input_sams ) {
			call filter {
				input: sam = sam
			}
		}
	}
}
