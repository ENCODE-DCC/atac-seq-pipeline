task align {
	File tsv
	Array[Array[File]] fastqs = read_tsv(tsv)
	command {
		zcat -f ${sep=' ' fastqs} > "test.sam"
		#zcat -f ${fastqs} > "test.sam"
	}
	output {
		File sam = glob('*.sam')[0]
	}
}

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
	#Array[Array[File]] fastqs
	String input_type
	Array[Array[Array[String]]?]? input_fastqs

	if ( input_type=='fastq' ) {
		scatter( input_list in input_fastqs ) {
			File input_files_tsv = write_tsv(input_fastqs)
			call align {
				input: tsv = input_list
			}
		}
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
