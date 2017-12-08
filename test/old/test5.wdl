task align {
	Array[Array[File]] fastqs
	command {
		echo "${sep=' ' fastqs}" > "test.sam"
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
	Array[Array[Array[String]]] fastqs

	if ( input_type=='fastq' ) {
		scatter( list in fastqs ) {
			call align {
				input: fastqs = list
			}
		}
	}
	Array[String] sams
	if (input_type=='fastq') {
		sams = align.sam
	}
	scatter( sam in sams ) {
 		call filter {
			input: sam = sam
		}
	}
}
