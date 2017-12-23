workflow test_cc {
	Array[File] fastqs = ["fastqs/fq1_1-1.fastq.gz","fastqs/fq1_1-2.fastq.gz"]

	scatter(fastq in fastqs){
		call t1 {
			input:
				fastq = fastq
		}
		call t2 {
			input:
				bam = t1.bam
		}
	}
}

task t1 {
	File fastq
	command {
		zcat -f ${fastq} > ${basename(fastq) + ".bam"}
	}
	output {
		File bam = basename(fastq) + ".bam"
	}
}

task t2 {
	File bam
	command {		
		cat ${bam} > ${basename(bam) + ".bed"}
	}
	output {
		File bed = basename(bam) + ".bed"
	}
}