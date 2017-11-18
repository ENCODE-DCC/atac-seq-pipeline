task inputs {
	Array[Array[String]] fastqs
	Array[String] bams
	command <<<
		python <<CODE
		name = ['fastq','bam']
		arr = [${length(fastqs)},
		       ${length(bams)}]
		num_rep = max(arr)
		type = name[arr.index(num_rep)]
		with open('num_rep.txt','w') as fp:
		    fp.write(str(num_rep)) 
		with open('type.txt','w') as fp:
		    fp.write(type)		    
		CODE
	>>>
	output {
		String type = read_string("type.txt")
		Int num_rep = read_int("num_rep.txt")
	}
}

task pair_gen { # returns every pair of true replicate
	Int num_rep
	command <<<
		python <<CODE
		for i in range(${num_rep}):
		    for j in range(i+1,${num_rep}):
		        print('{}\t{}'.format(i,j))
		CODE
	>>>
	output {
		Array[Array[Int]] pairs = read_tsv(stdout())
	}	
}

task t1 {
	File a
	command {
		zcat -f ${a} > "t1.txt"
	}
	output {
		File ret = "t1.txt"
	}
}

task t2 {
	Array[File] tas
	command {
		zcat -f "${sep=' ' tas}" > "t2.txt"
	}
	output {
		File ret = "t2.txt"
	}		
}

workflow test_wf {
	Array[Array[String]] fastqs = [] # [ ["fastqs/fq1_1-1.fastq.gz", "fastqs/fq1_1-2.fastq.gz"], ["fq1_2-1.fastq.gz", "fq1_2-2.fastq.gz"] ]
	Array[String] bams = ["fastqs/fq1_1-1.fastq.gz", "fastqs/fq1_1-2.fastq.gz"]
	call inputs { input : fastqs = fastqs, bams = bams }
	call pair_gen { input : num_rep = inputs.num_rep }

	scatter( i in range(inputs.num_rep) ) {
		call t1 { input: a = bams[i] }
	}
	call t2 { input: tas = t1.ret }
}
