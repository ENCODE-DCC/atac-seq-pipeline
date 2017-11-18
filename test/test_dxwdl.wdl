task inputs {
	Array[Array[Array[String]]] fastqs
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
	File? a
	command {
		zcat -f ${a}
	}
	output {
		String b = read_string(stdout())
	}
}

task t2 {
	File? p1
	File? p2
	command {
		zcat -f ${p1} ${p2}
	}
	output {
		String b = read_string(stdout())
	}		
}

workflow test_wf {
	Array[Array[Array[String]]] fastqs # = [ [["fastqs/fq1_1-1.fastq.gz"], ["fastqs/fq1_1-2.fastq.gz"]], [["fq1_2-1.fastq.gz"], ["fq1_2-2.fastq.gz"]] ]
	Array[String] bams # = ["fastqs/fq1_1-1.fastq.gz", "fastqs/fq1_1-2.fastq.gz"]
	String? hello
	Boolean? paired_end
	call inputs { input : fastqs = fastqs, bams = bams }
	call pair_gen { input : num_rep = inputs.num_rep }

	if ( select_first([paired_end,false]) ) {
		#scatter( i in range(inputs.num_rep) ) {
			call t1 { input : a = bams[0] }
		#}

		#scatter( i in range(pair_gen.pairs) ) {
	#		if defined(hello) {
			call t2 { 
				input : 
					#p1 = bams[(pair_gen.pairs[i][0])], 
					#p2 = bams[(pair_gen.pairs[i][0])],
					p1 = fastqs[0][0][0], 
					p2 = fastqs[1][0][0], 
			}
		#}
	}

	#scatter( pair in pair_gen.pairs ) {
	#	call t2 { input : peaks = bams, i1 = pair[0], i2 = pair[1] }
	#}
	#if ( 1>2 ) {
	#		call t1 { input : a="merong" }
	#}
	#call t1 as t1_ { input : a=t1.b }
}
