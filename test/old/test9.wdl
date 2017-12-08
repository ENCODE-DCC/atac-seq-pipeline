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
	String? a
	command {
		echo "t1:${a}"
	}
	output {
		String b = read_string(stdout())
	}
}

task t2 {
	String? p1
	String? p2
	command {
		echo "p1:${p1},p2:${p2}"
	}
	output {
		String b = read_string(stdout())
	}		
}

workflow test_wf {
	Array[Array[String]] fastqs = [ ["fastqs/fq1_1-1.fastq.gz", "fastqs/fq1_1-2.fastq.gz"], ["fq1_2-1.fastq.gz", "fq1_2-2.fastq.gz"] ]
	Array[File] bams = ["fastqs/fq1_1-1.fastq.gz", "fastqs/fq1_1-2.fastq.gz"]
	call inputs { input : fastqs = fastqs, bams = bams }
	call pair_gen { input : num_rep = inputs.num_rep }

	scatter( i in range(inputs.num_rep) ) {
		call t1 { input : a = bams[i] }
	}

	scatter( pair in pair_gen.pairs ) {
		call t2 { input : p1 = bams[ pair[0] ], p2 = pair[1] }
	}

	#scatter( pair in pair_gen.pairs ) {
	#	call t2 { input : peaks = bams, i1 = pair[0], i2 = pair[1] }
	#}
	#if ( 1>2 ) {
	#		call t1 { input : a="merong" }
	#}
	#call t1 as t1_ { input : a=t1.b }
}
