task t1 {
	String? test
	File chrsz
	Array[Array[File]] i
	Int res_cpu

	Int? mem_mb
	Int? time_hr
	#String? queue
	command {		
		echo ${chrsz} > "test.bam"
	}
		#echo "${sep=' ' i}" "${chrsz}" "${test}" > "test.bam"
	output {
		File o1 = glob('*.bam')[0]
	}
	runtime {
		cpu : "${res_cpu}"
		memory : "${mem_mb} MB"
		time : "${time_hr}"
		#queue : "${queue}"
	}
}

task gather_outputs {
	File o1
	File o2
	command {
		ln -s ${o1}
		ln -s ${o2}
	}
}

workflow test_wf {
	#String genome_tsv
	String tsv_path 
	Map[String,String] tsv = read_map(tsv_path)
	String txt_path 
	String txt = read_string(txt_path)
	String tsv_path2
	Object? tsv_obj # = read_object(tsv_path2)
	Array[Array[Array[File]]]+ fastqs
	Int res_cpu

	scatter( i in range(length(fastqs)) ) {
		call t1 {
			input: i = fastqs[i],
			res_cpu = res_cpu/length(fastqs),
			#chrsz = tsv["blacklist"]
			chrsz = sub(tsv["chrsz"],"[\r\n]+","")
			#chrsz = tsv_obj.chrsz,
			#chrsz = txt,
		}
	}
	#call gather_outputs { input: o1 = t1.o1, o2 = test_sub_wf.o2 }
}
