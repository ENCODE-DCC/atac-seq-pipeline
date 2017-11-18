task align {
	Boolean? make_signal
	command {		
		echo ${make_signal} > test1.txt
		echo ${make_signal} > test2.txt
	}
	output {
		Array[File] ret = []
	}
}

task filter {
	File f
	command {		
		cat ${f} > "test3.txt"
	}
	output {
		String ret = read_string("test3.txt")
	}
}

workflow test_wf {
	call align as align1 {
	}
	call align as align2 {
		input:
			make_signal = align1.input.make_signal,
	}
	#if ( defined(align.ret) ) {
	#	call filter {
	#		input:
	#			f = align.ret,
	#	}
	#}	
}
