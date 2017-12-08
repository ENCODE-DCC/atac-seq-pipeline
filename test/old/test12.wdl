task align {
	String a
	command {		
		echo ${a} > out.txt
	}
	output {
		String ret = read_string("out.txt")
	}
}

task filter {
	String a
	command {		
		echo ${a} > out.txt
	}
	output {
		String ret = read_string("out.txt")
	}
}

workflow test_wf {
	Array[String] arr = ["1","2","3","4"]
	Array[String] arr2 = ["a","b","c","d"]

	#String input_type = "1"
	String input_type = "2"

	scatter( i in range(length(arr)) ) {
		if ( input_type=="1" ) {
			call align { 
				input : a = arr[i]
			}
		}
		if ( input_type=="1" || input_type=="2" ) {
			call filter {
				input : a = if defined(align.ret) then align.ret else arr2[i]
			}
		}
	}
}
