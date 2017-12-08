task align {
	Array[Array[File]] fastqs
	Array[Array[String]] adapters
	command {		
		echo "${write_tsv(adapters)} ${sep=' ' fastqs[0]}" > ret.txt
	}
	output {
		File ret = "ret.txt"
	}
}

workflow test_wf {
	Array[Array[Array[String]]] fastqs
	Array[Array[Array[String]]] adapters

	scatter(i in range(length(fastqs))) {
		call align {
			input: 
				fastqs = fastqs[i],
				adapters = if (length(adapters)>0) 
							then adapters[i] else [],
				#adapters = adapters[i],
		}
	}
}
