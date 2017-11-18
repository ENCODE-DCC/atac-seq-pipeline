task align {
	File a
	command {
		#zcat -f ${a} | gzip -nc > (basename ${a})
		#readlink -f ${a} > tmp.txt
		ln ${a} ${basename(a)}
	}
	output {
		File ret = basename(a)
	}
	#output {
	#	File ret2 = read_string('tmp.txt')
	#}
}

workflow test_wf {
	Map[String,File] map

	scatter( val in map ) {
		call align {
			input: a = val.right
		}
	}
}
