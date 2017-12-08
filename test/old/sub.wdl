task t2 {
	String var2
	command {
		echo ${var2} > 'o2.txt'
	}
	output {
		File o2 = 'o2.txt'
	}
}

workflow test_sub_wf {
	String var2

	call t2 { input: var2 = var2 }
	output {
		File o2 = t2.o2
	}
}

