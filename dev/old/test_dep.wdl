task t1 {
	command {		
		echo $(which run_spp.R) > "o.txt"
	}
	output {
		File o = "o.txt"
	}
}

workflow test_wf {
	call t1 {
	}
}
