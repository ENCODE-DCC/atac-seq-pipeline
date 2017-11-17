task pair_gen { # returns every pair of true replicate
	Int num_rep
	command {
		echo ${num_rep}
	}
	output {
		String pairs = stdout()
	}	
}

workflow test_wf {
	Int num_rep = 4
	call pair_gen { input : num_rep = num_rep }
}
