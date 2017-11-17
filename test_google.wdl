task pair_gen { # returns every pair of true replicate
	Int num_rep
	command {
	}
	output {
		String pairs = num_rep
	}	
}

workflow test_wf {
	Int num_rep = 4
	call pair_gen { input : num_rep = num_rep }
}
