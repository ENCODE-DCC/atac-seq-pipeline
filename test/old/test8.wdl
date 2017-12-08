task pair_gen { # returns every pair of true replicate
	Int num_rep
	command {
		for i in $(seq 0 $((${num_rep}-1)));
		do
		  for j in $(seq $((i+1)) $((${num_rep}-1)));
		    do
		    echo -e "$i\t$j"
		  done
		done
		sleep 10
	}
	output {
		Array[Array[Int]] pairs = read_tsv(stdout())
	}
	runtime {
		cpu : 2
	}
}

workflow test_wf {
	Int num_rep = 4
	call pair_gen { input : num_rep = num_rep }
}
