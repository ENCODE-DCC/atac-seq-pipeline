#java -jar $(which cromwell-30.jar) run *.wdl

workflow test {
	Boolean b0
	Boolean b1
	Boolean b2
	scatter( i in range(3) ) {
		if ( b0 ) {
			call t0 as t1 { input: i=i }
		}
	}
	if ( b1 ) {
		scatter( i in range(3) ) {
			call t0 as t2 { input: i=t1.out[i] }
		}
	}
	if ( b1 && b2 ) {
		scatter( i in range(3) ) {
			call t0 as t3 { input: i=t1.out[i] }
		}
	}
}

task t0 {
	Int? i
	command {
		echo ${i}
	}
	output {
		Int out = read_int(stdout())
	}
}
