# not working
#java -jar $(which cromwell-29.jar) run test_nested_if_in_scatter_in_cromwell-30.wdl
# working
#java -jar $(which cromwell-30.jar) run test_nested_if_in_scatter_in_cromwell-30.wdl

workflow test {
	Boolean b0 = true
	Boolean b1 = true
	Boolean b2 = true
	if ( b0 ) {
		scatter( i in range(2) ) {
			if ( b1 ) {
				if ( b2 ) {
					call t0 {}
				}
			}
		}
	}
	call t1 { input: arr=t0.out }
}

task t0 {
	command {
		echo test > out.txt
	}
	output {
		File out = glob('out.txt')[0]
	}
}

task t1 {
	Array[File?] arr
	command {
		echo ${sep=' ' arr}
	}
	output {
		String out = stdout()
	}
}
