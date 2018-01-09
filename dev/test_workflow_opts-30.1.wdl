#java -jar $(which cromwell-30.jar) run *.wdl

workflow test {
	Int num
	scatter( i in range(num) ) {
		call t0 { input: i=i }
	}
}

task t0 {
	Int? i
	command {
		echo ${i} > "out.txt"
	}
	output {
		File out = glob("out.txt")[0]
	}
}
