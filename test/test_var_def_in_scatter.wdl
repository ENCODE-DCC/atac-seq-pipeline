#java -jar $(which cromwell-30.jar) run *.wdl

workflow test {
	scatter(i in range(3)) {
		Int ctl_ta = i
	}
	output {
		ctl_ta
	}
}

