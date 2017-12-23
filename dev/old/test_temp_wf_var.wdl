workflow test {
	String i1 = 'i1'

	call t1 { input : i1 = i1 }

	if ( false ) {
		call t1 as t0 { input : i1 = i1 }
	}

	String i2 = t1.out

	call t1 as t2 { input : i1 = i2 }
}

task t1 {
	String i1
	command {
		echo "${i1},babo" > 'txt.txt'
	}
	output {
		String out = read_string('txt.txt')
	}
}
