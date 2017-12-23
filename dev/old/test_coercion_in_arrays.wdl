task t1 {
	command {
		echo 1 > "t1.txt"
		echo 2 >> "t1.txt"
	}
	output {
		Array[Int] idx = read_lines("t1.txt")
	}
}

workflow test_wf {
	call t1 {}
}
