workflow test {
	call t1 { input : flag1 = false, }
}

task t1 {
	Boolean? flag1
	command {
		echo test1 > test1.txt
		echo test2 > test2.txt
	}
	output {
		#File out_false = if false then glob('test1.txt')[0] else glob('test2.txt')[0]
		File out = if select_first([flag1,false]) then glob('test1.txt')[0] else glob('test2.txt')[0]
	}
}
