workflow test {
	#Array[Array[Array[String]]] fastqs = [ [  ['1-1-1','1-1-2'] ], [  ['2-1-1','2-1-2'] ] ]	
	Array[Array[Array[String]]] fastqs = [ [  [] ], [  [] ] ]	
	call t1 { input : fastqs = fastqs }
}

task t1 {
	Array[Array[Array[String]]] fastqs
	command <<<
		python <<CODE
		import json
		print( json.loads('[${sep=',' fastqs}]') )
		CODE
	>>>
	output {
		String out = stdout()
	}
}
