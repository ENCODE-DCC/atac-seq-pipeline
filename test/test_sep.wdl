task t1 {
	Array[String] arr1
	Array[String] arr2
	command {
		echo "t1 ${sep=' ' arr1}" > "t1.txt"
		echo "t2 ${sep=' ' arr2}" > "t2.txt"
	}
	output {
		Array[File] o1 = ["t1.txt"]
		Array[File] o2 = []
	}
}

workflow test_wf {
	Array[String] arr1 = ['a','b','c']
	Array[String] arr2 = []
	call t1 { input: arr1=arr1, arr2=arr2 }
}
