# not working
#java -jar $(which cromwell-29.jar) run test_select_first_bug_in_task_output_FIXED_IN_CROMWELL-29.wdl
# working
#java -jar $(which cromwell-30.jar) run test_select_first_bug_in_task_output_FIXED_IN_CROMWELL-30.wdl

workflow test {
	call filter { input:
		no_dup_removal = true,
	}
}

task filter {
	# parameters from workflow
	Boolean? no_dup_removal 	# no dupe reads removal when filtering BAM
								# dup.qc and pbc.qc will be emptry files
								# and nodup_bam in the output is 
	command {
		${if select_first([no_dup_removal,false]) then ""
			else "echo temp > null.dup.qc"}
		touch null
	}
	output {
		File dup_qc = if select_first([no_dup_removal,false]) then glob("null")[0] else glob("*.dup.qc")[0]
	}
}

