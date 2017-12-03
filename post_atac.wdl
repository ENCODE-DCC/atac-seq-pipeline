# ENCODE DCC Post processing for metadata output JSON
# Author: Jin Lee (leepc12@gmail.com)

workflow post_atac {
	Array[File] meta_out_jsons
	scatter( meta_out_json in meta_out_json ) {
		call make_input_json_for_resuming { input:
			meta_out_json = meta_out_json,
		}
	}
	call make_task_dependency_tree { input:
		meta_out_jsons = meta_out_jsons,
	}
	call final_qc_report { input:
		meta_out_jsons = meta_out_jsons,
	}
}

task make_input_json_for_resuming {
	Array[File] meta_out_jsons

	command {
		python $(which encode_make_input_json.py)
	}
	output {

	}
}

task make_qc_report {
	Array[File] meta_out_jsons

	command {
		python $(which encode_parse_meta_out_jsons.py)
	}
	output {

	}
}
