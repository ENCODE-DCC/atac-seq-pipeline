# ENCODE DCC ATAC-Seq/DNase-Seq pipeline tester for task trim_adapter
# Author: Jin Lee (leepc12@gmail.com)
import "../../atac.wdl" as atac

workflow test_trim_adapter {
	Int cutadapt_min_trim_len = 5	# minimum trim length for cutadapt -m
	Float cutadapt_err_rate = 0.1	# Maximum allowed adapter error rate for cutadapt -e	
	
	Array[Array[String]] pe_adapters
	Array[Array[String]] pe_fastqs
	Array[Array[String]] se_adapters
	Array[Array[String]] se_fastqs

	Array[String] ref_pe_trimmed_fastqs
	Array[String] ref_se_trimmed_fastqs

	Int trim_adapter_cpu = 1
	Int trim_adapter_mem_mb = 12000
	Int trim_adapter_time_hr = 24
	String trim_adapter_disks = "local-disk 100 HDD"

	# pe: with adapters input, w/o auto detection
	call atac.trim_adapter as pe_trim_adapter { input :
		fastqs = pe_fastqs,
		adapters = pe_adapters,
		auto_detect_adapter = false,
		paired_end = true,
		min_trim_len = cutadapt_min_trim_len,
		err_rate = cutadapt_err_rate,

		cpu = trim_adapter_cpu,
		mem_mb = trim_adapter_mem_mb,
		time_hr = trim_adapter_time_hr,
		disks = trim_adapter_disks,
	}
	# pe: w/o adapters input, with auto detection
	call atac.trim_adapter as pe_trim_adapter_auto { input :
		fastqs = pe_fastqs,
		adapters = [],
		auto_detect_adapter = true,
		paired_end = true,
		min_trim_len = cutadapt_min_trim_len,
		err_rate = cutadapt_err_rate,

		cpu = trim_adapter_cpu,
		mem_mb = trim_adapter_mem_mb,
		time_hr = trim_adapter_time_hr,
		disks = trim_adapter_disks,
	}
	# se: with adapters input, w/o auto detection
	call atac.trim_adapter as se_trim_adapter { input :
		fastqs = se_fastqs,
		adapters = se_adapters,
		auto_detect_adapter = false,
		paired_end = false,
		min_trim_len = cutadapt_min_trim_len,
		err_rate = cutadapt_err_rate,

		cpu = trim_adapter_cpu,
		mem_mb = trim_adapter_mem_mb,
		time_hr = trim_adapter_time_hr,
		disks = trim_adapter_disks,
	}
	# se: w/o adapters input, with auto detection
	call atac.trim_adapter as se_trim_adapter_auto { input :
		fastqs = se_fastqs,
		adapters = [],
		auto_detect_adapter = true,
		paired_end = false,
		min_trim_len = cutadapt_min_trim_len,
		err_rate = cutadapt_err_rate,

		cpu = trim_adapter_cpu,
		mem_mb = trim_adapter_mem_mb,
		time_hr = trim_adapter_time_hr,
		disks = trim_adapter_disks,
	}

	call atac.compare_md5sum { input :
		labels = [
			'pe_trim_adapter_R1',
			'pe_trim_adapter_R2',
			'pe_trim_adapter_auto_R1',
			'pe_trim_adapter_auto_R2',

			'se_trim_adapter',
			'se_trim_adapter_auto',
		],
		files = [
			pe_trim_adapter.trimmed_merged_fastqs[0], 
			pe_trim_adapter.trimmed_merged_fastqs[1],
			pe_trim_adapter_auto.trimmed_merged_fastqs[0], 
			pe_trim_adapter_auto.trimmed_merged_fastqs[1],

			se_trim_adapter.trimmed_merged_fastqs[0], 
			se_trim_adapter_auto.trimmed_merged_fastqs[0], 
		],
		ref_files = [
			ref_pe_trimmed_fastqs[0],
			ref_pe_trimmed_fastqs[1],
			ref_pe_trimmed_fastqs[0],
			ref_pe_trimmed_fastqs[1],

			ref_se_trimmed_fastqs[0],
			ref_se_trimmed_fastqs[0],
		],
	}
}
