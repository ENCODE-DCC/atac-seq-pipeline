# ENCODE DCC ATAC-Seq/DNase-Seq pipeline tester for task trim_adapter
# Author: Jin Lee (leepc12@gmail.com)
import "../../atac.wdl" as atac

workflow test_trim_adapter {
	String cutadapt_param
	
	Array[String] pe_adapters_R1
	Array[String] pe_adapters_R2
	Array[String] pe_fastqs_R1
	Array[String] pe_fastqs_R2
	Array[String] se_adapters_R1
	Array[String] se_fastqs_R1

	Array[String] ref_pe_trimmed_fastqs
	Array[String] ref_se_trimmed_fastqs

	Int trim_adapter_cpu = 1
	Int trim_adapter_mem_mb = 12000
	Int trim_adapter_time_hr = 24
	String trim_adapter_disks = "local-disk 100 HDD"

	# pe: with adapters input, w/o auto detection
	call atac.trim_adapter as pe_trim_adapter { input :
		fastqs_R1 = pe_fastqs_R1,
		fastqs_R2 = pe_fastqs_R2,
		adapters_R1 = pe_adapters_R1,
		adapters_R2 = pe_adapters_R2,
		auto_detect_adapter = false,
		paired_end = true,
		cutadapt_param = cutadapt_param,

		cpu = trim_adapter_cpu,
		mem_mb = trim_adapter_mem_mb,
		time_hr = trim_adapter_time_hr,
		disks = trim_adapter_disks,
	}
	# pe: w/o adapters input, with auto detection
	call atac.trim_adapter as pe_trim_adapter_auto { input :
		fastqs_R1 = pe_fastqs_R1,
		fastqs_R2 = pe_fastqs_R2,
		adapters_R1 = [],
		adapters_R2 = [],
		auto_detect_adapter = true,
		paired_end = true,
		cutadapt_param = cutadapt_param,

		cpu = trim_adapter_cpu,
		mem_mb = trim_adapter_mem_mb,
		time_hr = trim_adapter_time_hr,
		disks = trim_adapter_disks,
	}
	# se: with adapters input, w/o auto detection
	call atac.trim_adapter as se_trim_adapter { input :
		fastqs_R1 = se_fastqs_R1,
		fastqs_R2 = [],
		adapters_R1 = se_adapters_R1,
		adapters_R2 = [],
		auto_detect_adapter = false,
		paired_end = false,
		cutadapt_param = cutadapt_param,

		cpu = trim_adapter_cpu,
		mem_mb = trim_adapter_mem_mb,
		time_hr = trim_adapter_time_hr,
		disks = trim_adapter_disks,
	}
	# se: w/o adapters input, with auto detection
	call atac.trim_adapter as se_trim_adapter_auto { input :
		fastqs_R1 = se_fastqs_R1,
		fastqs_R2 = [],
		adapters_R1 = [],
		adapters_R2 = [],
		auto_detect_adapter = true,
		paired_end = false,
		cutadapt_param = cutadapt_param,

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
		files = select_all([
			pe_trim_adapter.trim_merged_fastq_R1, 
			pe_trim_adapter.trim_merged_fastq_R2,
			pe_trim_adapter_auto.trim_merged_fastq_R1, 
			pe_trim_adapter_auto.trim_merged_fastq_R2,

			se_trim_adapter.trim_merged_fastq_R1, 
			se_trim_adapter_auto.trim_merged_fastq_R1, 
		]),
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
