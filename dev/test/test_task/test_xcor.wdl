# ENCODE DCC ATAC-Seq/DNase-Seq pipeline tester
# Author: Jin Lee (leepc12@gmail.com)
import '../../../atac.wdl' as atac
import 'compare_md5sum.wdl' as compare_md5sum

workflow test_xcor {
	Int xcor_subsample
	Int xcor_subsample_default = 25000000

	String pe_ta
	String se_ta

	String ref_pe_xcor_log
	String ref_pe_xcor_log_subsample
	String ref_se_xcor_log
	String ref_se_xcor_log_subsample
	String mito_chr_name = 'chrM'

	Int xcor_cpu = 1
	Int xcor_mem_mb = 16000
	Int xcor_time_hr = 6
	String xcor_disks = 'local-disk 100 HDD'

	call atac.xcor as pe_xcor { input :
		ta = pe_ta,
		subsample = xcor_subsample_default,
		paired_end = true,
		mito_chr_name = mito_chr_name,

		cpu = xcor_cpu,
		mem_mb = xcor_mem_mb,
		time_hr = xcor_time_hr,
		disks = xcor_disks,
	}
	call atac.xcor as pe_xcor_subsample { input :
		ta = pe_ta,
		subsample = xcor_subsample,
		paired_end = true,
		mito_chr_name = mito_chr_name,

		cpu = xcor_cpu,
		mem_mb = xcor_mem_mb,
		time_hr = xcor_time_hr,
		disks = xcor_disks,
	}
	call atac.xcor as se_xcor { input :
		ta = se_ta,
		subsample = xcor_subsample_default,
		paired_end = false,
		mito_chr_name = mito_chr_name,

		cpu = xcor_cpu,
		mem_mb = xcor_mem_mb,
		time_hr = xcor_time_hr,
		disks = xcor_disks,
	}
	call atac.xcor as se_xcor_subsample { input :
		ta = se_ta,
		subsample = xcor_subsample,
		paired_end = false,
		mito_chr_name = mito_chr_name,

		cpu = xcor_cpu,
		mem_mb = xcor_mem_mb,
		time_hr = xcor_time_hr,
		disks = xcor_disks,
	}

	call compare_md5sum.compare_md5sum { input :
		labels = [
			'pe_xcor',
			'pe_xcor_subsample',
			'se_xcor',
			'se_xcor_subsample',
		],
		files = [
			pe_xcor.score,
			pe_xcor_subsample.score,
			se_xcor.score,
			se_xcor_subsample.score,
		],
		ref_files = [
			ref_pe_xcor_log,
			ref_pe_xcor_log_subsample,
			ref_se_xcor_log,
			ref_se_xcor_log_subsample,
		],
	}
}
