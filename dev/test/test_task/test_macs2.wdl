# ENCODE DCC ATAC-Seq/DNase-Seq pipeline tester for task macs2
# Author: Jin Lee (leepc12@gmail.com)
import '../../../atac.wdl' as atac
import 'compare_md5sum.wdl' as compare_md5sum

workflow test_macs2 {
	Int cap_num_peak
	Float pval_thresh
	Int smooth_win

	# test macs2 for SE set only
	String se_ta

	String ref_se_macs2_npeak # raw narrow-peak
	String ref_se_macs2_bfilt_npeak # blacklist filtered narrow-peak
	String ref_se_macs2_frip_qc 

	String se_blacklist
	String se_chrsz
	String se_gensz

	Int macs2_mem_mb = 16000
	Int macs2_time_hr = 24
	String macs2_disks = 'local-disk 100 HDD'

	call atac.call_peak as se_macs2 { input :
		peak_caller = 'macs2',
		peak_type = 'narrowPeak',
		ta = se_ta,
		gensz = se_gensz,
		chrsz = se_chrsz,
		cap_num_peak = cap_num_peak,
		pval_thresh = pval_thresh,
		smooth_win = smooth_win,
		blacklist = se_blacklist,
		keep_irregular_chr_in_bfilt_peak = false,

		cpu = 1,
		mem_mb = macs2_mem_mb,
		time_hr = macs2_time_hr,
		disks = macs2_disks,
	}

	call compare_md5sum.compare_md5sum { input :
		labels = [
			'se_macs2_npeak',
			'se_macs2_bfilt_npeak',
			'se_macs2_frip_qc',
		],
		files = [
			se_macs2.peak,
			se_macs2.bfilt_peak,
			se_macs2.frip_qc,
		],
		ref_files = [
			ref_se_macs2_npeak,
			ref_se_macs2_bfilt_npeak,
			ref_se_macs2_frip_qc,
		],
	}
}
