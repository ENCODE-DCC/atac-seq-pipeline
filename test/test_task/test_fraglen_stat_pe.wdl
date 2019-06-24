# ENCODE DCC ATAC-Seq/DNase-Seq pipeline tester
# Author: Jin Lee (leepc12@gmail.com)
import "../../atac.wdl" as atac
import "compare_md5sum.wdl" as compare_md5sum

workflow test_fraglen_stat_pe {
	File nodup_bam

	File ref_nucleosomal_qc
	File ref_fraglen_dist_plot

	call atac.fraglen_stat_pe { input : 
		nodup_bam = nodup_bam,
	}

	call compare_md5sum.compare_md5sum { input :
		labels = [
			'test_nucleosomal_qc',
			'test_fraglen_dist_plot',
		],
		files = [
			fraglen_stat_pe.nucleosomal_qc,
			fraglen_stat_pe.fraglen_dist_plot,
		],
		ref_files = [
			ref_nucleosomal_qc,
			ref_fraglen_dist_plot,
		],
	}
}
