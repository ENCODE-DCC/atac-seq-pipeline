# ENCODE DCC ATAC-Seq/DNase-Seq pipeline tester
# Author: Jin Lee (leepc12@gmail.com)
import '../../../atac.wdl' as atac
import 'compare_md5sum.wdl' as compare_md5sum

workflow test_fraglen_stat_pe {
	File nodup_bam

	File ref_nucleosomal_qc

	call atac.fraglen_stat_pe { input : 
		nodup_bam = nodup_bam,
	}

	call compare_md5sum.compare_md5sum { input :
		labels = [
			'test_nucleosomal_qc',
		],
		files = [
			fraglen_stat_pe.nucleosomal_qc,
		],
		ref_files = [
			ref_nucleosomal_qc,
		],
	}
}
