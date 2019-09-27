# ENCODE DCC ATAC-Seq/DNase-Seq pipeline tester
# Author: Jin Lee (leepc12@gmail.com)
import '../../../atac.wdl' as atac
import 'compare_md5sum.wdl' as compare_md5sum

workflow test_preseq {
	File bam
	Boolean paired_end

	#File ref_qc_html
	File ref_picard_est_lib_size_qc
	File ref_preseq_log

	call atac.preseq { input : 
		paired_end = paired_end,
		bam = bam,

		mem_mb = 4000,
	}

	call compare_md5sum.compare_md5sum { input :
		labels = [
			'test_picard_est_lib_size_qc',
			'test_preseq_log',
		],
		files = select_all([
			preseq.picard_est_lib_size_qc,
			preseq.preseq_log,
		]),
		ref_files = [
			ref_picard_est_lib_size_qc,
			ref_preseq_log,
		],
	}
}
