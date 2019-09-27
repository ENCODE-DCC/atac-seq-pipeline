# ENCODE DCC ATAC-Seq/DNase-Seq pipeline tester
# Author: Jin Lee (leepc12@gmail.com)
import '../../../atac.wdl' as atac
import 'compare_md5sum.wdl' as compare_md5sum

workflow test_compare_signal_to_roadmap {
	File pval_bw
	File reg2map_bed
	File reg2map
	File roadmap_meta

	File ref_roadmap_compare_log

	call atac.compare_signal_to_roadmap { input : 
		pval_bw = pval_bw,

		reg2map_bed = reg2map_bed,
		reg2map = reg2map,
		roadmap_meta = roadmap_meta,
	}

	call compare_md5sum.compare_md5sum { input :
		labels = [
			'ref_roadmap_compare_log',
		],
		files = [
			compare_signal_to_roadmap.roadmap_compare_log,
		],
		ref_files = [
			ref_roadmap_compare_log,
		],
	}
}
