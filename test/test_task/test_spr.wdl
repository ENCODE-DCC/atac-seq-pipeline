# ENCODE DCC ATAC-Seq/DNase-Seq pipeline tester
# Author: Jin Lee (leepc12@gmail.com)
import "../../atac.wdl" as atac

workflow test_spr {
	String pe_ta
	String se_ta

	String ref_pe_ta_pr1
	String ref_pe_ta_pr2
	String ref_se_ta_pr1
	String ref_se_ta_pr2

	call atac.spr as pe_spr { input :
		ta = pe_ta,
		paired_end = true,
	}	
	call atac.spr as se_spr { input :
		ta = se_ta,
		paired_end = false,
	}

	call atac.compare_md5sum { input :
		labels = [
			'pe_spr_pr1',
			'pe_spr_pr2',
			'se_spr_pr1',
			'se_spr_pr2',
		],
		files = [
			pe_spr.ta_pr1,
			pe_spr.ta_pr2,
			se_spr.ta_pr1,
			se_spr.ta_pr2,
		],
		ref_files = [
			ref_pe_ta_pr1,
			ref_pe_ta_pr2,
			ref_se_ta_pr1,
			ref_se_ta_pr2,
		],
	}
}
