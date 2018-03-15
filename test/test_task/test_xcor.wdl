# ENCODE DCC ATAC-Seq/DNase-Seq pipeline tester
# Author: Jin Lee (leepc12@gmail.com)
import "../../atac.wdl" as atac

workflow test_xcor {
	Int xcor_subsample

	String pe_ta
	String se_ta

	String ref_pe_xcor_log
	String ref_pe_xcor_log_subsample
	String ref_se_xcor_log
	String ref_se_xcor_log_subsample

	call atac.xcor as pe_xcor { input :
		ta = pe_ta,
		paired_end = true,
	}
	call atac.xcor as pe_xcor_subsample { input :
		ta = pe_ta,
		subsample = xcor_subsample,
		paired_end = true,
	}
	call atac.xcor as se_xcor { input :
		ta = se_ta,
		paired_end = false,
	}
	call atac.xcor as se_xcor_subsample { input :
		ta = se_ta,
		subsample = xcor_subsample,
		paired_end = false,
	}

	call atac.compare_md5sum { input :
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
