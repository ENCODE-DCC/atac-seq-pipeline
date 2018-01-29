# ENCODE DCC ATAC-Seq/DNase-Seq pipeline tester for task filter
# Author: Jin Lee (leepc12@gmail.com)
import "../../atac.wdl" as atac

workflow test_filter {
	Int multimapping

	String pe_bam
	String pe_bam_no_multimapping
	String se_bam
	String se_bam_no_multimapping

	String ref_pe_nodup_bam
	String ref_pe_nodup_bam_no_multimapping
	String ref_pe_filt_bam
	String ref_se_nodup_bam
	String ref_se_nodup_bam_no_multimapping
	String ref_se_filt_bam

	call atac.filter as pe_filter { input :
		bam = pe_bam,
		multimapping = multimapping,
		paired_end = true,
		cpu = 1,
	}
	call atac.filter as pe_filter_no_multimapping { input :
		bam = pe_bam_no_multimapping,
		paired_end = true,
		cpu = 1,
	}
	call atac.filter as pe_filter_no_dup_removal { input :
		bam = pe_bam,
		multimapping = multimapping,
		no_dup_removal = true,
		paired_end = true,
		cpu = 1,
	}
	call atac.filter as se_filter { input :
		bam = se_bam,
		multimapping = multimapping,
		paired_end = false,
		cpu = 1,
	}
	call atac.filter as se_filter_no_multimapping { input :
		bam = se_bam_no_multimapping,
		paired_end = false,
		cpu = 1,
	}
	call atac.filter as se_filter_no_dup_removal { input :
		bam = se_bam,
		multimapping = multimapping,
		no_dup_removal = true,
		paired_end = false,
		cpu = 1,
	}

	call atac.compare_md5sum { input :
		labels = [
			'pe_filter',
			'pe_filter_no_multimapping',
			'pe_filter_no_dup_removal',
			'se_filter',
			'se_filter_no_multimapping',
			'se_filter_no_dup_removal',
		],
		files = [
			pe_filter.nodup_bam,
			pe_filter_no_multimapping.nodup_bam,
			pe_filter_no_dup_removal.nodup_bam,
			se_filter.nodup_bam,
			se_filter_no_multimapping.nodup_bam,
			se_filter_no_dup_removal.nodup_bam,
		],
		ref_files = [
			ref_pe_nodup_bam,
			ref_pe_nodup_bam_no_multimapping,
			ref_pe_filt_bam,
			ref_se_nodup_bam,
			ref_se_nodup_bam_no_multimapping,
			ref_se_filt_bam,
		],
	}
}
