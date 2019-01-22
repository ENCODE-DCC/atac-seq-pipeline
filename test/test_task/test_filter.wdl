# ENCODE DCC ATAC-Seq/DNase-Seq pipeline tester for task filter
# Author: Jin Lee (leepc12@gmail.com)
import "../../atac.wdl" as atac

workflow test_filter {
	Int multimapping

	String dup_marker = 'picard'
	Int mapq_thresh = 30

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
	String mito_chr_name = 'chrM'

	Int filter_cpu = 1
	Int filter_mem_mb = 20000
	Int filter_time_hr = 24
	String filter_disks = "local-disk 100 HDD"

	call atac.filter as pe_filter { input :
		bam = pe_bam,
		multimapping = multimapping,
		paired_end = true,
		mito_chr_name = mito_chr_name,

		dup_marker = dup_marker,
		mapq_thresh = mapq_thresh,
		no_dup_removal = false,

		cpu = filter_cpu,
		mem_mb = filter_mem_mb,
		time_hr = filter_time_hr,
		disks = filter_disks,
	}
	call atac.filter as pe_filter_no_multimapping { input :
		bam = pe_bam_no_multimapping,
		multimapping = 0,
		paired_end = true,
		mito_chr_name = mito_chr_name,

		dup_marker = dup_marker,
		mapq_thresh = mapq_thresh,
		no_dup_removal = false,

		cpu = filter_cpu,
		mem_mb = filter_mem_mb,
		time_hr = filter_time_hr,
		disks = filter_disks,
	}
	call atac.filter as pe_filter_no_dup_removal { input :
		bam = pe_bam,
		multimapping = multimapping,
		paired_end = true,
		mito_chr_name = mito_chr_name,

		dup_marker = dup_marker,
		mapq_thresh = mapq_thresh,
		no_dup_removal = true,

		cpu = filter_cpu,
		mem_mb = filter_mem_mb,
		time_hr = filter_time_hr,
		disks = filter_disks,
	}
	call atac.filter as se_filter { input :
		bam = se_bam,
		multimapping = multimapping,
		paired_end = false,
		mito_chr_name = mito_chr_name,

		dup_marker = dup_marker,
		mapq_thresh = mapq_thresh,
		no_dup_removal = false,

		cpu = filter_cpu,
		mem_mb = filter_mem_mb,
		time_hr = filter_time_hr,
		disks = filter_disks,
	}
	call atac.filter as se_filter_no_multimapping { input :
		bam = se_bam_no_multimapping,
		multimapping = 0,
		paired_end = false,
		mito_chr_name = mito_chr_name,

		dup_marker = dup_marker,
		mapq_thresh = mapq_thresh,
		no_dup_removal = false,

		cpu = filter_cpu,
		mem_mb = filter_mem_mb,
		time_hr = filter_time_hr,
		disks = filter_disks,
	}
	call atac.filter as se_filter_no_dup_removal { input :
		bam = se_bam,
		multimapping = multimapping,
		paired_end = false,
		mito_chr_name = mito_chr_name,

		dup_marker = dup_marker,
		mapq_thresh = mapq_thresh,
		no_dup_removal = true,

		cpu = filter_cpu,
		mem_mb = filter_mem_mb,
		time_hr = filter_time_hr,
		disks = filter_disks,
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
