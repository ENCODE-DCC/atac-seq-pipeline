# ENCODE DCC ATAC-Seq/DNase-Seq pipeline tester for task bowtie2
# Author: Jin Lee (leepc12@gmail.com)
import "../../atac.wdl" as atac

workflow test_bowtie2 {
	Int multimapping

	Array[String] pe_trimmed_fastqs
	Array[String] se_trimmed_fastqs

	# we don't compare BAM because BAM's header includes date
	# hence md5sums don't match all the time
	String ref_pe_flagstat
	String ref_pe_flagstat_no_multimapping

	String ref_se_flagstat
	String ref_se_flagstat_no_multimapping

	String pe_bowtie2_idx_tar
	String se_bowtie2_idx_tar

	call atac.bowtie2 as pe_bowtie2 { input :
		idx_tar = pe_bowtie2_idx_tar,
		fastqs = pe_trimmed_fastqs,
		multimapping = multimapping,
		paired_end = true,
		cpu = 1,
	}
	call atac.bowtie2 as pe_bowtie2_no_multimapping { input :
		idx_tar = pe_bowtie2_idx_tar,
		fastqs = pe_trimmed_fastqs,
		paired_end = true,
		cpu = 1,
	}
	call atac.bowtie2 as se_bowtie2 { input :
		idx_tar = se_bowtie2_idx_tar,
		fastqs = se_trimmed_fastqs,
		multimapping = multimapping,
		paired_end = false,
		cpu = 1,
	}
	call atac.bowtie2 as se_bowtie2_no_multimapping { input :
		idx_tar = se_bowtie2_idx_tar,
		fastqs = se_trimmed_fastqs,
		paired_end = false,
		cpu = 1,
	}

	call atac.compare_md5sum { input :
		labels = [
			'pe_bowtie2',
			'pe_bowtie2_no_multimapping',
			'se_bowtie2',
			'se_bowtie2_no_multimapping',
		],
		files = [
			pe_bowtie2.flagstat_qc,
			pe_bowtie2_no_multimapping.flagstat_qc,
			se_bowtie2.flagstat_qc,
			se_bowtie2_no_multimapping.flagstat_qc,
		],
		ref_files = [
			ref_pe_flagstat,
			ref_pe_flagstat_no_multimapping,
			ref_se_flagstat,
			ref_se_flagstat_no_multimapping,
		],
	}
}
