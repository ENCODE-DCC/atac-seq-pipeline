# ENCODE DCC ATAC-Seq/DNase-Seq pipeline tester
# Author: Jin Lee (leepc12@gmail.com)

# DO NOT REMOVE ANY COMMENTS STARTING WITH #@
# This WDL script is not stand-alone.
# You will need to append TASK_DEF block from ../atac.wdl

#@WORKFLOW_DEF_BEGIN
workflow test_atac {
	### pe samples (replicate 1)
	Array[Array[String]] pe_adapters
	Array[Array[String]] pe_fastqs
	String pe_trimmed_fastq_R1
	String pe_trimmed_fastq_R2
	String pe_bam
	String pe_bam_no_multimapping
	String pe_nodup_bam
	String pe_nodup_bam_no_multimapping
	String pe_ta

	### se samples (replicate 1), _rep2 is for replicate 2
	Array[Array[String]] se_adapters
	Array[Array[String]] se_fastqs
	String se_trimmed_fastq
	String se_bam
	String se_bam_no_multimapping
	String se_nodup_bam
	String se_nodup_bam_no_multimapping
	String se_ta

	String se_ta_rep2
	String se_ta_pooled
	String se_peak_rep1 # test overlap,idr for SE set only
	String se_peak_rep2
	String se_peak_pooled
	String se_overlap_peak_rep1_vs_rep2
	String se_overlap_peak_rep1_pr
	String se_overlap_peak_rep2_pr
	String se_overlap_peak_ppr

	### pe reference
	String ref_pe_trimmed_fastq_R1
	String ref_pe_trimmed_fastq_R2
	String ref_pe_flagstat
	String ref_pe_flagstat_no_multimapping
	String ref_pe_nodup_bam
	String ref_pe_nodup_bam_no_multimapping
	String ref_pe_filt_bam # with flag no_dup_removal on
	String ref_pe_ta
	String ref_pe_ta_disable_tn5_shift
	String ref_pe_ta_subsample
	String ref_pe_xcor_log
	String ref_pe_ta_pr1

	### se reference
	String ref_se_trimmed_fastq
	String ref_se_flagstat
	String ref_se_flagstat_no_multimapping
	String ref_se_nodup_bam
	String ref_se_nodup_bam_no_multimapping
	String ref_se_filt_bam # with flag no_dup_removal on
	String ref_se_ta
	String ref_se_ta_disable_tn5_shift
	String ref_se_ta_subsample
	String ref_se_xcor_log
	String ref_se_ta_pr1

	String ref_se_pooled_ta
	String ref_se_macs2_frip_qc # test macs2 for SE set only
	String ref_se_macs2_sig_pval
	String ref_se_overlap_frip_qc_rep1_vs_rep2
	String ref_se_idr_frip_qc_rep1_vs_rep2
	String ref_se_reproducibility_overlap_qc

	# genome data and parameters
	File pe_genome_tsv
	File se_genome_tsv

	# common parameters for both samples	
	Int multimapping
	Int cap_num_peak
	Float pval_thresh
	Int smooth_win
	Float idr_thresh
	Int bam2ta_subsample
	Int xcor_subsample

	# read genome data from tsv file for se and pe
	call read_genome_tsv as pe_read_genome_tsv { input: genome_tsv = pe_genome_tsv }
	call read_genome_tsv as se_read_genome_tsv { input: genome_tsv = se_genome_tsv }
	String pe_bowtie2_idx_tar = pe_read_genome_tsv.genome['bowtie2_idx_tar']
	String pe_blacklist = pe_read_genome_tsv.genome['blacklist']
	String pe_chrsz = pe_read_genome_tsv.genome['chrsz']
	String pe_gensz = pe_read_genome_tsv.genome['gensz']
	String se_bowtie2_idx_tar = se_read_genome_tsv.genome['bowtie2_idx_tar']
	String se_blacklist = se_read_genome_tsv.genome['blacklist']
	String se_chrsz = se_read_genome_tsv.genome['chrsz']
	String se_gensz = se_read_genome_tsv.genome['gensz']

	### test pe
	call trim_adapter as pe_trim_adapter { input :
		fastqs = pe_fastqs,
		adapters = pe_adapters,
		auto_detect_adapter = false,
		paired_end = true,
		cpu = 1,
	}
	call trim_adapter as pe_trim_adapter_auto { input :
		fastqs = pe_fastqs,
		adapters = [],
		auto_detect_adapter = true,
		paired_end = true,
		cpu = 1,
	}
	call bowtie2 as pe_bowtie2 { input :
		idx_tar = pe_bowtie2_idx_tar,
		fastqs = [pe_trimmed_fastq_R1, pe_trimmed_fastq_R2],
		multimapping = multimapping,
		paired_end = true,
		cpu = 2,
	}
	call bowtie2 as pe_bowtie2_no_multimapping { input :
		idx_tar = pe_bowtie2_idx_tar,
		fastqs = [pe_trimmed_fastq_R1, pe_trimmed_fastq_R2],
		paired_end = true,
		cpu = 2,
	}
	call filter as pe_filter { input :
		bam = pe_bam,
		multimapping = multimapping,
		paired_end = true,
		cpu = 1,
	}
	call filter as pe_filter_no_multimapping { input :
		bam = pe_bam_no_multimapping,
		paired_end = true,
		cpu = 1,
	}
	call filter as pe_filter_no_dup_removal { input :
		bam = pe_bam
		multimapping = multimapping,
		no_dup_removal = true,
		paired_end = true,
		cpu = 1,
	}
	call bam2ta as pe_bam2ta { input :
		bam = pe_nodup_bam
		disable_tn5_shift = false,
		paired_end = true,
	}
	call bam2ta as pe_bam2ta_disable_tn5_shift { input :
		bam = pe_nodup_bam
		disable_tn5_shift = true,
		paired_end = true,
	}
	call bam2ta as pe_bam2ta_subsample { input :
		bam = pe_nodup_bam
		disable_tn5_shift = false,
		subsample = bam2ta_subsample,
		paired_end = true,
	}
	call xcor as pe_xcor { input :
		ta = pe_ta,
		subsample = xcor_subsample,
		paired_end = true,
	}
	call spr as pe_spr { input :
		ta = pe_ta,
		paired_end = true,
	}	
	### test se
	call trim_adapter as se_trim_adapter { input :
		fastqs = se_fastqs,
		adapters = se_adapters,
		auto_detect_adapter = false,
		paired_end = false,
		cpu = 1,
	}
	call trim_adapter as se_trim_adapter_auto { input :
		fastqs = se_fastqs,
		adapters = [],
		auto_detect_adapter = true,
		paired_end = false,
		cpu = 1,
	}
	call bowtie2 as se_bowtie2 { input :
		idx_tar = se_bowtie2_idx_tar,
		fastqs = [se_trimmed_fastq],
		multimapping = multimapping,
		paired_end = false,
		cpu = 2,
	}
	call bowtie2 as se_bowtie2_no_multimapping { input :
		idx_tar = se_bowtie2_idx_tar,
		fastqs = [se_trimmed_fastq],
		paired_end = false,
		cpu = 2,
	}
	call filter as se_filter { input :
		bam = se_bam,
		paired_end = false,
		cpu = 1,
	}
	call filter as se_filter_no_multimapping { input :
		bam = se_bam_no_multimapping,
		paired_end = false,
		cpu = 1,
	}
	call filter as se_filter_no_dup_removal { input :
		bam = se_bam
		multimapping = multimapping,
		no_dup_removal = true,
		paired_end = false,
		cpu = 1,
	}
	call bam2ta as se_bam2ta { input :
		bam = se_nodup_bam
		disable_tn5_shift = false,
		paired_end = false,
	}
	call bam2ta as se_bam2ta_disable_tn5_shift { input :
		bam = se_nodup_bam
		disable_tn5_shift = true,
		paired_end = false,
	}
	call bam2ta as se_bam2ta_subsample { input :
		bam = se_nodup_bam
		disable_tn5_shift = false,
		subsample = bam2ta_subsample,
		paired_end = false,
	}
	call xcor as se_xcor { input :
		ta = se_ta,
		subsample = xcor_subsample,
		paired_end = false,
	}
	call spr as se_spr { input :
		ta = se_ta,
		paired_end = false,
	}
	call pool_ta as se_pool_ta { input :
		tas = [se_ta, se_ta_rep2],
	}

	# test peak-caling
	call macs2 as se_macs2 { input :
		ta = se_ta,
		gensz = se_gensz,
		chrsz = se_chrsz,
		cap_num_peak = cap_num_peak,
		pval_thresh = pval_thresh,
		smooth_win = smooth_win,
		make_signal = true,
		blacklist = se_blacklist,
	}
	call overlap as se_overlap { input :
		prefix = "rep1-rep2",
		peak1 = se_peak_rep1,
		peak2 = se_peak_rep2,
		peak_pooled = se_peak_pooled,
		peak_type = 'narrowPeak',
		blacklist = se_blacklist,
		ta = se_ta_pooled,
	}
	call idr as se_idr { input : 
		prefix = "rep1-rep2",
		peak1 = se_peak_rep1,
		peak2 = se_peak_rep2,
		peak_pooled = se_peak_pooled,
		idr_thresh = idr_thresh,
		peak_type = 'narrowPeak',
		rank = 'p.value',
		blacklist = se_blacklist,
		ta = se_ta_pooled,
	}
	call reproducibility as se_reproducibility_overlap { input :
		prefix = 'overlap',
		peaks = [se_overlap_peak_rep1_vs_rep2],
		peaks_pr = [se_overlap_peak_rep1_pr, se_overlap_peak_rep2_pr],
		peak_ppr = se_overlap_peak_ppr,
	}

	call compare_md5sum { input :
		labels = [
			'pe_trim_adapter_R1',
			'pe_trim_adapter_R2',
			'pe_trim_adapter_auto_R1',
			'pe_trim_adapter_auto_R2',
			'pe_bowtie2',
			'pe_bowtie2_no_multimapping',
			'pe_filter',
			'pe_filter_no_multimapping',
			'pe_filter_no_dup_removal',
			'pe_bam2ta',
			'pe_bam2ta_disable_tn5_shift',
			'pe_bam2ta_subsample',
			'pe_xcor',
			'pe_spr',

			'se_trim_adapter',
			'se_trim_adapter_auto',
			'se_bowtie2',
			'se_bowtie2_no_multimapping',
			'se_filter',
			'se_filter_no_multimapping',
			'se_filter_no_dup_removal',
			'se_bam2ta',
			'se_bam2ta_disable_tn5_shift',
			'se_bam2ta_subsample',
			'se_xcor',
			'se_spr',

			'se_pool_ta',

			'se_macs2',
			'se_macs2_sig_pval',
			'se_overlap',
			'se_idr',
			'se_reproducibility_overlap',
		],

		files = [
			pe_trim_adapter.trimmed_merged_fastqs[0], 
			pe_trim_adapter.trimmed_merged_fastqs[1],
			pe_trim_adapter_auto.trimmed_merged_fastqs[0], 
			pe_trim_adapter_auto.trimmed_merged_fastqs[1],
			pe_bowtie2.flagstat_qc,
			pe_bowtie2_no_multimapping.flagstat_qc,
			pe_filter.nodup_bam,
			pe_filter_no_multimapping.nodup_bam,
			pe_filter_no_dup_removal.nodup_bam,
			pe_bam2ta.ta,
			pe_bam2ta_disable_tn5_shift.ta,
			pe_bam2ta_subsample.ta,
			pe_xcor.score,
			pe_spr.ta_pr1,

			se_trim_adapter.trimmed_merged_fastqs[0], 
			se_trim_adapter_auto.trimmed_merged_fastqs[0], 
			se_bowtie2.flagstat_qc,
			se_bowtie2_no_multimapping.flagstat_qc,
			se_filter.nodup_bam,
			se_filter_no_multimapping.nodup_bam,
			se_filter_no_dup_removal.nodup_bam,
			se_bam2ta.ta,
			se_bam2ta_disable_tn5_shift.ta,
			se_bam2ta_subsample.ta,
			se_xcor.score,			
			se_spr.ta_pr1,

			se_pool_ta.ta_pooled,

			se_macs2.frip_qc,
			se_macs2.sig_pval,
			se_overlap.frip_qc,
			se_idr.frip_qc,
			se_reproducibility_overlap.reproducibility_qc,
		],

		ref_files = [
			ref_pe_trimmed_fastq_R1,
			ref_pe_trimmed_fastq_R2,
			ref_pe_trimmed_fastq_R1,
			ref_pe_trimmed_fastq_R2,
			ref_pe_flagstat,
			ref_pe_flagstat_no_multimapping,
			ref_pe_nodup_bam,
			ref_pe_nodup_bam_no_multimapping,
			ref_pe_filt_bam,
			ref_pe_ta,
			ref_pe_ta_disable_tn5_shift,
			ref_pe_ta_subsample,
			ref_pe_xcor_log,
			ref_pe_ta_pr1,

			ref_se_trimmed_fastq,
			ref_se_trimmed_fastq,
			ref_se_flagstat,
			ref_se_flagstat_no_multimapping,
			ref_se_nodup_bam,
			ref_se_nodup_bam_no_multimapping,
			ref_se_filt_bam,
			ref_se_ta,
			ref_se_ta_disable_tn5_shift,
			ref_se_ta_subsample,
			ref_se_xcor_log,
			ref_se_ta_pr1,

			ref_se_pooled_ta,

			ref_se_macs2_frip_qc,
			ref_se_macs2_sig_pval,
			ref_se_overlap_frip_qc_rep1_vs_rep2,
			ref_se_idr_frip_qc_rep1_vs_rep2,
			ref_se_reproducibility_overlap_qc,
		],
	}
}
#@WORKFLOW_DEF_END

#@TASK_DEF_BEGIN
task compare_md5sum {
	Array[String] labels
	Array[File] files
	Array[File] ref_files

	command <<<
		python <<CODE	
		from collections import OrderedDict
		from json import dumps
		import os
		import json
		import hashlib

		def md5sum(filename, blocksize=65536):
		    hash = hashlib.md5()
		    with open(filename, 'rb') as f:
		        for block in iter(lambda: f.read(blocksize), b""):
		            hash.update(block)
		    return hash.hexdigest()

		with open('${write_lines(labels)}','r') as fp:
			labels = fp.read().splitlines()
		with open('${write_lines(files)}','r') as fp:
			files = fp.read().splitlines()
		with open('${write_lines(ref_files)}','r') as fp:
			ref_files = fp.read().splitlines()

		result = []
		match = OrderedDict()
		match_overall = True
		for i, label in enumerate(labels):
			f = files[i]
			ref_f = ref_files[i]
			md5 = md5sum(f)
			ref_md5 = md5sum(ref_f)
			# if text file, read in contents
			if f.endswith('.qc') or f.endswith('.txt') or \
				f.endswith('.log') or f.endswith('.out'):
				with open(f,'r') as fp:
					contents = fp.read()
				with open(ref_f,'r') as fp:
					ref_contents = fp.read()
			else:
				contents = ''
				ref_contents = ''
			matched = md5==ref_md5
			result.append(OrderedDict([
				('label', label),
				('match', matched),
				('md5sum', md5),
				('ref_md5sum', ref_md5),
				('basename', os.path.basename(f)),
				('ref_basename', os.path.basename(ref_f)),
				('contents', contents),
				('ref_contents', ref_contents),
				]))
			match[label] = matched
			match_overall &= matched
		with open('result.json','w') as fp:
			fp.write(json.dumps(result, indent=4))
		match_tmp = []
		for key in match:
			val = match[key]
			match_tmp.append('{}\t{}'.format(key, val))
		with open('match.tsv','w') as fp:
			fp.writelines('\n'.join(match_tmp))
		with open('match_overall.txt','w') as fp:
			fp.write(str(match_overall))
		CODE
	>>>
	output {
		Map[String,String] match = read_map('match.tsv') # key:label, val:match
		Boolean match_overall = read_boolean('match_overall.txt')
		File json = glob('result.json')[0] # details (json file)
		String json_str =read_string('result.json') # details (string)
	}
}
#@TASK_DEF_END