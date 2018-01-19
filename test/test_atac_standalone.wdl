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
		File json = glob('result.json')[0] # details
	}
}
#@TASK_DEF_END#@TASK_DEF_BEGIN
### genomic tasks
task trim_adapter { # trim adapters and merge trimmed fastqs
	# parameters from workflow
	Array[Array[File]] fastqs 		# [merge_id][end_id]
	Array[Array[String]] adapters 	# [merge_id][end_id]
	Boolean paired_end
	# mandatory
	Boolean auto_detect_adapter		# automatically detect/trim adapters
	# optional
	Int? min_trim_len 		# minimum trim length for cutadapt -m
	Float? err_rate			# Maximum allowed adapter error rate 
							# for cutadapt -e	
	# resource
	Int? cpu
	Int? mem_mb
	Int? time_hr
	String? disks

	command {
		python $(which encode_trim_adapter.py) \
			${write_tsv(fastqs)} \
			--adapters ${write_tsv(adapters)} \
			${if paired_end then "--paired-end" else ""} \
			${if auto_detect_adapter then "--auto-detect-adapter" else ""} \
			${"--min-trim-len " + min_trim_len} \
			${"--err-rate " + err_rate} \
			${"--nth " + select_first([cpu,2])}
	}
	output {
		# WDL glob() globs in an alphabetical order
		# so R1 and R2 can be switched, which results in an
		# unexpected behavior of a workflow
		# so we prepend merge_fastqs_'end'_ (R1 or R2)
		# to the basename of original filename
		# this prefix will be later stripped in bowtie2 task
		Array[File] trimmed_merged_fastqs = glob("merge_fastqs_R?_*.fastq.gz")
	}
	runtime {
		cpu : select_first([cpu,2])
		memory : "${select_first([mem_mb,'10000'])} MB"
		time : select_first([time_hr,24])
		disks : select_first([disks,"local-disk 100 HDD"])
	}
}

task bowtie2 {
	# parameters from workflow
	File idx_tar 		# reference bowtie2 index tar
	Array[File] fastqs 	# [end_id]
	Boolean paired_end
	Int? multimapping
	# optional
	String? score_min 	# min acceptable alignment score func
						# w.r.t read length
	# resource
	Int? cpu
	Int? mem_mb
	Int? time_hr
	String? disks

	command {
		python $(which encode_bowtie2.py) \
			${idx_tar} \
			${sep=' ' fastqs} \
			${if paired_end then "--paired-end" else ""} \
			${"--multimapping " + multimapping} \
			${"--score-min " + score_min} \
			${"--nth " + select_first([cpu,4])}
	}
	output {
		File bam = glob("*.bam")[0]
		File bai = glob("*.bai")[0]
		File align_log = glob("*.align.log")[0]
		File flagstat_qc = glob("*.flagstat.qc")[0]
	}
	runtime {
		cpu : select_first([cpu,4])
		memory : "${select_first([mem_mb,'20000'])} MB"
		time : select_first([time_hr,48])
		disks : select_first([disks,"local-disk 100 HDD"])
		preemptible: 0
	}
}

task xcor {
	# parameters from workflow
	File ta
	Boolean paired_end
	# optional
	Int? subsample 		# number of reads to subsample TAGALIGN
						# this will be used for xcor only
						# will not affect any downstream analysis
	# resource
	Int? cpu
	Int? mem_mb	
	Int? time_hr
	String? disks

	command {
		python $(which encode_xcor.py) \
			${ta} \
			${if paired_end then "--paired-end" else ""} \
			${"--subsample " + select_first([subsample,25000000])} \
			--speak=0 \
			${"--nth " + select_first([cpu,2])}
	}
	output {
		File plot_pdf = glob("*.cc.plot.pdf")[0]
		File plot_png = glob("*.cc.plot.png")[0]
		File score = glob("*.cc.qc")[0]
		Int fraglen = read_int(glob("*.cc.fraglen.txt")[0])
	}
	runtime {
		cpu : select_first([cpu,2])
		memory : "${select_first([mem_mb,'10000'])} MB"
		time : select_first([time_hr,6])
		disks : select_first([disks,"local-disk 100 HDD"])
	}
}

task macs2 {
	# parameters from workflow
	File ta
	String gensz		# Genome size (sum of entries in 2nd column of 
                        # chr. sizes file, or hs for human, ms for mouse)
	File chrsz			# 2-col chromosome sizes file
	Int? cap_num_peak	# cap number of raw peaks called from MACS2
	Float? pval_thresh	# p.value threshold
	Int? smooth_win		# size of smoothing window
	Boolean? make_signal
	File? blacklist 	# blacklist BED to filter raw peaks
	# fixed var
	String peak_type = "narrowPeak"
	# resource
	Int? mem_mb
	Int? time_hr
	String? disks

	command {
		python $(which encode_macs2_atac.py) \
			${ta} \
			${"--gensz "+ gensz} \
			${"--chrsz " + chrsz} \
			${"--cap-num-peak " + select_first([cap_num_peak,300000])} \
			${"--pval-thresh "+ pval_thresh} \
			${"--smooth-win "+ smooth_win} \
			${if select_first([make_signal,false]) then "--make-signal" else ""} \
			${"--blacklist "+ blacklist}
		
		# ugly part to deal with optional outputs with Google JES backend
		${if select_first([make_signal,false]) then "" 
			else "touch null.pval.signal.bigwig null.fc.signal.bigwig"}
		${if defined(blacklist) then "" 
			else "touch null.bfilt."+peak_type+".gz"}
		touch null 
	}
	output {
		File npeak = glob("*[!.][!b][!f][!i][!l][!t]."+peak_type+".gz")[0]
		File bfilt_npeak = if defined(blacklist) then glob("*.bfilt."+peak_type+".gz")[0] else npeak
		File sig_pval = if select_first([make_signal,false]) then glob("*.pval.signal.bigwig")[0] else glob("null")[0]
		File sig_fc = if select_first([make_signal,false]) then glob("*.fc.signal.bigwig")[0] else glob("null")[0]
		File frip_qc = glob("*.frip.qc")[0]
	}
	runtime {
		memory : "${select_first([mem_mb,'16000'])} MB"
		time : select_first([time_hr,24])
		disks : select_first([disks,"local-disk 100 HDD"])
	}
}

task filter {
	# parameters from workflow
	File bam
	Boolean paired_end
	Int? multimapping
	# optional
	String? dup_marker 			# picard.jar MarkDuplicates (picard) or 
								# sambamba markdup (sambamba)
	Int? mapq_thresh			# threshold for low MAPQ reads removal
	Boolean? no_dup_removal 	# no dupe reads removal when filtering BAM
								# dup.qc and pbc.qc will be empty files
								# and nodup_bam in the output is 
	# resource					# filtered bam with dupes	
	Int? cpu
	Int? mem_mb
	Int? time_hr
	String? disks

	command {
		python $(which encode_filter.py) \
			${bam} \
			${if paired_end then "--paired-end" else ""} \
			${"--multimapping " + multimapping} \
			${"--dup-marker " + dup_marker} \
			${"--mapq-thresh " + mapq_thresh} \
			${if select_first([no_dup_removal,false]) then "--no-dup-removal" else ""} \
			${"--nth " + cpu}
		# ugly part to deal with optional outputs with Google JES backend
		${if select_first([no_dup_removal,false]) then 
			"touch null.dup.qc null.pbc.qc; " else ""}
		touch null
	}
	output {
		File nodup_bam = glob("*.bam")[0]
		File nodup_bai = glob("*.bai")[0]
		File flagstat_qc = glob("*.flagstat.qc")[0]
		File dup_qc = if select_first([no_dup_removal,false]) then glob("null")[0] else glob("*.dup.qc")[0]
		File pbc_qc = if select_first([no_dup_removal,false]) then glob("null")[0] else glob("*.pbc.qc")[0]
	}
	runtime {
		cpu : select_first([cpu,2])
		memory : "${select_first([mem_mb,'20000'])} MB"
		time : select_first([time_hr,24])
		disks : select_first([disks,"local-disk 100 HDD"])
	}
}

task bam2ta {
	# parameters from workflow
	File bam
	Boolean paired_end
	Boolean disable_tn5_shift 	# no tn5 shifting (it's for dnase-seq)
	# optional
	String? regex_grep_v_ta 	# Perl-style regular expression pattern 
                        		# to remove matching reads from TAGALIGN
	Int? subsample 				# number of reads to subsample TAGALIGN
								# this affects all downstream analysis
	# resource
	Int? cpu
	Int? mem_mb
	Int? time_hr
	String? disks

	command {
		python $(which encode_bam2ta.py) \
			${bam} \
			${if paired_end then "--paired-end" else ""} \
			${if disable_tn5_shift then "--disable-tn5-shift" else ""} \
			${"--regex-grep-v-ta " +"'"+regex_grep_v_ta+"'"} \
			${"--subsample " + subsample} \
			${"--nth " + cpu}
	}
	output {
		File ta = glob("*.tagAlign.gz")[0]
	}
	runtime {
		cpu : select_first([cpu,2])
		memory : "${select_first([mem_mb,'10000'])} MB"
		time : select_first([time_hr,6])
		disks : select_first([disks,"local-disk 100 HDD"])
	}
}

task spr { # make two self pseudo replicates
	# parameters from workflow
	File ta
	Boolean paired_end

	# resource
	Int? mem_mb

	command {
		python $(which encode_spr.py) \
			${ta} \
			${if paired_end then "--paired-end" else ""}
	}
	output {
		File ta_pr1 = glob("*.pr1.tagAlign.gz")[0]
		File ta_pr2 = glob("*.pr2.tagAlign.gz")[0]
	}
	runtime {
		memory : "${select_first([mem_mb,'12000'])} MB"
	}
}

task pool_ta {
	# parameters from workflow
	Array[File] tas

	command {
		python $(which encode_pool_ta.py) \
			${sep=' ' tas}
	}
	output {
		File ta_pooled = glob("*.tagAlign.gz")[0]
	}
}

task idr {
	# parameters from workflow
	String? prefix 		# prefix for IDR output file
	File? peak1 			
	File? peak2
	File? peak_pooled
	Float? idr_thresh
	File? blacklist 	# blacklist BED to filter raw peaks
	# parameters to compute FRiP
	File? ta			# to calculate FRiP
	Int? fraglen 		# fragment length from xcor
	File? chrsz			# 2-col chromosome sizes file
	String peak_type
	String rank

	command {
		python $(which encode_idr.py) \
			${peak1} ${peak2} ${peak_pooled} \
			${"--prefix " + prefix} \
			${"--idr-thresh " + idr_thresh} \
			${"--peak-type " + peak_type} \
			--idr-rank ${rank} \
			${"--fraglen " + fraglen} \
			${"--chrsz " + chrsz} \
			${"--blacklist "+ blacklist} \
			${"--ta " + ta}

		# ugly part to deal with optional outputs with Google backend
		${if defined(blacklist) then "" 
			else "touch null.bfilt."+peak_type+".gz"}
		${if defined(ta) then "" 
			else "touch null.frip.qc"}			
		touch null 
	}
	output {
		File idr_peak = glob("*[!.][!b][!f][!i][!l][!t]."+peak_type+".gz")[0]
		File bfilt_idr_peak = if defined(blacklist) then 
							glob("*.bfilt."+peak_type+".gz")[0] else idr_peak
		File idr_plot = glob("*.txt.png")[0]
		File idr_unthresholded_peak = glob("*.txt.gz")[0]
		File idr_log = glob("*.log")[0]
		File frip_qc = if defined(ta) then glob("*.frip.qc")[0] else glob("null")[0]
	}
}

task overlap {
	# parameters from workflow
	String prefix 		# prefix for IDR output file
	File? peak1
	File? peak2
	File? peak_pooled
	File? blacklist 	# blacklist BED to filter raw peaks
	# parameters to compute FRiP
	File? ta			# to calculate FRiP
	Int? fraglen 		# fragment length from xcor
	File? chrsz			# 2-col chromosome sizes file
	String peak_type

	command {
		python $(which encode_naive_overlap.py) \
			${peak1} ${peak2} ${peak_pooled} \
			${"--prefix " + prefix} \
			${"--peak-type " + peak_type} \
			${"--fraglen " + fraglen} \
			${"--chrsz " + chrsz} \
			${"--blacklist "+ blacklist} \
			${"--ta " + ta}

		# ugly part to deal with optional outputs with Google backend
		${if defined(blacklist) then "" 
			else "touch null.bfilt."+peak_type+".gz"}
		${if defined(ta) then "" 
			else "touch null.frip.qc"}			
		touch null 
	}
	output {
		File overlap_peak = glob("*[!.][!b][!f][!i][!l][!t]."+peak_type+".gz")[0]
		File bfilt_overlap_peak = if defined(blacklist) then 
							glob("*.bfilt."+peak_type+".gz")[0] else overlap_peak
		File frip_qc = if defined(ta) then glob("*.frip.qc")[0] else glob("null")[0]
	}
}

task reproducibility {
	# parameters from workflow
	String prefix
	Array[File] peaks # peak files from pair of true replicates
						# in a sorted order. for example of 4 replicates,
						# 1,2 1,3 1,4 2,3 2,4 3,4.
                        # x,y means peak file from rep-x vs rep-y
	Array[File]? peaks_pr	# peak files from pseudo replicates
	File? peak_ppr			# Peak file from pooled pseudo replicate.

	command {
		python $(which encode_reproducibility_qc.py) \
			${sep=' ' peaks} \
			--peaks-pr ${sep=' ' peaks_pr} \
			${"--peak-ppr "+ peak_ppr} \
			--prefix ${prefix}
	}
	output {
		File reproducibility_qc = glob("*reproducibility.qc")[0]
	}
}

# gather all outputs and generate 
# - qc.html		: organized final HTML report
# - qc.json		: all QCs
task qc_report {
	# optional metadata
 	String? name # name of sample
	String? desc # description for sample
	#String? encode_accession_id	# ENCODE accession ID of sample
	# workflow params
	Boolean paired_end
	String pipeline_type
	String peak_caller
	Float idr_thresh
	# QCs
	Array[File?] flagstat_qcs
	Array[File?] nodup_flagstat_qcs
	Array[File?] dup_qcs
	Array[File?] pbc_qcs
	Array[File?] xcor_plots
	Array[File?] xcor_scores
	Array[File?] idr_plots
	Array[File?] idr_plots_pr
	File? idr_plot_ppr
	Array[File?] frip_qcs
	Array[File?] frip_qcs_pr1
	Array[File?] frip_qcs_pr2
	File? frip_qc_pooled
	File? frip_qc_ppr1 
	File? frip_qc_ppr2 
	Array[File?] frip_idr_qcs
	Array[File?] frip_idr_qcs_pr
	File? frip_idr_qc_ppr 
	Array[File?] frip_overlap_qcs
	Array[File?] frip_overlap_qcs_pr
	File? frip_overlap_qc_ppr
	File? idr_reproducibility_qc
	File? overlap_reproducibility_qc

	command {
		python $(which encode_qc_report.py) \
			${"--name '" + name + "'"} \
			${"--desc '" + desc + "'"} \
			${if paired_end then "--paired-end" else ""} \
			--pipeline-type ${pipeline_type} \
			--peak-caller ${peak_caller} \
			--idr-thresh ${idr_thresh} \
			--flagstat-qcs ${sep=' ' flagstat_qcs} \
			--nodup-flagstat-qcs ${sep=' ' nodup_flagstat_qcs} \
			--dup-qcs ${sep=' ' dup_qcs} \
			--pbc-qcs ${sep=' ' pbc_qcs} \
			--xcor-plots ${sep=' ' xcor_plots} \
			--xcor-scores ${sep=' ' xcor_scores} \
			--idr-plots ${sep=' ' idr_plots} \
			--idr-plots-pr ${sep=' ' idr_plots_pr} \
			${"--idr-plot-ppr " + idr_plot_ppr} \
			--frip-qcs ${sep=' ' frip_qcs} \
			--frip-qcs-pr1 ${sep=' ' frip_qcs_pr1} \
			--frip-qcs-pr2 ${sep=' ' frip_qcs_pr2} \
			${"--frip-qc-pooled " + frip_qc_pooled} \
			${"--frip-qc-ppr1 " + frip_qc_ppr1} \
			${"--frip-qc-ppr2 " + frip_qc_ppr2} \
			--frip-idr-qcs ${sep=' ' frip_idr_qcs} \
			--frip-idr-qcs-pr ${sep=' ' frip_idr_qcs_pr} \
			${"--frip-idr-qc-ppr " + frip_idr_qc_ppr} \
			--frip-overlap-qcs ${sep=' ' frip_overlap_qcs} \
			--frip-overlap-qcs-pr ${sep=' ' frip_overlap_qcs_pr} \
			${"--frip-overlap-qc-ppr " + frip_overlap_qc_ppr} \
			${"--idr-reproducibility-qc " + idr_reproducibility_qc} \
			${"--overlap-reproducibility-qc " + overlap_reproducibility_qc} \
			--out-qc-html qc.html \
			--out-qc-json qc.json
	}
	output {
		File report = glob('*qc.html')[0]
		File qc_json = glob('*qc.json')[0]
		#File encode_accession_json= glob('*encode_accession.json')[0]
	}
}

### workflow system tasks
task read_genome_tsv {
	File genome_tsv
	command {
		echo "Reading genome_tsv ${genome_tsv} ..."
	}
	output {
		Map[String,String] genome = read_map(genome_tsv)
	}
}

task pair_gen {
	Int num_rep
	command <<<
		python <<CODE
		for i in range(${num_rep}):
		    for j in range(i+1,${num_rep}):
		        print('{}\t{}'.format(i,j))
		CODE
	>>>
	output {
		Array[Array[Int]] pairs = if num_rep>1 then read_tsv(stdout()) else [[]]
	}
}
#@TASK_DEF_END
