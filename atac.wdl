# ENCODE DCC ATAC-Seq/DNase-Seq pipeline
# Author: Jin Lee (leepc12@gmail.com)

workflow atac {
	# mandatory pipeline type 
	String pipeline_type  		# atac or dnase
	String peak_caller = "macs2"

	# mandatory input files
	Array[Array[Array[String]]] fastqs 
								# [rep_id][merge_id][end_id]
								# 	after merging, it will reduce to 
								# 	[rep_id][end_id]
	Array[String] bams 			# [rep_id] if starting from bams
	Array[String] nodup_bams 	# [rep_id] if starting from filtered bams
	Array[String] tas 			# [rep_id] if starting from tag-aligns
	Array[String] peaks			# [rep_id] if starting from peaks
	Array[String] peaks_pr1		# [rep_id] if starting from peaks
	Array[String] peaks_pr2		# [rep_id] if starting from peaks
	File? peak_ppr1				# if starting from peaks
	File? peak_ppr2				# if starting from peaks
	File? peak_pooled			# if starting from peaks

	# mandaory adapters
	Array[Array[Array[String]]] adapters 
								# [rep_id][merge_id][end_id]
	# mandatory genome param
	String genome_tsv 		# reference genome data TSV file including
							# all important genome specific data file paths
							# and parameters
	Boolean paired_end 		# endedness of sample

	# optional but important
	Boolean? align_only 	# disable downstream analysis (peak calling, ...)
							# after alignment
	Boolean? true_rep_only 	# disable all analyses for pseudo replicates
							# naive-overlap and IDR will also be disabled
	Int? multimapping 		# multimapping reads

	# optional for MACS2
	Int? cap_num_peak 		# cap number of raw peaks called from MACS2
	Float? pval_thresh 		# p.value threshold
	Int? smooth_win 		# size of smoothing window
	Int? macs2_mem_mb 		# resource (memory in MB)
	Int? macs2_time_hr		# resource (walltime in hour)
	String? macs2_disks 	# resource disks for cloud platforms

	# optional for IDR
	Boolean? enable_idr		# enable IDR analysis on raw peaks
	Float? idr_thresh		# IDR threshold

	# OTHER IMPORTANT mandatory/optional parameters are declared in a task level

	# 1) determine input file type and num_rep (number of replicates)
	# 2) generate pairs of all replicates for later use
	# 3) read from genome_tsv
	call inputs {
		input :
			pipeline_type = pipeline_type,
			peak_caller = peak_caller,
			fastqs = fastqs,
			bams = bams,
			nodup_bams = nodup_bams,
			tas = tas,
			peaks = peaks,
			genome_tsv = genome_tsv,
			enable_idr_ = enable_idr,
			align_only_ = align_only,
			true_rep_only_ = true_rep_only,
	}

	# pipeline starts here (parallelized for each replicate)
	scatter(i in range(inputs.num_rep)) {
		if ( inputs.is_before_bam ) {
			# trim adapters and merge trimmed fastqs
			call trim_adapter { input :
				fastqs = fastqs[i],
				adapters = if length(adapters)>0 then adapters[i] else [],
				paired_end = paired_end,
			}
			# align trimmed/merged fastqs with bowtie2
			call bowtie2 { input :
				idx_tar = inputs.bowtie2_idx_tar,
				fastqs = trim_adapter.trimmed_merged_fastqs, #[R1,R2]
				paired_end = paired_end,
				multimapping = multimapping,
			}
		}
		if ( inputs.is_before_nodup_bam ) {
			# filter/dedup bam
			call filter { input :
				bam = if defined(bowtie2.bam) then bowtie2.bam else bams[i],
				paired_end = paired_end,
				multimapping = multimapping,
			}
		}
		if ( inputs.is_before_ta ) {
			# convert bam to tagalign and subsample it if necessary
			call bam2ta { input :
				bam = if defined(filter.nodup_bam) 
						then filter.nodup_bam else nodup_bams[i],
				disable_tn5_shift = inputs.disable_tn5_shift,
				paired_end = paired_end,
			}
		}
		if ( !inputs.align_only && inputs.is_before_peak ) {
			# subsample tagalign (non-mito) and cross-correlation analysis
			call xcor { input :
				ta = if defined(bam2ta.ta) then bam2ta.ta else tas[i],
				paired_end = paired_end,
			}
			# call peaks on tagalign
			call macs2 { input :
				ta = if defined(bam2ta.ta) then bam2ta.ta else tas[i],
				gensz = inputs.gensz,
				chrsz = inputs.chrsz,
				cap_num_peak = cap_num_peak,
				pval_thresh = pval_thresh,
				smooth_win = smooth_win,
				make_signal = true,
				blacklist = if inputs.has_blacklist then [inputs.blacklist] else [],
				mem_mb = macs2_mem_mb,
				disks = macs2_disks,
				time_hr = macs2_time_hr,
			}
		}
		if ( !inputs.align_only && inputs.is_before_peak && !inputs.true_rep_only ) {
			# make two self pseudo replicates per true replicate
			call spr { input :
				ta = if defined(bam2ta.ta) then bam2ta.ta else tas[i],
				paired_end = paired_end,
			}
		}
	}

	if ( !inputs.align_only && inputs.is_before_peak && inputs.num_rep>1 ) {
		# pool tagaligns from true replicates
		call pool_ta { input :
			tas = if defined(bam2ta.ta[0]) then bam2ta.ta else tas,
		}
		# call peaks on pooled replicate
		call macs2 as macs2_pooled { input :
			ta = pool_ta.ta_pooled,
			gensz = inputs.gensz,
			chrsz = inputs.chrsz,
			cap_num_peak = cap_num_peak,
			pval_thresh = pval_thresh,
			smooth_win = smooth_win,
			make_signal = true,
			blacklist = if inputs.has_blacklist then [inputs.blacklist] else [],
			mem_mb = macs2_mem_mb,
			disks = macs2_disks,
			time_hr = macs2_time_hr,
		}
	}
	if ( !inputs.align_only && !inputs.true_rep_only ) {		
		scatter(i in range(inputs.num_rep)) {
			if ( inputs.is_before_peak ) {
				# call peaks on 1st pseudo replicated tagalign 
				call macs2 as macs2_pr1 { input :
					ta = spr.ta_pr1[i],
					gensz = inputs.gensz,
					chrsz = inputs.chrsz,
					cap_num_peak = cap_num_peak,
					pval_thresh = pval_thresh,
					smooth_win = smooth_win,
					blacklist = if inputs.has_blacklist then [inputs.blacklist] else [],
					mem_mb = macs2_mem_mb,
					disks = macs2_disks,
					time_hr = macs2_time_hr,
				}
				call macs2 as macs2_pr2 { input :
					ta = spr.ta_pr2[i],
					gensz = inputs.gensz,
					chrsz = inputs.chrsz,
					cap_num_peak = cap_num_peak,
					pval_thresh = pval_thresh,
					smooth_win = smooth_win,
					blacklist = if inputs.has_blacklist then [inputs.blacklist] else [],
					mem_mb = macs2_mem_mb,
					disks = macs2_disks,
					time_hr = macs2_time_hr,
				}
			}
		}
		if ( inputs.is_before_peak && inputs.num_rep>1 ) {
			# pool tagaligns from pseudo replicates
			call pool_ta as pool_ta_pr1 { input :
				tas = spr.ta_pr1,
			}
			call pool_ta as pool_ta_pr2 { input :
				tas = spr.ta_pr2,
			}
			# call peaks on 1st pooled pseudo replicates
			call macs2 as macs2_ppr1 { input :
				ta = pool_ta_pr1.ta_pooled,
				gensz = inputs.gensz,
				chrsz = inputs.chrsz,
				cap_num_peak = cap_num_peak,
				pval_thresh = pval_thresh,
				smooth_win = smooth_win,
				blacklist = if inputs.has_blacklist then [inputs.blacklist] else [],
				mem_mb = macs2_mem_mb,
				disks = macs2_disks,
				time_hr = macs2_time_hr,
			}
			# call peaks on 2nd pooled pseudo replicates
			call macs2 as macs2_ppr2 { input :
				ta = pool_ta_pr2.ta_pooled,
				gensz = inputs.gensz,
				chrsz = inputs.chrsz,
				cap_num_peak = cap_num_peak,
				pval_thresh = pval_thresh,
				smooth_win = smooth_win,
				blacklist = if inputs.has_blacklist then [inputs.blacklist] else [],
				mem_mb = macs2_mem_mb,
				disks = macs2_disks,
				time_hr = macs2_time_hr,
			}
		}
	}

	# Naive overlap on every pair of true replicates
	if ( !inputs.align_only && inputs.num_rep>1 ) {
		scatter( pair in inputs.pairs ) {
			call overlap { input :
				prefix = "rep"+(pair[0]+1)+"-rep"+(pair[1]+1),
				peak1 = if inputs.is_before_peak then macs2.npeak[(pair[0])] 
						else peaks[(pair[0])],
				peak2 = if inputs.is_before_peak then macs2.npeak[(pair[1])]
						else peaks[(pair[1])],
				peak_pooled = if inputs.is_before_peak then macs2_pooled.npeak
						else peak_pooled,
				peak_type = inputs.peak_type,
				blacklist = if inputs.has_blacklist then [inputs.blacklist] else [],
				ta = if inputs.is_before_peak then [pool_ta.ta_pooled] else [],
			}
		}
	}
	if ( !inputs.align_only && !inputs.true_rep_only ) {
		# Naive overlap on pseduo replicates
		scatter( i in range(inputs.num_rep) ) {
			call overlap as overlap_pr { input : 
				prefix = "rep"+(i+1)+"-pr",
				peak1 = if inputs.is_before_peak then macs2_pr1.npeak[i]
						else peaks_pr1[i],
				peak2 = if inputs.is_before_peak then macs2_pr2.npeak[i]
						else peaks_pr2[i],
				peak_pooled = if inputs.is_before_peak then macs2.npeak[i]
						else peak_pooled,
				peak_type = inputs.peak_type,
				blacklist = if inputs.has_blacklist then [inputs.blacklist] else [],
				ta = if inputs.is_before_ta then [bam2ta.ta[i]]
						else if inputs.is_before_peak then [tas[i]]
						else [],
			}
		}
	}
	if ( !inputs.align_only && !inputs.true_rep_only && inputs.num_rep>1 ) {
		# Naive overlap on pooled pseudo replicates
		call overlap as overlap_ppr { input : 
			prefix = "ppr",
			peak1 = if inputs.is_before_peak then macs2_ppr1.npeak
					else peak_ppr1,
			peak2 = if inputs.is_before_peak then macs2_ppr2.npeak
					else peak_ppr2,
			peak_pooled = if inputs.is_before_peak then macs2_pooled.npeak
					else peak_pooled,
			peak_type = inputs.peak_type,
			blacklist = if inputs.has_blacklist then [inputs.blacklist] else [],
			ta = if inputs.is_before_peak then [pool_ta.ta_pooled]
					else [],
		}
	}
	if ( !inputs.align_only && !inputs.true_rep_only ) {
		# reproducibility QC for overlapping peaks
		call reproducibility as reproducibility_overlap { input :
			prefix = 'overlap',
			peaks = if inputs.num_rep>1	then overlap.bfilt_overlap_peak else [],
			peaks_pr = overlap_pr.bfilt_overlap_peak,
			peak_ppr = overlap_ppr.bfilt_overlap_peak,
		}
	}

	if ( !inputs.align_only && inputs.num_rep>1 && inputs.enable_idr ) {
		# IDR on every pair of true replicates
		scatter( pair in inputs.pairs ) {
			call idr { input : 
				prefix = "rep"+(pair[0]+1)+"-rep"+(pair[1]+1),
				peak1 = if inputs.is_before_peak then macs2.npeak[(pair[0])]
						else peaks[(pair[0])],
				peak2 = if inputs.is_before_peak then macs2.npeak[(pair[1])]
						else peaks[(pair[1])],
				peak_pooled = if inputs.is_before_peak then macs2_pooled.npeak 
						else peak_pooled,
				idr_thresh = select_first([idr_thresh,0.1]),
				peak_type = inputs.peak_type,
				blacklist = if inputs.has_blacklist then [inputs.blacklist] else [],
				ta = if inputs.is_before_peak then [pool_ta.ta_pooled] else [],
			}
		}
	}
	if ( !inputs.align_only && !inputs.true_rep_only && inputs.enable_idr ) {
		# IDR on pseduo replicates
		scatter( i in range(inputs.num_rep) ) {
			call idr as idr_pr { input : 
				prefix = "rep"+(i+1)+"-pr",
				peak1 = if inputs.is_before_peak then macs2_pr1.npeak[i]
						else peaks_pr1[i],
				peak2 = if inputs.is_before_peak then macs2_pr2.npeak[i]
						else peaks_pr2[i],
				peak_pooled = if inputs.is_before_peak then macs2.npeak[i]
						else peak_pooled,
				idr_thresh = select_first([idr_thresh,0.1]),
				peak_type = inputs.peak_type,
				blacklist = if inputs.has_blacklist then [inputs.blacklist] else [],
				ta = if inputs.is_before_ta then [bam2ta.ta[i]]
						else if inputs.is_before_peak then [tas[i]]
						else [],
			}
		}
	}
	if ( !inputs.align_only && !inputs.true_rep_only && inputs.num_rep>1 && inputs.enable_idr ) {
		# IDR on pooled pseduo replicates
		call idr as idr_ppr { input : 
			prefix = "ppr",
			peak1 = if inputs.is_before_peak then macs2_ppr1.npeak
					else peak_ppr1,
			peak2 = if inputs.is_before_peak then macs2_ppr2.npeak
					else peak_ppr2,
			peak_pooled = if inputs.is_before_peak then macs2_pooled.npeak
					else peak_pooled,
			idr_thresh = select_first([idr_thresh,0.1]),
			peak_type = inputs.peak_type,
			blacklist = if inputs.has_blacklist then [inputs.blacklist] else [],
			ta = if inputs.is_before_peak then [pool_ta.ta_pooled] else [],
		}
	}
	if ( !inputs.align_only && !inputs.true_rep_only && inputs.enable_idr ) {
		# reproducibility QC for IDR peaks
		call reproducibility as reproducibility_idr { input :
			prefix = 'idr',
			peaks = if inputs.num_rep>1	then idr.bfilt_idr_peak else [],
			peaks_pr = idr_pr.bfilt_idr_peak,
			peak_ppr = idr_ppr.bfilt_idr_peak,
		}
	}

	call qc_report { input :
		paired_end = paired_end,
		peak_caller = peak_caller,
		idr_thresh = select_first([idr_thresh,0.1]),
		flagstat_qcs = if inputs.is_before_bam 
						then bowtie2.flagstat_qc else [],
		nodup_flagstat_qcs = if inputs.is_before_nodup_bam 
						then filter.flagstat_qc	else [],
		dup_qcs = if inputs.is_before_nodup_bam
						then filter.dup_qc else [],
		pbc_qcs = if inputs.is_before_nodup_bam
						then filter_pbc_qc else [],
		xcor_plots = if !inputs.align_only && inputs.is_before_peak
						then xcor.plot_png else [],
		xcor_scores = if !inputs.align_only && inputs.is_before_peak
						then xcor.score else [],
		frip_qcs = if !inputs.align_only && inputs.is_before_peak
						then macs2.frip_qc else [],
		frip_qcs_pr1 = if !inputs.align_only && inputs.is_before_peak && !inputs.true_rep_only 
						then macs2_pr1.frip_qc else [],
		frip_qcs_pr2 = if !inputs.align_only && inputs.is_before_peak && !inputs.true_rep_only 
						then macs2_pr2.frip_qc else [],
		frip_qc_pooled = if !inputs.align_only && inputs.is_before_peak && inputs.num_rep>1 
						then [macs2_pooled.frip_qc] else [],
		frip_qc_ppr1 = if !inputs.align_only && inputs.is_before_peak && !inputs.true_rep_only && inputs.num_rep>1 
						then [macs2_ppr1.frip_qc] else [],
		frip_qc_ppr2 = if !inputs.align_only && inputs.is_before_peak && !inputs.true_rep_only && inputs.num_rep>1 
						then [macs2_ppr2.frip_qc] else [],
		idr_plots = if !inputs.align_only && inputs.num_rep>1 && inputs.enable_idr 
						then idr.idr_plot else [],
		idr_plots_pr = if !inputs.align_only && !inputs.true_rep_only && inputs.enable_idr
						then idr_pr.idr_plot else [],
		idr_plot_ppr = if !inputs.align_only && !inputs.true_rep_only && inputs.num_rep>1 && inputs.enable_idr
						then [idr_ppr.idr_plot] else [],
		frip_idr_qcs = if !inputs.align_only && inputs.is_before_peak && inputs.num_rep>1 && inputs.enable_idr
						then idr.frip_qc else [],
		frip_idr_qcs_pr = if !inputs.align_only && inputs.is_before_peak && !inputs.true_rep_only && inputs.enable_idr
						then idr_pr.frip_qc else [],
		frip_idr_qc_ppr = if !inputs.align_only && inputs.is_before_peak && !inputs.true_rep_only && inputs.num_rep>1 && inputs.enable_idr
						then [idr_ppr.frip_qc] else [],
		frip_overlap_qcs = if !inputs.align_only && inputs.is_before_peak && inputs.num_rep>1
						then overlap.frip_qc else [],
		frip_overlap_qcs_pr = if !inputs.align_only && inputs.is_before_peak && !inputs.true_rep_only
						then overlap_pr.frip_qc else [],
		frip_overlap_qc_ppr = if !inputs.align_only && inputs.is_before_peak && !inputs.true_rep_only && inputs.num_rep>1
						then [overlap_ppr.frip_qc] else [],
		idr_reproducibility_qc =
				if !inputs.align_only && !inputs.true_rep_only && inputs.enable_idr &&
					defined(reproducibility_idr.reproducibility_qc)
						then [reproducibility_idr.reproducibility_qc] else [],
		overlap_reproducibility_qc = 
				if !inputs.align_only && !inputs.true_rep_only && 
					defined(reproducibility_overlap.reproducibility_qc)
						then [reproducibility_overlap.reproducibility_qc] else [],
	}

	output {
		qc_report.report
		qc_report.qc_json
	}
}

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
			${"--adapters " + write_tsv(adapters)} \
			${if paired_end then "--paired-end" else ""} \
			${if auto_detect_adapter then "--auto-detect-adapter" else ""} \
			${"--min-trim-len " + min_trim_len} \
			${"--err-rate " + err_rate} \
			${"--nth " + select_first([cpu,4])}
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
								# dup.qc and pbc.qc will be emptry files
								# and nodup_bam in the output is 
								# filtered bam with dupes
	# resource
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
		touch null
	}
	output {
		File nodup_bam = glob("*.bam")[0]
		File nodup_bai = glob("*.bai")[0]
		File flagstat_qc = glob("*.flagstat.qc")[0]
		# optional (if no_dup_removal then empty qc files)
		# need temporary boolean in output (cromwell if bug in output)
		Boolean no_dup_removal_ = select_first([no_dup_removal,false])
		File dup_qc = if no_dup_removal_ then glob("null")[0] else glob("*.dup.qc")[0]
		File pbc_qc = if no_dup_removal_ then glob("null")[0] else glob("*.pbc.qc")[0]
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

	command {
		python $(which encode_spr.py) \
			${ta} \
			${if paired_end then "--paired-end" else ""}
	}
	output {
		File ta_pr1 = glob("*.pr1.tagAlign.gz")[0]
		File ta_pr2 = glob("*.pr2.tagAlign.gz")[0]
	}
}

task pool_ta {
	# parameters from workflow
	Array[File?] tas

	command {
		python $(which encode_pool_ta.py) \
			${sep=' ' tas}
	}
	output {
		File ta_pooled = glob("*.tagAlign.gz")[0]
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
			${"--nth " + cpu}
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
	File? ta
	String gensz		# Genome size (sum of entries in 2nd column of 
                        # chr. sizes file, or hs for human, ms for mouse)
	File chrsz			# 2-col chromosome sizes file
	Int? cap_num_peak	# cap number of raw peaks called from MACS2
	Float? pval_thresh	# p.value threshold
	Int? smooth_win		# size of smoothing window
	Boolean? make_signal
	Array[File] blacklist 	# blacklist BED to filter raw peaks
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
			${"--p-val-thresh "+ pval_thresh} \
			${"--smooth-win "+ smooth_win} \
			${if select_first([make_signal,false]) then "--make-signal" else ""} \
			${if length(blacklist)>0 then "--blacklist "+ blacklist[0] else ""}
		touch null
	}
	output {
		File npeak = glob("*[!.][!b][!f][!i][!l][!t]."+peak_type+".gz")[0]
		File bfilt_npeak = if length(blacklist)>0 then 
						glob("*.bfilt."+peak_type+".gz")[0] else npeak
		# need temporary boolean in output (cromwell if bug in output)
		Boolean make_signal_ = select_first([make_signal,false])
		File sig_pval = if make_signal_ then glob("*.pval.signal.bigwig")[0] else glob("null")[0]
		File sig_fc = if make_signal_ then glob("*.fc.signal.bigwig")[0] else glob("null")[0]
		File frip_qc = glob("*.frip.qc")[0]
	}
	runtime {
		memory : "${select_first([mem_mb,'16000'])} MB"
		time : select_first([time_hr,24])
		disks : select_first([disks,"local-disk 100 HDD"])
	}
}

task idr {
	# parameters from workflow
	String? prefix 		# prefix for IDR output file
	File? peak1 			
	File? peak2
	File? peak_pooled
	Float? idr_thresh
	Array[File] blacklist 	# blacklist BED to filter raw peaks
	Array[File?] ta		# to calculate FRiP
	String peak_type

	command {
		python $(which encode_idr.py) \
			${peak1} ${peak2} ${peak_pooled} \
			${"--prefix " + prefix} \
			${"--idr-thresh " + idr_thresh} \
			${"--peak-type " + peak_type} \
			--idr-rank p.value \
			${if length(blacklist)>0 then "--blacklist "+ blacklist[0] else ""} \
			${if length(ta)>0 then "--ta "+ ta[0] else ""}
		touch null
	}
	output {
		File idr_peak = glob("*[!.][!b][!f][!i][!l][!t]."+peak_type+".gz")[0]
		File bfilt_idr_peak = if length(blacklist)>0 then 
							glob("*.bfilt."+peak_type+".gz")[0] else idr_peak
		File idr_plot = glob("*.txt.png")[0]
		File idr_unthresholded_peak = glob("*.txt.gz")[0]
		File idr_log = glob("*.log")[0]
		# need temporary boolean in output (cromwell if bug in output)
		#Boolean has_ta = length(ta)>0
		#File frip_qc = if has_ta then glob("*.frip.qc")[0] else glob("null")[0]
		File frip_qc = if length(ta)>0 then glob("*.frip.qc")[0] else glob("null")[0]
	}
}

task overlap {
	# parameters from workflow
	String prefix 		# prefix for IDR output file
	File? peak1
	File? peak2
	File? peak_pooled
	Array[File] blacklist 	# blacklist BED to filter raw peaks
	Array[File?] ta 		# to calculate FRiP	
	String peak_type

	command {
		python $(which encode_naive_overlap.py) \
			${peak1} ${peak2} ${peak_pooled} \
			${"--prefix " + prefix} \
			${"--peak-type " + peak_type} \
			${if length(blacklist)>0 then "--blacklist "+ blacklist[0] else ""} \
			${if length(ta)>0 then "--ta "+ ta[0] else ""}
		touch null
	}
	output {
		File overlap_peak = glob("*[!.][!b][!f][!i][!l][!t]."+peak_type+".gz")[0]
		File bfilt_overlap_peak = if length(blacklist)>0 then 
							glob("*.bfilt."+peak_type+".gz")[0] else overlap_peak
		# need temporary boolean in output (cromwell if bug in output)
		#Boolean has_ta = length(ta)>0
		#File frip_qc = if has_ta then glob("*.frip.qc")[0] else glob("null")[0]
		File frip_qc = if length(ta)>0 then glob("*.frip.qc")[0] else glob("null")[0]
	}
}

task reproducibility {
	# parameters from workflow
	String prefix
	Array[File?] peaks # peak files from pair of true replicates
						# in a sorted order. for example of 4 replicates,
						# 1,2 1,3 1,4 2,3 2,4 3,4.
                        # x,y means peak file from rep-x vs rep-y
	Array[File?] peaks_pr	# peak files from pseudo replicates
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
	String peak_caller
	Float idr_thresh
	# QCs
	Array[File] flagstat_qcs
	Array[File] nodup_flagstat_qcs
	Array[File] dup_qcs
	Array[File] pbc_qcs
	Array[File?] xcor_plots
	Array[File?] xcor_scores
	Array[File?] idr_plots
	Array[File?] idr_plots_pr
	Array[File?] idr_plot_ppr # actually not an array
	Array[File?] frip_qcs
	Array[File?] frip_qcs_pr1
	Array[File?] frip_qcs_pr2
	Array[File?] frip_qc_pooled # actually not an array
	Array[File?] frip_qc_ppr1 # actually not an array
	Array[File?] frip_qc_ppr2 # actually not an array
	Array[File?] frip_idr_qcs
	Array[File?] frip_idr_qcs_pr
	Array[File?] frip_idr_qc_ppr # actually not an array
	Array[File?] frip_overlap_qcs
	Array[File?] frip_overlap_qcs_pr
	Array[File?] frip_overlap_qc_ppr # actually not an array
	Array[File?] idr_reproducibility_qc # actually not an array
	Array[File?] overlap_reproducibility_qc # actually not an array

	command {
		python $(which encode_qc_report.py) \
			${"--name '" + name + "'"} \
			${"--desc '" + desc + "'"} \
			${if paired_end then "--paired-end" else ""} \
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
			--idr-plot-ppr ${sep=' ' idr_plot_ppr} \
			--frip-qcs ${sep=' ' frip_qcs} \
			--frip-qcs-pr1 ${sep=' ' frip_qcs_pr1} \
			--frip-qcs-pr2 ${sep=' ' frip_qcs_pr2} \
			--frip-qc-pooled ${sep=' ' frip_qc_pooled} \
			--frip-qc-ppr1 ${sep=' ' frip_qc_ppr1} \
			--frip-qc-ppr1 ${sep=' ' frip_qc_ppr2} \
			--frip-idr-qcs ${sep=' ' frip_idr_qcs} \
			--frip-idr-qcs-pr ${sep=' ' frip_idr_qcs_pr} \
			--frip-idr-qc-ppr ${sep=' ' frip_idr_qc_ppr} \
			--frip-overlap-qcs ${sep=' ' frip_overlap_qcs} \
			--frip-overlap-qcs-pr ${sep=' ' frip_overlap_qcs_pr} \
			--frip-overlap-qc-ppr ${sep=' ' frip_overlap_qc_ppr} \
			--idr-reproducibility-qc ${sep=' ' idr_reproducibility_qc} \
			--overlap-reproducibility-qc ${sep=' ' overlap_reproducibility_qc} \
			--out-qc-html qc.html \
			--out-qc-json qc.json
	}
}

### workflow system tasks

# to reduce overhead of provisioning extra nodes
# we have only one task to
# 	1) determine input type and number of replicates	
# 	2) generate pair (rep-x_vs_rep-y) of all true replicate
# 	3) read genome_tsv and get ready to download files 
# 		On Google Cloud Platform 
#		files are downloaded from gs://atac-seq-pipeline-genome-data/
#	4) have output booleans for optional flags in workflow
#		to simplify if-else statements
task inputs {
	# parameters from workflow
	String pipeline_type
	String peak_caller
	Array[Array[Array[String]]] fastqs 
	Array[String] bams
	Array[String] nodup_bams
	Array[String] tas
	Array[String] peaks
	File genome_tsv
	Boolean? enable_idr_
	Boolean? align_only_
	Boolean? true_rep_only_

	command <<<
		python <<CODE
		name = ['fastq','bam','nodup_bam','ta','peak']
		arr = [${length(fastqs)},${length(bams)},
		       ${length(nodup_bams)},${length(tas)},
		       ${length(peaks)}]
		num_rep = max(arr)
		type = name[arr.index(num_rep)]
		with open('num_rep.txt','w') as fp:
		    fp.write(str(num_rep)) 
		with open('type.txt','w') as fp:
		    fp.write(type)		    
		for i in range(num_rep):
		    for j in range(i+1,num_rep):
		        print('{}\t{}'.format(i,j))
		CODE
	>>>
	output {
		# peak file type
		String peak_type = if peak_caller=='macs2' then 'narrowPeak'
						else if peak_caller=='spp' then 'regionPeak'
						else 'narrowPeak'
		# read genome TSV
		Map[String,String] genome = read_map(genome_tsv)
		String ref_fa = genome['ref_fa']
		String bowtie2_idx_tar = genome['bowtie2_idx_tar']
		#String bwa_idx_tar = genome['bwa_idx_tar']
		String blacklist = genome['blacklist']
		String chrsz = genome['chrsz']
		String gensz = genome['gensz']

		# input types and usefule booleans
		String type = read_string("type.txt")
		Int num_rep = read_int("num_rep.txt")
		Boolean is_before_bam =
			type=='fastq'
		Boolean is_before_nodup_bam =
			type=='fastq' || type=='bam'
		Boolean is_before_ta =
			type=='fastq' || type=='bam' ||	type=='nodup_bam'
		Boolean is_before_peak = 
			type=='fastq' || type=='bam' ||	type=='nodup_bam' || type=='ta'
		Boolean is_peak = type=='peak'
		Boolean	enable_idr = select_first([enable_idr_,false])
		Boolean align_only = select_first([align_only_,false])
		Boolean true_rep_only = select_first([true_rep_only_,false])
		Boolean disable_tn5_shift = pipeline_type!='atac'
		Boolean has_blacklist = blacklist!='/dev/null'

		# pair of true replicates
		Array[Array[Int]] pairs = if num_rep>1 then read_tsv(stdout()) else [[]]
	}
}
 