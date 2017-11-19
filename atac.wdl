# ENCODE DCC ATAC-Seq/DNase-Seq pipeline
# Author: Jin Lee (leepc12@gmail.com)

workflow atac {
	# mandatory input files
	Array[Array[Array[String]]] fastqs 
								# [rep_id][merge_id][end_id]
								# 	after merging, it will reduce to 
								# 	[rep_id][end_id]
	Array[String] bams 			# [rep_id] if starting from bams
	Array[String] nodup_bams 		# [rep_id] if starting from filtered bams
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

	# optional metadata
 	String? name 			# name of sample
	String? desc 			# description for sample
	String? accession_id 	# ENCODE accession ID of sample

	# OTHER IMPORTANT mandatory/optional parameters are declared in a task level

	# 1) determine input file type and num_rep (number of replicates)
	# 2) generate pairs of all replicates for later use
	# 3) read from genome_tsv
	call inputs {
		input :
			fastqs = fastqs,
			bams = bams,
			nodup_bams = nodup_bams,
			tas = tas,
			peaks = peaks,
			genome_tsv = genome_tsv,
	}

	# pipeline starts here (parallelized for each replicate)
	scatter(i in range(inputs.num_rep)) {
		if ( inputs.is_before_bam ) {
			# trim adapters and merge trimmed fastqs
			call trim_adapter {
				input:
					fastqs = fastqs[i],
					adapters = if length(adapters)>0 
							then adapters[i] else [],
					paired_end = paired_end,
			}
			# align trimmed/merged fastqs with bowtie2
			call bowtie2 {
				input:
					idx_tar = inputs.bowtie2_idx_tar,
					fastqs = trim_adapter.trimmed_merged_fastqs, #[R1,R2]
					paired_end = paired_end,
					multimapping = multimapping,
			}
		}
		if ( inputs.is_before_nodup_bam ) {
			# filter/dedup bam
			call filter {
				input:
					bam = if defined(bowtie2.bam) 
							then bowtie2.bam else bams[i],
					paired_end = paired_end,
					multimapping = multimapping,
			}
		}
		if ( inputs.is_before_ta ) {
			# convert bam to tagalign and subsample it if necessary
			call bam2ta {
				input:
					bam = if defined(filter.nodup_bam) 
							then filter.nodup_bam else nodup_bams[i],
					paired_end = paired_end,
			}
		}
		if ( inputs.is_before_peak ) {
			# subsample tagalign (non-mito) and cross-correlation analysis
			call xcor {
				input:
					ta = if defined(bam2ta.ta) 
							then bam2ta.ta else tas[i],
					paired_end = paired_end,
			}
			# call peaks on tagalign
			call macs2 {
				input:
					ta = if defined(bam2ta.ta) 
							then bam2ta.ta else tas[i],
					gensz = inputs.gensz,
					chrsz = inputs.chrsz,
					cap_num_peak = cap_num_peak,
					pval_thresh = pval_thresh,
					smooth_win = smooth_win,
					make_signal = true,
					mem_mb = macs2_mem_mb,
					disks = macs2_disks,
					time_hr = macs2_time_hr,
			}
		}
		if ( !select_first([true_rep_only,false]) && 
			inputs.is_before_peak) {
			# make two self pseudo replicates per true replicate
			call spr {
				input:
					ta = if defined(bam2ta.ta) 
							then bam2ta.ta else tas[i],
					paired_end = paired_end,
			}
		}
		# filter out peaks with blacklist
		call blacklist_filter as bfilt_macs2 {
			input:
				peak = if defined(macs2.npeak)
						then macs2.npeak else peaks[i],
				blacklist = inputs.blacklist,				
		}
	}

	if ( inputs.num_rep>1 ) {
		if (is_before_peak) {
			# pool tagaligns from true replicates
			call pool_ta {
				input :
					tas = if defined(bam2ta.ta[0])
							then bam2ta.ta else tas,
			}
			# call peaks on pooled replicate
			call macs2 as macs2_pooled {
				input:
					ta = pool_ta.ta_pooled,
					gensz = inputs.gensz,
					chrsz = inputs.chrsz,
					cap_num_peak = cap_num_peak,
					pval_thresh = pval_thresh,
					smooth_win = smooth_win,
					make_signal = true,
					mem_mb = macs2_mem_mb,
					disks = macs2_disks,
					time_hr = macs2_time_hr,
			}
		}
		call blacklist_filter as bfilt_macs2_pooled {
			input:
				peak = if defined(macs2_pooled.npeak)
						then macs2_pooled.npeak else peak_pooled,
				blacklist = inputs.blacklist,
		}
	}
	if ( !select_first([true_rep_only,false]) ) {		
		scatter(i in range(inputs.num_rep)) {
			if ( inputs.is_before_peak ) {
				# call peaks on 1st pseudo replicated tagalign 
				call macs2 as macs2_pr1 {
					input:
						ta = spr.ta_pr1[i],
						gensz = inputs.gensz,
						chrsz = inputs.chrsz,
						cap_num_peak = cap_num_peak,
						pval_thresh = pval_thresh,
						smooth_win = smooth_win,
						mem_mb = macs2_mem_mb,
						disks = macs2_disks,
						time_hr = macs2_time_hr,
				}
				call macs2 as macs2_pr2 {
					input:
						ta = spr.ta_pr2[i],
						gensz = inputs.gensz,
						chrsz = inputs.chrsz,
						cap_num_peak = cap_num_peak,
						pval_thresh = pval_thresh,
						smooth_win = smooth_win,
						mem_mb = macs2_mem_mb,
						disks = macs2_disks,
						time_hr = macs2_time_hr,
				}
			}
			call blacklist_filter as bfilt_macs2_pr1 {
				input:
					peak = if defined(macs2_pr1.npeak)
						then macs2_pr1.npeak else peaks_pr1[i],
					blacklist = inputs.blacklist,
			}
			call blacklist_filter as bfilt_macs2_pr2 {
				input:
					peak = if defined(macs2_pr2.npeak)
						then macs2_pr2.npeak else peaks_pr2[i],
					blacklist = inputs.blacklist,
			}
		}
		if ( inputs.num_rep>1 ) {
			if ( inputs.is_before_peak ) {
				# pool tagaligns from pseudo replicates
				call pool_ta as pool_ta_pr1 {
					input :
						tas = spr.ta_pr1,
				}
				call pool_ta as pool_ta_pr2 {
					input :
						tas = spr.ta_pr2,
				}
				# call peaks on 1st pooled pseudo replicates
				call macs2 as macs2_ppr1 {
					input:
						ta = pool_ta_pr1.ta_pooled,
						gensz = inputs.gensz,
						chrsz = inputs.chrsz,
						cap_num_peak = cap_num_peak,
						pval_thresh = pval_thresh,
						smooth_win = smooth_win,
						mem_mb = macs2_mem_mb,
						disks = macs2_disks,
						time_hr = macs2_time_hr,
				}
				# call peaks on 2nd pooled pseudo replicates
				call macs2 as macs2_ppr2 {
					input:
						ta = pool_ta_pr2.ta_pooled,
						gensz = inputs.gensz,
						chrsz = inputs.chrsz,
						cap_num_peak = cap_num_peak,
						pval_thresh = pval_thresh,
						smooth_win = smooth_win,
						mem_mb = macs2_mem_mb,
						disks = macs2_disks,
						time_hr = macs2_time_hr,
				}
			}
			call blacklist_filter as bfilt_macs2_ppr1 {
				input:
					peak = if defined(macs2_ppr1.npeak)
							then macs2_ppr1.npeak else peak_ppr1,
					blacklist = inputs.blacklist,
			}
			call blacklist_filter as bfilt_macs2_ppr2 {
				input:
					peak = if defined(macs2_ppr2.npeak)
							then macs2_ppr2.npeak else peak_ppr2,
					blacklist = inputs.blacklist,
			}
		}
	}

	# Naive overlap on every pair of true replicates
	if ( inputs.num_rep>1 ) {
		scatter( pair in inputs.pairs ) {
			call overlap {
				input :
					prefix = "rep"+(pair[0]+1)+
							"-rep"+(pair[1]+1),
					peak1 = if defined(macs2.npeak[0])
						then macs2.npeak[(pair[0])] 
						else peaks[(pair[0])],
					peak2 = if defined(macs2.npeak[0])
						then macs2.npeak[(pair[1])]
						else peaks[(pair[1])],
					peak_pooled = if defined(macs2_pooled.npeak)
						then macs2_pooled.npeak
						else peak_pooled,
			}
			call blacklist_filter as bfilt_overlap {
				input:
					peak = overlap.overlap_peak,
					blacklist = inputs.blacklist,
			}
		}
	}
	if ( !select_first([true_rep_only,false]) ) {
		# Naive overlap on pseduo replicates
		scatter( i in range(inputs.num_rep) ) {
			call overlap as overlap_pr {
				input : 
					prefix = "rep"+(i+1)+"-pr",
					peak1 = if defined(macs2_pr1.npeak[0])
						then macs2_pr1.npeak[i]
						else peaks_pr1[i],
					peak2 = if defined(macs2_pr2.npeak[0])
						then macs2_pr2.npeak[i]
						else peaks_pr2[i],
					peak_pooled = if defined(macs2.npeak[i])
						then macs2.npeak[i]
						else peak_pooled,
			}
			call blacklist_filter as bfilt_overlap_pr {
				input:
					peak = overlap_pr.overlap_peak,
					blacklist = inputs.blacklist,
			}
		}
		if ( inputs.num_rep>1 ) {
			# Naive overlap on pooled pseudo replicates
			call overlap as overlap_ppr {
				input : 
					prefix = "ppr",
					peak1 = if defined(macs2_ppr1.npeak)
						then macs2_ppr1.npeak
						else peak_ppr1,
					peak2 = if defined(macs2_ppr2.npeak)
						then macs2_ppr2.npeak
						else peak_ppr2,
					peak_pooled = if defined(macs2_pooled.npeak)
						then macs2_pooled.npeak
						else peak_pooled,
			}
			call blacklist_filter as bfilt_overlap_ppr {
				input:
					peak = overlap_ppr.overlap_peak,
					blacklist = inputs.blacklist,
			}
		}
		# reproducibility QC for overlapping peaks
		call reproducibility as reproducibility_overlap {
			input:
				prefix = 'overlap',
				peaks = if defined(bfilt_overlap.filtered_peak[0])
				then bfilt_overlap.filtered_peak else [],
				peaks_pr = bfilt_overlap_pr.filtered_peak,
				peak_ppr = bfilt_overlap_ppr.filtered_peak,
		}
	}

	if ( select_first([enable_idr,false]) ) {
		if ( inputs.num_rep>1 ) {
			scatter( pair in inputs.pairs ) {
				# IDR on every pair of true replicates
				call idr {
					input : 
						prefix = "rep"+(pair[0]+1)
								+"-rep"+(pair[1]+1),
						peak1 = if defined(macs2.npeak[0])
							then macs2.npeak[(pair[0])]
							else peaks[(pair[0])],
						peak2 = if defined(macs2.npeak[0])
							then macs2.npeak[(pair[1])]
							else peaks[(pair[1])],
						peak_pooled = if defined(macs2_pooled.npeak)
							then macs2_pooled.npeak 
							else peak_pooled,
						idr_thresh = idr_thresh,
				}
				call blacklist_filter as bfilt_idr {
					input:
						peak = idr.idr_peak,
						blacklist = inputs.blacklist,
				}
			}
		}
		if ( !select_first([true_rep_only,false]) ) {
			# IDR on pseduo replicates
			scatter( i in range(inputs.num_rep) ) {
				call idr as idr_pr {
					input : 
						prefix = "rep"+(i+1)+"-pr",
						peak1 = if defined(macs2_pr1.npeak[i])
							then macs2_pr1.npeak[i]
							else peaks_pr1[i],
						peak2 = if defined(macs2_pr2.npeak[i])
							then macs2_pr2.npeak[i]
							else peaks_pr2[i],
						peak_pooled = if defined(macs2.npeak[i])
							then macs2.npeak[i]
							else peak_pooled,
						idr_thresh = idr_thresh,
				}
				call blacklist_filter as bfilt_idr_pr {
					input:
						peak = idr_pr.idr_peak,
						blacklist = inputs.blacklist,
				}
			}
			if ( inputs.num_rep>1 ) {
				call idr as idr_ppr {
					input : 
						prefix = "ppr",
						peak1 = if defined(macs2_ppr1.npeak)
							then macs2_ppr1.npeak
							else peak_ppr1,
						peak2 = if defined(macs2_ppr2.npeak)
							then macs2_ppr2.npeak
							else peak_ppr2,
						peak_pooled = if defined(macs2_pooled.npeak)
							then macs2_pooled.npeak
							else peak_pooled,
						idr_thresh = idr_thresh,
				}
				call blacklist_filter as bfilt_idr_ppr {
					input:
						peak = idr_ppr.idr_peak,
						blacklist = inputs.blacklist,
				}
			}
			# reproducibility QC for IDR peaks
			call reproducibility as reproducibility_idr {
				input:
					prefix = 'idr',
					peaks = if defined(bfilt_idr.filtered_peak[0])
						then bfilt_idr.filtered_peak else [],
					peaks_pr = bfilt_idr_pr.filtered_peak,
					peak_ppr = bfilt_idr_ppr.filtered_peak,
			}
		}
	}
}

# genomic tasks
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
			${if auto_detect_adapter
				then "--auto-detect-adapter" else ""} \
			${"--min-trim-len " + min_trim_len} \
			${"--err-rate " + err_rate} \
			${"--nth " + select_first([cpu,4])}
	}
	output {
		Array[File] trimmed_merged_fastqs = if paired_end then 
			[glob("*.R1.fastq.gz")[0], glob("*.R2.fastq.gz")[0]]
			else glob("*.R1.fastq.gz")
	}
	runtime {
		cpu : select_first([cpu,4])
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
			${if select_first([no_dup_removal,false])
				then "--no-dup-removal" else ""} \
			${"--nth " + cpu}			
	}
	output {
		File nodup_bam = glob("*.bam")[0]
		File nodup_bai = glob("*.bai")[0]
		File flagstat_qc = glob("*.flagstat.qc")[0]
		# optional (if no_dup_removall then empty qc files)		
		File dup_qc = glob("*.dup.qc")[0]
		File pbc_qc = glob("*.pbc.qc")[0]
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
	# optional
	Boolean? disable_tn5_shift 	# no tn5 shifting (it's for dnase-seq)
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
			${if select_first([disable_tn5_shift,false])
				then "--disable-tn5-shift" else ""} \
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
			${"--subsample " + select_first(
								[subsample,25000000])} \
			--speak=0 \
			${"--nth " + cpu}
	}
	output {
		File plot = glob("*.cc.plot.pdf")[0]
		File score = glob("*.cc.qc")[0]
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
	# resource
	Int? mem_mb
	Int? time_hr
	String? disks

	command {
		python $(which encode_macs2.py) \
			${ta} \
			${"--gensz "+ gensz} \
			${"--chrsz " + chrsz} \
			${"--cap-num-peak " + select_first(
								[cap_num_peak,300000])} \
			${"--p-val-thresh "+ pval_thresh} \
			${"--smooth-win "+ smooth_win} \
			${if select_first([make_signal,false])
				then "--make-signal" else ""}
	}
	output {
		File npeak = glob("*.narrowPeak.gz")[0]
		# optional (if not make_signal then empty signal files)
		File sig_pval = glob("*.pval.signal.bigwig")[0]
		File sig_fc = glob("*.fc.signal.bigwig")[0]
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

	command {
		python $(which encode_idr.py) \
			${peak1} ${peak2} ${peak_pooled} \
			${"--prefix " + prefix} \
			${"--idr-thresh " + idr_thresh} \
			--idr-rank p.value
	}
	output {
		File idr_peak = glob("*eak.gz")[0]
		File idr_plot = glob("*.txt.png")[0]
		File idr_unthresholded_peak = glob("*.txt.gz")[0]
		File idr_log = glob("*.log")[0]
	}
}

task overlap {
	# parameters from workflow
	String? prefix 		# prefix for IDR output file
	File? peak1
	File? peak2
	File? peak_pooled

	command {
		python $(which encode_naive_overlap.py) \
			${peak1} ${peak2} ${peak_pooled} \
			${"--prefix " + prefix}
	}
	output {
		File overlap_peak = glob("*eak.gz")[0]
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
		File reproducibility_qc = 
			glob("*reproducibility.qc")[0]
	}
}

task blacklist_filter {
	# parameters from workflow
	File? peak
	File? blacklist

	command {
		python $(which encode_blacklist_filter.py) \
			${peak} \
			--blacklist ${blacklist}
	}
	output {
		File filtered_peak = glob('*.gz')[0]
	}
}

task frip {
	# parameters from workflow
	File? peak
	File? ta

	command {
		python $(which encode_frip.py) \
			${peak} \
			${ta}
	}
	output {
		File frip_qc = glob('*.frip.qc')[0]
	}
}

# workflow system tasks
# to reduce overhead of provisioning extra nodes
# we have only one task to
# 	1) determine input type and number of replicates	
# 	2) generate pair (rep-x_vs_rep-y) of all true replicate
# 	3) read genome_tsv and get ready to download files 
# 		On Google Cloud Platform 
#		files are downloaded from gs://atac-seq-pipeline-genome-data/
task inputs {
	# parameters from workflow
	Array[Array[Array[String]]] fastqs 
	Array[String] bams
	Array[String] nodup_bams
	Array[String] tas
	Array[String] peaks
	File genome_tsv

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
		String type = read_string("type.txt")
		Int num_rep = read_int("num_rep.txt")

		Boolean is_before_bam =
			type=='fastq'
		Boolean is_before_nodup_bam =
			type=='fastq' || 
			type=='bam'
		Boolean is_before_ta =
			type=='fastq' || 
			type=='bam' ||
			type=='nodup_bam'
		Boolean is_before_peak = 
			type=='fastq' || 
			type=='bam' ||
			type=='nodup_bam' || 
			type=='ta'

		Array[Array[Int]] pairs = if num_rep>1 then 
									read_tsv(stdout()) else [[]]

		String ref_fa = read_map(genome_tsv)['ref_fa']
		String bowtie2_idx_tar = read_map(genome_tsv)['bowtie2_idx_tar']
		String bwa_idx_tar = read_map(genome_tsv)['bwa_idx_tar']
		String blacklist = read_map(genome_tsv)['blacklist']
		String chrsz = read_map(genome_tsv)['chrsz']
		String gensz = read_map(genome_tsv)['gensz']
	}
}
 
task gather_and_report {
	# fastqs
	Array[Array[Array[File]]] fastqs
	# raw bams
	Array[File] bams
	# filtered bams
	Array[File] nodup_bams
	# tag-aligns
	Array[File] tas
	# MACS2 raw peaks
	Array[File] peaks
	Array[File] peaks_pr1
	Array[File] peaks_pr2
	File? peak_ppr1
	File? peak_ppr2
	File? peak_pooled
	# blacklist filtered MACS2 raw peaks
	Array[File] bfilt_peaks
	Array[File] bfilt_peaks_pr1
	Array[File] bfilt_peaks_pr2
	File? bfilt_peak_ppr1
	File? bfilt_peak_ppr2
	File? bfilt_peak_pooled
	# adapters
	Array[Array[Array[String]]] adapters
	# blacklist filtered IDR peaks
	Array[File] bfilt_idr_peaks
	Array[File] bfilt_idr_peaks_pr
	File? bfilt_idr_peak_ppr
	# blacklist filtered overlapping peaks
	Array[File] bfilt_overlap_peaks
	Array[File] bfilt_overlap_peaks_pr
	File? bfilt_overlap_peak_ppr

	File? idr_reproducibility_qc
	File? overlap_reproducibility_qc

	command {
		python $(which encode_gather_atac.py) \
			"--bams " + ${sep=' ' bams} \
			"--nodup-bams " + ${sep=' ' nodup_bams} \
			"--tas " + ${sep=' ' tas} \
			"--peaks " + ${sep=' ' peaks} \
			"--bfilt-peaks " + ${sep=' ' bfilt_peaks}
	}
	output {
		File qc_json = glob('qc.json')[0]
		File report = glob('report.html')[0]
	}
}