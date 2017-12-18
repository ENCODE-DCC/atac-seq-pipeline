# ENCODE DCC ATAC-Seq/DNase-Seq pipeline
# Author: Jin Lee (leepc12@gmail.com)

workflow atac {
	###### pipeline inputs ######

	String pipeline_type  		# atac or dnase

	# input files (choose one of the input types)
	Array[Array[Array[String]]]? fastqs 
								# [rep_id][merge_id][end_id] if starting from fastqs
								# 	after merging, it will reduce to 
								# 	[rep_id][end_id]	
	Array[String]? bams 		# [rep_id] if starting from bams	
	Array[String]? nodup_bams 	# [rep_id] if starting from filtered bams
	Array[String]? tas 			# [rep_id] if starting from tag-aligns

	Array[String]? peaks		# [rep_id] if starting from peaks
	Array[String]? peaks_pr1	# [rep_id] if starting from peaks
	Array[String]? peaks_pr2	# [rep_id] if starting from peaks
	File? peak_ppr1				# if starting from peaks
	File? peak_ppr2				# if starting from peaks
	File? peak_pooled			# if starting from peaks

	# adapters (define if starting from fastqs)
	# if there are no adapters to be trimmed, do not define
	# you can selectively trim adapter for each fastq
	#  by keeping the same structure as "fastqs" and only fill in known adapters
	# activate "atac.trim_adapter.auto_detect_adapter" 
	#  if you want auto adapter detection/removal for non-empty entries in "adapters"
	Array[Array[Array[String]]]? adapters 
								# [rep_id][merge_id][end_id]

	# mandatory genome param
	File genome_tsv 		# reference genome data TSV file including
							# all important genome specific data file paths
							# and parameters
	Boolean paired_end 		# endedness of sample

	# optional but important
	Boolean? align_only 	# disable downstream analysis (peak calling, ...)
							# after alignment
	Boolean? true_rep_only 	# disable all analyses for pseudo replicates
							# naive-overlap and IDR will also be disabled
	Boolean? disable_xcor 	# disable cross-correlation analysis
	Int? multimapping 		# multimapping reads

	# task-specific variables but defined in workflow level (limit of WDL)
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

	###### initialization for pipeline ######

	# temp null variable for optional File/String
	String? null

	# read genome data and paths
	call read_genome_tsv { input: genome_tsv = genome_tsv }	# For Google JES backend
	String bowtie2_idx_tar = read_genome_tsv.genome['bowtie2_idx_tar']
	String? blacklist = if read_genome_tsv.genome['blacklist']=='/dev/null' then null 
					else read_genome_tsv.genome['blacklist']
	String chrsz = read_genome_tsv.genome['chrsz']
	String gensz = read_genome_tsv.genome['gensz']

	# simplified variables for optional flags
	Boolean align_only_ = select_first([align_only, false])
	Boolean true_rep_only_ = select_first([true_rep_only, false])
	Boolean enable_idr_ = select_first([enable_idr, false])
	Boolean disable_xcor_ = select_first([disable_xcor, false])

	###### pipeline starts here ######

	Array[Array[Array[String]]] fastqs_ = select_first([fastqs, []])
	Array[Array[Array[String]]] adapters_ = select_first([adapters, []])
	Int fastqs_len = length(fastqs_)
	Int adapters_len = length(adapters_)
	if ( fastqs_len>0 ) {
		scatter(i in range(fastqs_len)) {
			# trim adapters and merge trimmed fastqs
			call trim_adapter { input :
				fastqs = fastqs_[i],
				adapters = if adapters_len>0 then adapters_[i] else [],
				paired_end = paired_end,
			}
			# align trimmed/merged fastqs with bowtie2
			call bowtie2 { input :
				idx_tar = bowtie2_idx_tar,
				fastqs = trim_adapter.trimmed_merged_fastqs, #[R1,R2]
				paired_end = paired_end,
				multimapping = multimapping,
			}
		}
	}
	Array[String] bams_ = select_first([bams, bowtie2.bam, []])
	if ( length(bams_)>0 ) {
		scatter(bam in bams_) {
			# filter/dedup bam
			call filter { input :
				bam = bam,
				paired_end = paired_end,
				multimapping = multimapping,
			}
		}
	}
	Array[String] nodup_bams_ = select_first([nodup_bams, filter.nodup_bam, []])
	if ( length(nodup_bams_)>0 ) {
		scatter(bam in nodup_bams_) {
			# convert bam to tagalign and subsample it if necessary
			call bam2ta { input :
				bam = bam,
				disable_tn5_shift = if pipeline_type=='atac' then false else true,
				paired_end = paired_end,
			}
		}
	}
	Array[String] tas_ = select_first([tas,bam2ta.ta,[]])
	Int tas_len = length(tas_)
	if ( tas_len>0 ) {
		if ( !disable_xcor_ ) {
			# subsample tagalign (non-mito) and cross-correlation analysis
			scatter(ta in tas_) {
				call xcor { input :
					ta = ta,
					paired_end = paired_end,
				}
			}
		}
		if ( !align_only_ ) {
			# subsample tagalign (non-mito) and cross-correlation analysis
			scatter(ta in tas_) {
				# call peaks on tagalign
				call macs2 { input :
					ta = ta,
					gensz = gensz,
					chrsz = chrsz,
					cap_num_peak = cap_num_peak,
					pval_thresh = pval_thresh,
					smooth_win = smooth_win,
					make_signal = true,
					blacklist = blacklist,
					mem_mb = macs2_mem_mb,
					disks = macs2_disks,
					time_hr = macs2_time_hr,
				}
			}
			if ( !true_rep_only_ ) {
				scatter(ta in tas_) {
					# make two self pseudo replicates per true replicate
					call spr { input :
						ta = ta,
						paired_end = paired_end,
					}
					# call peaks on 1st pseudo replicated tagalign 
					call macs2 as macs2_pr1 { input :
						ta = spr.ta_pr1,
						gensz = gensz,
						chrsz = chrsz,
						cap_num_peak = cap_num_peak,
						pval_thresh = pval_thresh,
						smooth_win = smooth_win,
						blacklist = blacklist,
						mem_mb = macs2_mem_mb,
						disks = macs2_disks,
						time_hr = macs2_time_hr,
					}
					# call peaks on 2nd pseudo replicated tagalign 
					call macs2 as macs2_pr2 { input :
						ta = spr.ta_pr2,
						gensz = gensz,
						chrsz = chrsz,
						cap_num_peak = cap_num_peak,
						pval_thresh = pval_thresh,
						smooth_win = smooth_win,
						blacklist = blacklist,
						mem_mb = macs2_mem_mb,
						disks = macs2_disks,
						time_hr = macs2_time_hr,
					}				
				}
			}
			if ( tas_len>1 ) {
				# pool tagaligns from true replicates
				call pool_ta { input :
					tas = tas_,
				}
				# call peaks on pooled replicate
				call macs2 as macs2_pooled { input :
					ta = pool_ta.ta_pooled,
					gensz = gensz,
					chrsz = chrsz,
					cap_num_peak = cap_num_peak,
					pval_thresh = pval_thresh,
					smooth_win = smooth_win,
					make_signal = true,
					blacklist = blacklist,
					mem_mb = macs2_mem_mb,
					disks = macs2_disks,
					time_hr = macs2_time_hr,
				}
				if ( !true_rep_only_ ) {
					# pool tagaligns from pseudo replicates
					call pool_ta as pool_ta_pr1 { input :
						tas = select_first([spr.ta_pr1]),
					}
					call pool_ta as pool_ta_pr2 { input :
						tas = select_first([spr.ta_pr2]),
					}
					# call peaks on 1st pooled pseudo replicates
					call macs2 as macs2_ppr1 { input :
						ta = pool_ta_pr1.ta_pooled,
						gensz = gensz,
						chrsz = chrsz,
						cap_num_peak = cap_num_peak,
						pval_thresh = pval_thresh,
						smooth_win = smooth_win,
						blacklist = blacklist,
						mem_mb = macs2_mem_mb,
						disks = macs2_disks,
						time_hr = macs2_time_hr,
					}
					# call peaks on 2nd pooled pseudo replicates
					call macs2 as macs2_ppr2 { input :
						ta = pool_ta_pr2.ta_pooled,
						gensz = gensz,
						chrsz = chrsz,
						cap_num_peak = cap_num_peak,
						pval_thresh = pval_thresh,
						smooth_win = smooth_win,
						blacklist = blacklist,
						mem_mb = macs2_mem_mb,
						disks = macs2_disks,
						time_hr = macs2_time_hr,
					}
				}
			}
		}
	}

	Array[String] peaks_ = select_first([peaks, macs2.npeak, []])
	Array[String] peaks_pr1_ = select_first([peaks_pr1, macs2_pr1.npeak, []])
	Array[String] peaks_pr2_ = select_first([peaks_pr2, macs2_pr2.npeak, []])
	Int num_rep = length(peaks_)
	String? peak_pooled_ = select_first([peak_pooled, macs2_pooled.npeak, '/dev/null'])
	String? peak_ppr1_ = select_first([peak_ppr1, macs2_ppr1.npeak, '/dev/null'])
	String? peak_ppr2_ = select_first([peak_ppr2, macs2_ppr2.npeak, '/dev/null'])
	# determine peak_type
	String peak_type = 'narrowPeak'
	# determine idr ranking method
	String idr_rank = 'p.value'
	# generate all possible pairs of true replicates
	call pair_gen { input: num_rep = num_rep }

	if ( !align_only_ ) {
		if ( num_rep>1 ) {
			# Naive overlap on every pair of true replicates
			scatter( pair in pair_gen.pairs ) {
				call overlap { input :
					prefix = "rep"+(pair[0]+1)+"-rep"+(pair[1]+1),
					peak1 = peaks_[(pair[0])],
					peak2 = peaks_[(pair[1])],
					peak_pooled = peak_pooled_,
					peak_type = peak_type,
					blacklist = blacklist,
					ta = if tas_len>0 then pool_ta.ta_pooled else null,
				}
			}
			if ( enable_idr_ ) {
				# IDR on every pair of true replicates
				scatter( pair in pair_gen.pairs ) {
					call idr { input : 
						prefix = "rep"+(pair[0]+1)+"-rep"+(pair[1]+1),
						peak1 = peaks_[(pair[0])],
						peak2 = peaks_[(pair[1])],
						peak_pooled = peak_pooled_,
						idr_thresh = select_first([idr_thresh,0.1]),
						peak_type = peak_type,
						rank = idr_rank,
						blacklist = blacklist,
						ta = if tas_len>0 then pool_ta.ta_pooled else null,
					}
				}
			}
		}
		if ( !true_rep_only_ ) {
			# Naive overlap on pseduo replicates
			scatter( i in range(num_rep) ) {
				call overlap as overlap_pr { input : 
					prefix = "rep"+(i+1)+"-pr",
					peak1 = peaks_pr1_[i],
					peak2 = peaks_pr2_[i],
					peak_pooled = peaks_[i],
					peak_type = peak_type,
					blacklist = blacklist,
					ta = if tas_len>0 then tas_[i] else null,
				}
			}
			if ( enable_idr_ ) {
				# IDR on pseduo replicates
				scatter( i in range(num_rep) ) {
					call idr as idr_pr { input : 
						prefix = "rep"+(i+1)+"-pr",
						peak1 = peaks_pr1_[i],
						peak2 = peaks_pr2_[i],
						peak_pooled = peaks_[i],
						idr_thresh = select_first([idr_thresh,0.1]),
						peak_type = peak_type,
						rank = idr_rank,
						blacklist = blacklist,
						ta = if tas_len>0 then tas_[i] else null,
					}
				}
			}
			if ( num_rep>1 ) {
				# Naive overlap on pooled pseudo replicates
				call overlap as overlap_ppr { input : 
					prefix = "ppr",
					peak1 = peak_ppr1_,
					peak2 = peak_ppr2_,
					peak_pooled = peak_pooled_,
					peak_type = peak_type,
					blacklist = blacklist,
					ta = if tas_len>0 then pool_ta.ta_pooled else null,
				}
				if ( enable_idr_ ) {
					# IDR on pooled pseduo replicates
					call idr as idr_ppr { input : 
						prefix = "ppr",
						peak1 = peak_ppr1_,
						peak2 = peak_ppr2_,
						peak_pooled = peak_pooled_,
						idr_thresh = select_first([idr_thresh,0.1]),
						peak_type = peak_type,
						rank = idr_rank,
						blacklist = blacklist,
						ta = if tas_len>0 then pool_ta.ta_pooled else null,
					}
				}
			}
			# reproducibility QC for overlapping peaks
			call reproducibility as reproducibility_overlap { input :
				prefix = 'overlap',
				peaks = select_first([overlap.bfilt_overlap_peak, []]),
				peaks_pr = overlap_pr.bfilt_overlap_peak,
				peak_ppr = overlap_ppr.bfilt_overlap_peak,
			}
			if ( enable_idr_ ) {
				# reproducibility QC for IDR peaks
				call reproducibility as reproducibility_idr { input :
					prefix = 'idr',
					peaks = select_first([idr.bfilt_idr_peak, []]),
					peaks_pr = idr_pr.bfilt_idr_peak,
					peak_ppr = idr_ppr.bfilt_idr_peak,
				}
			}
		}
	}

	call qc_report { input :
		paired_end = paired_end,
		pipeline_type = pipeline_type,
		peak_caller = 'macs2',
		idr_thresh = select_first([idr_thresh,0.1]),
		flagstat_qcs = select_first([bowtie2.flagstat_qc, []]),
		nodup_flagstat_qcs = select_first([filter.flagstat_qc, []]),
		dup_qcs = select_first([filter.dup_qc, []]),
		pbc_qcs = select_first([filter.pbc_qc, []]),
		xcor_plots = select_first([xcor.plot_png, []]),
		xcor_scores = select_first([xcor.score, []]),

		frip_qcs = select_first([macs2.frip_qc, []]),
		frip_qcs_pr1 = select_first([macs2_pr1.frip_qc, []]),
		frip_qcs_pr2 = select_first([macs2_pr2.frip_qc, []]),
		frip_qc_pooled = if defined(macs2_pooled.frip_qc) then macs2_pooled.frip_qc else null,
		frip_qc_ppr1 = if defined(macs2_ppr1.frip_qc) then macs2_ppr1.frip_qc else null,
		frip_qc_ppr2 = if defined(macs2_ppr2.frip_qc) then macs2_ppr2.frip_qc else null,

		idr_plots = select_first([idr.idr_plot, []]),
		idr_plots_pr = select_first([idr_pr.idr_plot, []]),
		idr_plot_ppr = if defined(idr_ppr.idr_plot) then idr_ppr.idr_plot else null,
		frip_idr_qcs = select_first([idr.frip_qc, []]),
		frip_idr_qcs_pr = select_first([idr_pr.frip_qc, []]),
		frip_idr_qc_ppr = if defined(idr_ppr.frip_qc) then idr_ppr.frip_qc else null,
		frip_overlap_qcs = select_first([overlap.frip_qc, []]),
		frip_overlap_qcs_pr = select_first([overlap_pr.frip_qc, []]),
		frip_overlap_qc_ppr = if defined(overlap_ppr.frip_qc) then overlap_ppr.frip_qc else null,
		idr_reproducibility_qc = if defined(reproducibility_idr.reproducibility_qc) 
								then reproducibility_idr.reproducibility_qc else null,
		overlap_reproducibility_qc = if defined(reproducibility_overlap.reproducibility_qc) 
								then reproducibility_overlap.reproducibility_qc else null,
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
			--adapters ${write_tsv(adapters)} \
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
			${"--p-val-thresh "+ pval_thresh} \
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
								# dup.qc and pbc.qc will be emptry files
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
		touch null # ugly part to deal with optional outputs
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
