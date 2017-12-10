# ENCODE DCC ATAC-Seq/DNase-Seq pipeline for DNANexus
# Author: Jin Lee (leepc12@gmail.com)

workflow atac {
	String pipeline_type  		# atac or dnase
	Array[Array[Array[File]]] fastqs 
								# [rep_id][merge_id][end_id]
								# 	after merging, it will reduce to 
								# 	[rep_id][end_id]
	Array[Array[Array[String]]] adapters 
								# [rep_id][merge_id][end_id]
	Boolean paired_end 		# endedness of sample
	Int? multimapping 		# multimapping reads
	Float? idr_thresh		# IDR threshold
	Int? smooth_win			# smoothing window for MACS2

	# genome ref. data files
	File ref_fa
	File bowtie2_idx_tar
	File blacklist
	File chrsz
	String gensz	

	# init
	Int num_rep = length(fastqs)
	Float idr_thresh_ = select_first([idr_thresh, 0.1])
	Int multimapping_ = select_first([multimapping, 0])
	Int smooth_win_	= select_first([smooth_win, 150])

	# pipeline starts here (parallelized for each replicate)
	scatter(i in range(num_rep)) {
		# trim adapters and merge trimmed fastqs
		call trim_adapter { input :
			fastqs = fastqs[i],
			adapters = if length(adapters)>0 then adapters[i] else [],
			paired_end = paired_end,
		}
		# align trimmed/merged fastqs with bowtie2
		call bowtie2 { input :
			idx_tar = bowtie2_idx_tar,
			fastqs = trim_adapter.trimmed_merged_fastqs, #[R1,R2]
			paired_end = paired_end,
			multimapping = multimapping_,
		}
		# filter/dedup bam
		call filter { input :
			bam = bowtie2.bam,
			paired_end = paired_end,
			multimapping = multimapping_,
		}
		# convert bam to tagalign and subsample it if necessary
		call bam2ta { input :
			bam = filter.nodup_bam,
			disable_tn5_shift = if pipeline_type=='atac' then false else true,
			paired_end = paired_end,
		}
		# subsample tagalign (non-mito) and cross-correlation analysis
		call xcor { input :
			ta = bam2ta.ta,
			paired_end = paired_end,
		}
		# call peaks on tagalign
		call macs2 { input :
			ta = bam2ta.ta,
			gensz = gensz,
			chrsz = chrsz,
			smooth_win = smooth_win_,
			make_signal = true,
			blacklist = blacklist,
		}
		# make two self pseudo replicates per true replicate
		call spr { input :
			ta = bam2ta.ta,
			paired_end = paired_end,
		}
	}

	# pool tagaligns from true replicates
	call pool_ta { input :
		tas = bam2ta.ta,
	}
	# call peaks on pooled replicate
	call macs2 as macs2_pooled { input :
		ta = pool_ta.ta_pooled,
		gensz = gensz,
		chrsz = chrsz,
		smooth_win = smooth_win_,
		make_signal = true,
		blacklist = blacklist,
	}
	scatter(i in range(num_rep)) {
		# call peaks on 1st pseudo replicated tagalign 
		call macs2 as macs2_pr1 { input :
			ta = spr.ta_pr1[i],
			gensz = gensz,
			chrsz = chrsz,
			smooth_win = smooth_win_,
			blacklist = blacklist,
		}
		call macs2 as macs2_pr2 { input :
			ta = spr.ta_pr2[i],
			gensz = gensz,
			chrsz = chrsz,
			smooth_win = smooth_win_,
			blacklist = blacklist,
		}
	}
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
		gensz = gensz,
		chrsz = chrsz,
		smooth_win = smooth_win_,
		blacklist = blacklist,
	}
	# call peaks on 2nd pooled pseudo replicates
	call macs2 as macs2_ppr2 { input :
		ta = pool_ta_pr2.ta_pooled,
		gensz = gensz,
		chrsz = chrsz,
		smooth_win = smooth_win_,
		blacklist = blacklist,
	}
	# generate every pair of true replicates
	call pair_gen { input : 
		num_rep = num_rep,
	}
	# Naive overlap on every pair of true replicates
	scatter( pair in pair_gen.pairs ) {
		call overlap { input :
			prefix = "rep"+(pair[0]+1)+"-rep"+(pair[1]+1),
			peak1 = macs2.npeak[(pair[0])],
			peak2 = macs2.npeak[(pair[1])],
			peak_pooled = macs2_pooled.npeak,
			peak_type = "narrowPeak",
			blacklist = blacklist,
			chrsz = chrsz,
			ta = pool_ta.ta_pooled,
		}
	}
	# Naive overlap on pseduo replicates
	scatter( i in range(num_rep) ) {
		call overlap as overlap_pr { input : 
			prefix = "rep"+(i+1)+"-pr",
			peak1 = macs2_pr1.npeak[i],
			peak2 = macs2_pr2.npeak[i],
			peak_pooled = macs2.npeak[i],
			peak_type = "narrowPeak",
			blacklist = blacklist,
			chrsz = chrsz,
			ta = bam2ta.ta[i],
		}
	}
	# Naive overlap on pooled pseudo replicates
	call overlap as overlap_ppr { input : 
		prefix = "ppr",
		peak1 = macs2_ppr1.npeak,
		peak2 = macs2_ppr2.npeak,
		peak_pooled = macs2_pooled.npeak,
		peak_type = "narrowPeak",
		blacklist = blacklist,
		chrsz = chrsz,
		ta = pool_ta.ta_pooled,
	}
	# reproducibility QC for overlapping peaks
	call reproducibility as reproducibility_overlap { input :
		prefix = 'overlap',
		peaks = overlap.bfilt_overlap_peak,
		peaks_pr = overlap_pr.bfilt_overlap_peak,
		peak_ppr = overlap_ppr.bfilt_overlap_peak,
	}

	# IDR on every pair of true replicates
	scatter( pair in pair_gen.pairs ) {
		call idr { input : 
			prefix = "rep"+(pair[0]+1)+"-rep"+(pair[1]+1),
			peak1 = macs2.npeak[(pair[0])],
			peak2 = macs2.npeak[(pair[1])],
			peak_pooled = macs2_pooled.npeak,
			idr_thresh = idr_thresh_,
			peak_type = "narrowPeak",
			rank = "p.value",
			blacklist = blacklist,
			chrsz = chrsz,
			ta = pool_ta.ta_pooled,
		}
	}
	# IDR on pseduo replicates
	scatter( i in range(num_rep) ) {
		call idr as idr_pr { input : 
			prefix = "rep"+(i+1)+"-pr",
			peak1 = macs2_pr1.npeak[i],
			peak2 = macs2_pr2.npeak[i],
			peak_pooled = macs2.npeak[i],
			idr_thresh = idr_thresh_,
			peak_type = "narrowPeak",
			rank = "p.value",
			blacklist = blacklist,
			chrsz = chrsz,
			ta = bam2ta.ta[i],
		}
	}
	# IDR on pooled pseduo replicates
	call idr as idr_ppr { input : 
		prefix = "ppr",
		peak1 = macs2_ppr1.npeak,
		peak2 = macs2_ppr2.npeak,
		peak_pooled = macs2_pooled.npeak,
		idr_thresh = idr_thresh_,
		peak_type = "narrowPeak",
		rank = "p.value",
		blacklist = blacklist,
		chrsz = chrsz,
		ta = pool_ta.ta_pooled,
	}
	# reproducibility QC for IDR peaks
	call reproducibility as reproducibility_idr { input :
		prefix = 'idr',
		peaks = idr.bfilt_idr_peak,
		peaks_pr = idr_pr.bfilt_idr_peak,
		peak_ppr = idr_ppr.bfilt_idr_peak,
	}

	call qc_report { input :
		paired_end = paired_end,
		pipeline_type = pipeline_type,
		peak_caller = 'macs2',
		idr_thresh = idr_thresh_,
		flagstat_qcs = bowtie2.flagstat_qc,
		nodup_flagstat_qcs = filter.flagstat_qc,
		dup_qcs = filter.dup_qc,
		pbc_qcs = filter.pbc_qc,
		xcor_plots = xcor.plot_png,
		xcor_scores = xcor.score,
		frip_qcs = macs2.frip_qc,
		frip_qcs_pr1 = macs2_pr1.frip_qc,
		frip_qcs_pr2 = macs2_pr2.frip_qc,
		frip_qc_pooled = [macs2_pooled.frip_qc],
		frip_qc_ppr1 = [macs2_ppr1.frip_qc],
		frip_qc_ppr2 = [macs2_ppr2.frip_qc],

		idr_plots = idr.idr_plot,
		idr_plots_pr = idr_pr.idr_plot,
		idr_plot_ppr = [idr_ppr.idr_plot],
		frip_idr_qcs = idr.frip_qc,
		frip_idr_qcs_pr = idr_pr.frip_qc,
		frip_idr_qc_ppr = [idr_ppr.frip_qc],
		frip_overlap_qcs = overlap.frip_qc,
		frip_overlap_qcs_pr = overlap_pr.frip_qc,
		frip_overlap_qc_ppr = [overlap_ppr.frip_qc],
		idr_reproducibility_qc = [reproducibility_idr.reproducibility_qc],
		overlap_reproducibility_qc = [reproducibility_overlap.reproducibility_qc],
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
	File blacklist 		# blacklist BED to filter raw peaks
	# fixed var
	String peak_type = "narrowPeak"
	# resource
	Int? mem_mb
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
		touch null 
	}
	output {
		File npeak = glob("*[!.][!b][!f][!i][!l][!t]."+peak_type+".gz")[0]
		File bfilt_npeak = glob("*.bfilt."+peak_type+".gz")[0]
		File sig_pval = if select_first([make_signal,false]) then glob("*.pval.signal.bigwig")[0] else glob("null")[0]
		File sig_fc = if select_first([make_signal,false]) then glob("*.fc.signal.bigwig")[0] else glob("null")[0]
		File frip_qc = glob("*.frip.qc")[0]
	}
	runtime {
		memory : "${select_first([mem_mb,'16000'])} MB"
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
	String? disks
	# temporary boolean (cromwell bug when if select_first is in output)
	Boolean no_dup_removal_ = select_first([no_dup_removal,false])

	command {
		python $(which encode_filter.py) \
			${bam} \
			${if paired_end then "--paired-end" else ""} \
			${"--multimapping " + multimapping} \
			${"--dup-marker " + dup_marker} \
			${"--mapq-thresh " + mapq_thresh} \
			${if no_dup_removal_ then "--no-dup-removal" else ""} \
			${"--nth " + cpu}

		# ugly part to deal with optional outputs
		${if no_dup_removal_ then "touch null.pbc.qc null.dup.qc" else ""}
		touch null
	}
	output {
		File nodup_bam = glob("*.bam")[0]
		File nodup_bai = glob("*.bai")[0]
		File flagstat_qc = glob("*.flagstat.qc")[0]
		File dup_qc = if no_dup_removal_ then glob("null")[0] else glob("*.dup.qc")[0]
		File pbc_qc = if no_dup_removal_ then glob("null")[0] else glob("*.pbc.qc")[0]
	}

	runtime {
		cpu : select_first([cpu,2])
		memory : "${select_first([mem_mb,'20000'])} MB"
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
	File peak1 			
	File peak2
	File peak_pooled
	Float? idr_thresh
	File blacklist 	# blacklist BED to filter raw peaks
	# parameters to compute FRiP
	File ta		# to calculate FRiP
	Int? fraglen 		# fragment length from xcor
	File chrsz			# 2-col chromosome sizes file
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
			${"--ta "+ ta}
	}
	output {
		File idr_peak = glob("*[!.][!b][!f][!i][!l][!t]."+peak_type+".gz")[0]
		File bfilt_idr_peak = glob("*.bfilt."+peak_type+".gz")[0]
		File idr_plot = glob("*.txt.png")[0]
		File idr_unthresholded_peak = glob("*.txt.gz")[0]
		File idr_log = glob("*.log")[0]
		File frip_qc = glob("*.frip.qc")[0]
	}
}

task overlap {
	# parameters from workflow
	String prefix 		# prefix for IDR output file
	File peak1
	File peak2
	File peak_pooled
	File blacklist 	# blacklist BED to filter raw peaks
	# parameters to compute FRiP
	File ta		# to calculate FRiP
	Int? fraglen 		# fragment length from xcor
	File chrsz			# 2-col chromosome sizes file
	String peak_type

	command {
		python $(which encode_naive_overlap.py) \
			${peak1} ${peak2} ${peak_pooled} \
			${"--prefix " + prefix} \
			${"--peak-type " + peak_type} \
			${"--fraglen " + fraglen} \
			${"--chrsz " + chrsz} \
			${"--blacklist "+ blacklist} \
			${"--ta "+ ta}
	}
	output {
		File overlap_peak = glob("*[!.][!b][!f][!i][!l][!t]."+peak_type+".gz")[0]
		File bfilt_overlap_peak = glob("*.bfilt."+peak_type+".gz")[0]
		File frip_qc = glob("*.frip.qc")[0]
	}
}

task reproducibility {
	# parameters from workflow
	String prefix
	Array[File] peaks # peak files from pair of true replicates
						# in a sorted order. for example of 4 replicates,
						# 1,2 1,3 1,4 2,3 2,4 3,4.
                        # x,y means peak file from rep-x vs rep-y
	Array[File] peaks_pr	# peak files from pseudo replicates
	File peak_ppr			# Peak file from pooled pseudo replicate.

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
	Array[File] flagstat_qcs
	Array[File] nodup_flagstat_qcs
	Array[File] dup_qcs
	Array[File] pbc_qcs
	Array[File] xcor_plots
	Array[File] xcor_scores
	Array[File] idr_plots
	Array[File] idr_plots_pr
	Array[File] idr_plot_ppr # not actually an array
	Array[File] frip_qcs
	Array[File] frip_qcs_pr1
	Array[File] frip_qcs_pr2
	Array[File] frip_qc_pooled # not actually an array
	Array[File] frip_qc_ppr1 # not actually an array
	Array[File] frip_qc_ppr2 # not actually an array
	Array[File] frip_idr_qcs
	Array[File] frip_idr_qcs_pr
	Array[File] frip_idr_qc_ppr # not actually an array
	Array[File] frip_overlap_qcs
	Array[File] frip_overlap_qcs_pr
	Array[File] frip_overlap_qc_ppr # not actually an array
	Array[File] idr_reproducibility_qc # not actually an array
	Array[File] overlap_reproducibility_qc # not actually an array

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
			--idr-plot-ppr ${sep=' ' idr_plot_ppr} \
			--frip-qcs ${sep=' ' frip_qcs} \
			--frip-qcs-pr1 ${sep=' ' frip_qcs_pr1} \
			--frip-qcs-pr2 ${sep=' ' frip_qcs_pr2} \
			--frip-qc-pooled ${sep=' ' frip_qc_pooled} \
			--frip-qc-ppr1 ${sep=' ' frip_qc_ppr1} \
			--frip-qc-ppr2 ${sep=' ' frip_qc_ppr2} \
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
	output {
		File report = glob('*qc.html')[0]
		File qc_json = glob('*qc.json')[0]
		#File encode_accession_json= glob('*encode_accession.json')[0]
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
		Array[Array[Int]] pairs = read_tsv(stdout())
	}
}
