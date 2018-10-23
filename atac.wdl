# ENCODE DCC ATAC-Seq/DNase-Seq pipeline
# Author: Jin Lee (leepc12@gmail.com)

workflow atac {
	#### input file definition
		# pipeline can start from any type of inputs and then leave all other types undefined
		# supported types: fastq, bam, nodup_bam (filtered bam), ta (tagAlign), peak
		# define up to 6 replicates
		# [rep_id] is for each replicate

 	### fastqs and adapters  	
	 	# define fastqs either with DNANexus style (1-dim array) or with default one (3-dim array)
	 	# [merge_id] is for pooing fastqs after trimming adapters
	 	# if adapters defined with any style, keep the same structure/dimension as fastq arrays
	 	# only defined adapters will be trimmed
	 	# or undefined adapters will be detected/trimmed by trim_adapter.auto_detect_adapter=true 
	 	# so you can selectively detect/trim adapters for a specific fastq
 	## DNANexus UI style fastq/adapter definition
	Array[File] fastqs_rep1_R1 = []	# [merge_id]
	Array[File] fastqs_rep1_R2 = [] # do not define _R2 array if your sample is not paired end
	Array[File] fastqs_rep2_R1 = [] # do not define if you have a single replicate
	Array[File] fastqs_rep2_R2 = []	# do not define _R2 array if your sample is not paired end
	Array[File] fastqs_rep3_R1 = [] # do not define if you have <=2 replicates
	Array[File] fastqs_rep3_R2 = []	# do not define _R2 array if your sample is not paired end
	Array[File] fastqs_rep4_R1 = [] # do not define if you have <=3 replicates
	Array[File] fastqs_rep4_R2 = []	# do not define _R2 array if your sample is not paired end
	Array[File] fastqs_rep5_R1 = [] # do not define if you have <=4 replicates
	Array[File] fastqs_rep5_R2 = []	# do not define _R2 array if your sample is not paired end
	Array[File] fastqs_rep6_R1 = [] # do not define if you have <=5 replicates
	Array[File] fastqs_rep6_R2 = []	# do not define _R2 array if your sample is not paired end
	Array[String] adapters_rep1_R1 = [] # [merge_id]
	Array[String] adapters_rep1_R2 = [] 
	Array[String] adapters_rep2_R1 = []
	Array[String] adapters_rep2_R2 = []
	Array[String] adapters_rep3_R1 = []
	Array[String] adapters_rep3_R2 = []
	Array[String] adapters_rep4_R1 = []
	Array[String] adapters_rep4_R2 = []
	Array[String] adapters_rep5_R1 = []
	Array[String] adapters_rep5_R2 = []
	Array[String] adapters_rep6_R1 = []
	Array[String] adapters_rep6_R2 = []
 	## default style fastq/adapter definition
 		# [read_end_id] is for fastq R1 or fastq R2
	Array[Array[Array[File]]] fastqs = [] 	# [rep_id][merge_id][read_end_id]
	Array[Array[Array[File]]] adapters = []	# [rep_id][merge_id][read_end_id]

	### other input types (bam, nodup_bam, ta)
	Array[File] bams = [] 		# [rep_id]
	Array[File] nodup_bams = [] # [rep_id]
	Array[File] tas = []		# [rep_id]

	### other input types (peak)
	Array[File] peaks = []		# [PAIR(rep_id1,rep_id2)]. example for 3 reps: [rep1_rep2, rep1_rep3, rep2_rep3]
	Array[File] peaks_pr1 = []	# [rep_id]. do not define if true_rep=true
	Array[File] peaks_pr2 = []	# [rep_id]. do not define if true_rep=true
	File? peak_ppr1				# do not define if you have a single replicate or true_rep=true
	File? peak_ppr2				# do not define if you have a single replicate or true_rep=true
	File? peak_pooled			# do not define if you have a single replicate or true_rep=true

	### pipeline type
	String pipeline_type  	# ATAC-Seq (atac) or DNase-Seq (dnase)
							# the only difference is that tn5 shiting is enabled for atac

	### mandatory genome param
	File genome_tsv 		# reference genome data TSV file including
							# all important genome specific data file paths and parameters
	Boolean paired_end

	### optional but important
	Boolean align_only = false		# disable all post-align analysis (peak-calling, overlap, idr, ...)
	Boolean true_rep_only = false 	# disable all analyses for pseudo replicates
									# overlap and idr will also be disabled
	Boolean disable_xcor = false 	# disable cross-correlation analysis
	Int multimapping = 0			# for multimapping reads
	Boolean disable_ataqc = false

	### task-specific variables but defined in workflow level (limit of WDL)
	## optional for MACS2 
	Int cap_num_peak = 300000	# cap number of raw peaks called from MACS2
	Float pval_thresh = 0.01	# p.value threshold
	Int smooth_win = 150		# size of smoothing window
	Int? macs2_mem_mb 			# resource (memory in MB)
	Int? macs2_time_hr			# resource (walltime in hour)
	String? macs2_disks 		# resource disks for cloud platforms
	## optional for IDR
	Boolean enable_idr = false 	# enable IDR analysis on raw peaks
	Float idr_thresh = 0.1		# IDR threshold

	### temp vars (do not define these)
	String peak_type = 'narrowPeak' # peak type for IDR and overlap
	String idr_rank = 'p.value' # IDR ranking method

	### read genome data and paths
	call read_genome_tsv { input:genome_tsv = genome_tsv }
	File bowtie2_idx_tar = read_genome_tsv.genome['bowtie2_idx_tar']
	File blacklist = read_genome_tsv.genome['blacklist']
	File chrsz = read_genome_tsv.genome['chrsz']
	String gensz = read_genome_tsv.genome['gensz']
	File ref_fa = read_genome_tsv.genome['ref_fa']
	# genome data for ATAQC
	File tss_enrich = read_genome_tsv.genome['tss_enrich']
	File dnase = read_genome_tsv.genome['dnase']
	File prom = read_genome_tsv.genome['prom']
	File enh = read_genome_tsv.genome['enh']
	File reg2map = read_genome_tsv.genome['reg2map']
	File reg2map_bed = read_genome_tsv.genome['reg2map_bed']
	File roadmap_meta = read_genome_tsv.genome['roadmap_meta']

	### pipeline starts here
	# temporary 2-dim arrays for DNANexus style fastqs and adapters	
	Array[Array[File]] fastqs_rep1 = if length(fastqs_rep1_R2)>0 then transpose([fastqs_rep1_R1,fastqs_rep1_R2])
									else transpose([fastqs_rep1_R1])
	Array[Array[File]] fastqs_rep2 = if length(fastqs_rep2_R2)>0 then transpose([fastqs_rep2_R1,fastqs_rep2_R2])
									else transpose([fastqs_rep2_R1])
	Array[Array[File]] fastqs_rep3 = if length(fastqs_rep3_R2)>0 then transpose([fastqs_rep3_R1,fastqs_rep3_R2])
									else transpose([fastqs_rep3_R1])
	Array[Array[File]] fastqs_rep4 = if length(fastqs_rep4_R2)>0 then transpose([fastqs_rep4_R1,fastqs_rep4_R2])
									else transpose([fastqs_rep4_R1])
	Array[Array[File]] fastqs_rep5 = if length(fastqs_rep5_R2)>0 then transpose([fastqs_rep5_R1,fastqs_rep5_R2])
									else transpose([fastqs_rep5_R1])
	Array[Array[File]] fastqs_rep6 = if length(fastqs_rep6_R2)>0 then transpose([fastqs_rep6_R1,fastqs_rep6_R2])
									else transpose([fastqs_rep6_R1])
	Array[Array[String]] adapters_rep1 = if length(adapters_rep1_R2)>0 then transpose([adapters_rep1_R1,adapters_rep1_R2])
									else transpose([adapters_rep1_R1])
	Array[Array[String]] adapters_rep2 = if length(adapters_rep2_R2)>0 then transpose([adapters_rep2_R1,adapters_rep2_R2])
									else transpose([adapters_rep2_R1])
	Array[Array[String]] adapters_rep3 = if length(adapters_rep3_R2)>0 then transpose([adapters_rep3_R1,adapters_rep3_R2])
									else transpose([adapters_rep3_R1])
	Array[Array[String]] adapters_rep4 = if length(adapters_rep4_R2)>0 then transpose([adapters_rep4_R1,adapters_rep4_R2])
									else transpose([adapters_rep4_R1])
	Array[Array[String]] adapters_rep5 = if length(adapters_rep5_R2)>0 then transpose([adapters_rep5_R1,adapters_rep5_R2])
									else transpose([adapters_rep5_R1])
	Array[Array[String]] adapters_rep6 = if length(adapters_rep6_R2)>0 then transpose([adapters_rep6_R1,adapters_rep6_R2])
									else transpose([adapters_rep6_R1])

	Array[Array[Array[File]]] fastqs_ = if length(fastqs_rep1)<1 then fastqs
		else if length(fastqs_rep2)<1 then [fastqs_rep1]
		else if length(fastqs_rep3)<1 then [fastqs_rep1,fastqs_rep2]
		else if length(fastqs_rep4)<1 then [fastqs_rep1,fastqs_rep2,fastqs_rep3]
		else if length(fastqs_rep5)<1 then [fastqs_rep1,fastqs_rep2,fastqs_rep3,fastqs_rep4]
		else if length(fastqs_rep6)<1 then [fastqs_rep1,fastqs_rep2,fastqs_rep3,fastqs_rep4,fastqs_rep5]
		else [fastqs_rep1,fastqs_rep2,fastqs_rep3,fastqs_rep4,fastqs_rep5,fastqs_rep6]
	Array[Array[Array[String]]] adapters_ = if length(adapters_rep1)<1 then adapters
		else if length(adapters_rep2)<1 then [adapters_rep1]
		else if length(adapters_rep3)<1 then [adapters_rep1,adapters_rep2]
		else if length(adapters_rep4)<1 then [adapters_rep1,adapters_rep2,adapters_rep3]
		else if length(adapters_rep5)<1 then [adapters_rep1,adapters_rep2,adapters_rep3,adapters_rep4]
		else if length(adapters_rep6)<1 then [adapters_rep1,adapters_rep2,adapters_rep3,adapters_rep4,adapters_rep5]
		else [adapters_rep1,adapters_rep2,adapters_rep3,adapters_rep4,adapters_rep5,adapters_rep6]

	scatter( i in range(length(fastqs_)) ) {
		# trim adapters and merge trimmed fastqs
		call trim_adapter { input :
			fastqs = fastqs_[i],
			adapters = if length(adapters_)>0 then adapters_[i] else [],
			paired_end = paired_end,
		}
		# align trimmed/merged fastqs with bowtie2s
		call bowtie2 { input :
			idx_tar = bowtie2_idx_tar,
			fastqs = trim_adapter.trimmed_merged_fastqs, #[R1,R2]
			paired_end = paired_end,
			multimapping = multimapping,
		}
	}

	Array[File] bams_ = flatten([bowtie2.bam, bams])
	scatter( bam in bams_ ) {
		# filter/dedup bam
		call filter { input :
			bam = bam,
			paired_end = paired_end,
			multimapping = multimapping,
		}
	}

	Array[File] nodup_bams_ = flatten([filter.nodup_bam, nodup_bams])
	scatter( bam in nodup_bams_ ) {
		# convert bam to tagalign and subsample it if necessary
		call bam2ta { input :
			bam = bam,
			disable_tn5_shift = if pipeline_type=='atac' then false else true,
			paired_end = paired_end,
		}
	}

	Array[File] tas_ = if align_only then [] else flatten([bam2ta.ta, tas])
	scatter( ta in tas_ ) {
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
			# resource
			mem_mb = macs2_mem_mb,
			disks = macs2_disks,
			time_hr = macs2_time_hr,
		}
	}
	if ( length(tas_)>1 ) {
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
			# resource
			mem_mb = macs2_mem_mb,
			disks = macs2_disks,
			time_hr = macs2_time_hr,
		}
	}
	if ( !disable_xcor ) {
		scatter( ta in tas_ ) {
			# subsample tagalign (non-mito) and cross-correlation analysis
			call xcor { input :
				ta = ta,
				paired_end = paired_end,
			}
		}
	}

	if ( !true_rep_only ) {
		scatter( ta in tas_ ) {
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
				# resource
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
				# resource
				mem_mb = macs2_mem_mb,
				disks = macs2_disks,
				time_hr = macs2_time_hr,
			}
		}
	}

	if ( !true_rep_only && length(tas_)>1 ) {
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
			cap_num_peak = cap_num_peak,
			pval_thresh = pval_thresh,
			smooth_win = smooth_win,
			blacklist = blacklist,
			# resource
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
			# resource
			mem_mb = macs2_mem_mb,
			disks = macs2_disks,
			time_hr = macs2_time_hr,
		}
	}

	# make peak arrays
	Array[File] peaks_ = flatten([macs2.npeak, peaks])

	# generate all possible pairs of true replicates (pair: left=prefix, right=[peak1,peak2])
	Array[Pair[String,Array[File]]] peak_pairs =  
		if length(peaks_)<=1 then [] # 1 rep
		else if length(peaks_)<=2 then # 2 reps
			 [('rep1-rep2',[peaks_[0],peaks_[1]])]
		else if length(peaks_)<=3 then # 3 reps
			 [('rep1-rep2',[peaks_[0],peaks_[1]]), ('rep1-rep3',[peaks_[0],peaks_[2]]),
			  ('rep2-rep3',[peaks_[1],peaks_[2]])]
		else if length(peaks_)<=4 then # 4 reps
			 [('rep1-rep2',[peaks_[0],peaks_[1]]), ('rep1-rep3',[peaks_[0],peaks_[2]]), ('rep1-rep4',[peaks_[0],peaks_[3]]),
			  ('rep2-rep3',[peaks_[1],peaks_[2]]), ('rep2-rep4',[peaks_[1],peaks_[3]]),
			  ('rep3-rep4',[peaks_[2],peaks_[3]])]
		else if length(peaks_)<=5 then # 5 reps
			 [('rep1-rep2',[peaks_[0],peaks_[1]]), ('rep1-rep3',[peaks_[0],peaks_[2]]), ('rep1-rep4',[peaks_[0],peaks_[3]]), ('rep1-rep5',[peaks_[0],peaks_[4]]),
			  ('rep2-rep3',[peaks_[1],peaks_[2]]), ('rep2-rep4',[peaks_[1],peaks_[3]]), ('rep2-rep5',[peaks_[1],peaks_[4]]),
			  ('rep3-rep4',[peaks_[2],peaks_[3]]), ('rep3-rep5',[peaks_[2],peaks_[4]]),
			  ('rep4-rep5',[peaks_[3],peaks_[4]])]
		else # 6 reps
			 [('rep1-rep2',[peaks_[0],peaks_[1]]), ('rep1-rep3',[peaks_[0],peaks_[2]]), ('rep1-rep4',[peaks_[0],peaks_[3]]), ('rep1-rep5',[peaks_[0],peaks_[4]]), ('rep1-rep6',[peaks_[0],peaks_[5]]),
			  ('rep2-rep3',[peaks_[1],peaks_[2]]), ('rep2-rep4',[peaks_[1],peaks_[3]]), ('rep2-rep5',[peaks_[1],peaks_[4]]), ('rep2-rep6',[peaks_[1],peaks_[5]]),
			  ('rep3-rep4',[peaks_[2],peaks_[3]]), ('rep3-rep5',[peaks_[2],peaks_[4]]), ('rep3-rep6',[peaks_[2],peaks_[5]]),
			  ('rep4-rep5',[peaks_[3],peaks_[4]]), ('rep4-rep6',[peaks_[3],peaks_[5]]),
			  ('rep5-rep6',[peaks_[4],peaks_[5]])]
	scatter( pair in peak_pairs ) {
		# Naive overlap on every pair of true replicates
		call overlap { input :
			prefix = pair.left,
			peak1 = pair.right[0],
			peak2 = pair.right[1],
			peak_pooled = select_first([macs2_pooled.npeak, peak_pooled]),
			peak_type = peak_type,
			blacklist = blacklist,
			chrsz = chrsz,
			ta = pool_ta.ta_pooled,
		}
	}
	if ( enable_idr ) {
		scatter( pair in peak_pairs ) {
			# IDR on every pair of true replicates
			call idr { input : 
				prefix = pair.left,
				peak1 = pair.right[0],
				peak2 = pair.right[1],
				peak_pooled = select_first([macs2_pooled.npeak, peak_pooled]),
				idr_thresh = idr_thresh,
				peak_type = peak_type,
				rank = idr_rank,
				blacklist = blacklist,
				chrsz = chrsz,
				ta = pool_ta.ta_pooled,
			}
		}
	}

	Array[File] peaks_pr1_ = select_first([macs2_pr1.npeak, peaks_pr1])
	Array[File] peaks_pr2_ = select_first([macs2_pr2.npeak, peaks_pr2])

	scatter( i in range(length(peaks_pr1_)) ) {
		# Naive overlap on pseduo replicates
		call overlap as overlap_pr { input : 
			prefix = "rep"+(i+1)+"-pr",
			peak1 = peaks_pr1_[i],
			peak2 = peaks_pr2_[i],
			peak_pooled = peaks_[i],
			peak_type = peak_type,
			blacklist = blacklist,
			chrsz = chrsz,
			ta = if length(tas_)>0 then tas_[i] else pool_ta.ta_pooled,
		}
	}
	if ( enable_idr ) {
		scatter( i in range(length(peaks_pr1_)) ) {
			# IDR on pseduo replicates
			call idr as idr_pr { input : 
				prefix = "rep"+(i+1)+"-pr",
				peak1 = peaks_pr1_[i],
				peak2 = peaks_pr2_[i],
				peak_pooled = peaks_[i],
				idr_thresh = idr_thresh,
				peak_type = peak_type,
				rank = idr_rank,
				blacklist = blacklist,
				chrsz = chrsz,
				ta = if length(tas_)>0 then tas_[i] else pool_ta.ta_pooled,
			}
		}
	}
	if ( length(peaks_pr1_)>1 ) {
		# Naive overlap on pooled pseudo replicates
		call overlap as overlap_ppr { input : 
			prefix = "ppr",
			peak1 = select_first([macs2_ppr1.npeak, peak_ppr1]),
			peak2 = select_first([macs2_ppr2.npeak, peak_ppr2]),
			peak_pooled = select_first([macs2_pooled.npeak, peak_pooled]),
			peak_type = peak_type,
			blacklist = blacklist,
			chrsz = chrsz,
			ta = pool_ta.ta_pooled,
		}
	}
	if ( enable_idr && length(peaks_pr1_)>1  ) {
		# IDR on pooled pseduo replicates
		call idr as idr_ppr { input : 
			prefix = "ppr",
			peak1 = select_first([macs2_ppr1.npeak, peak_ppr1]),
			peak2 = select_first([macs2_ppr2.npeak, peak_ppr2]),
			peak_pooled = select_first([macs2_pooled.npeak, peak_pooled]),
			idr_thresh = idr_thresh,
			peak_type = peak_type,
			rank = idr_rank,
			blacklist = blacklist,
			chrsz = chrsz,
			ta = pool_ta.ta_pooled,
		}
	}
	if ( !align_only && !true_rep_only ) {
		# reproducibility QC for overlapping peaks
		call reproducibility as reproducibility_overlap { input :
			prefix = 'overlap',
			peaks = overlap.bfilt_overlap_peak,
			peaks_pr = overlap_pr.bfilt_overlap_peak,
			peak_ppr = overlap_ppr.bfilt_overlap_peak,
			peak_type = peak_type,
			chrsz = chrsz,
		}
	}
	if ( !align_only && !true_rep_only && enable_idr ) {
		# reproducibility QC for IDR peaks
		call reproducibility as reproducibility_idr { input :
			prefix = 'idr',
			peaks = idr.bfilt_idr_peak,
			peaks_pr = idr_pr.bfilt_idr_peak,
			peak_ppr = idr_ppr.bfilt_idr_peak,
			peak_type = peak_type,
			chrsz = chrsz,
		}
	}
	# Generate final QC report and JSON		
	call qc_report { input :
		paired_end = paired_end,
		pipeline_type = pipeline_type,
		peak_caller = 'macs2',
		idr_thresh = idr_thresh,
		flagstat_qcs = bowtie2.flagstat_qc,
		nodup_flagstat_qcs = filter.flagstat_qc,
		dup_qcs = filter.dup_qc,
		pbc_qcs = filter.pbc_qc,
		xcor_plots = xcor.plot_png,
		xcor_scores = xcor.score,

		frip_macs2_qcs = macs2.frip_qc,
		frip_macs2_qcs_pr1 = macs2_pr1.frip_qc,
		frip_macs2_qcs_pr2 = macs2_pr2.frip_qc,
		frip_macs2_qc_pooled = macs2_pooled.frip_qc,
		frip_macs2_qc_ppr1 = macs2_ppr1.frip_qc,
		frip_macs2_qc_ppr2 = macs2_ppr2.frip_qc,

		idr_plots = idr.idr_plot,
		idr_plots_pr = idr_pr.idr_plot,
		idr_plot_ppr = idr_ppr.idr_plot,
		frip_idr_qcs = idr.frip_qc,
		frip_idr_qcs_pr = idr_pr.frip_qc,
		frip_idr_qc_ppr = idr_ppr.frip_qc,
		frip_overlap_qcs = overlap.frip_qc,
		frip_overlap_qcs_pr = overlap_pr.frip_qc,
		frip_overlap_qc_ppr = overlap_ppr.frip_qc,
		idr_reproducibility_qc = reproducibility_idr.reproducibility_qc,
		overlap_reproducibility_qc = reproducibility_overlap.reproducibility_qc,
	}

	# ATAQC is available only when pipeline starts from fastqs, take fastqs[] as base array for ataqc
	Array[Array[Array[File]]] fastqs_ataqc = 
		if basename(tss_enrich)=='null' || disable_ataqc || align_only || true_rep_only then [] else fastqs_

	scatter( i in range(length(fastqs_ataqc)) ) {
		call ataqc { input : 
			paired_end = paired_end,
			read_len_log = bowtie2.read_len_log[i],
			flagstat_log = bowtie2.flagstat_qc[i],
			bowtie2_log = bowtie2.align_log[i],
			pbc_log = filter.pbc_qc[i],
			dup_log = filter.dup_qc[i],
			bam = bams_[i],
			nodup_flagstat_log = filter.flagstat_qc[i],
			mito_dup_log = filter.mito_dup_log[i],
			nodup_bam = nodup_bams_[i],
			ta = tas_[i],
			peak = if enable_idr then select_first([idr_pr.bfilt_idr_peak])[i]
					else reproducibility_overlap.optimal_peak,
			idr_peak = reproducibility_idr.optimal_peak, #idr_peaks_ataqc[i],
			overlap_peak= reproducibility_overlap.optimal_peak, #overlap_peaks_ataqc[i],
			bigwig = macs2.sig_pval[i],
			ref_fa = ref_fa,
			chrsz = chrsz,
			tss_enrich = tss_enrich,
			blacklist = blacklist,
			dnase = dnase,
			prom = prom,
			enh = enh,
			reg2map_bed = reg2map_bed,
			reg2map = reg2map,
			roadmap_meta = roadmap_meta,
		}
	}

	output {
		File report = qc_report.report
		File qc_json = qc_report.qc_json
		Boolean qc_json_match = qc_report.qc_json_match
	}
}

task trim_adapter { # trim adapters and merge trimmed fastqs
	# parameters from workflow
	Array[Array[File]] fastqs 		# [merge_id][read_end_id]
	Array[Array[String]] adapters 	# [merge_id][read_end_id]
	Boolean paired_end
	# mandatory
	Boolean? auto_detect_adapter	# automatically detect/trim adapters
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
			${if select_first([auto_detect_adapter,false]) then "--auto-detect-adapter" else ""} \
			${"--min-trim-len " + select_first([min_trim_len,5])} \
			${"--err-rate " + select_first([err_rate,'0.1'])} \
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
		memory : "${select_first([mem_mb,'12000'])} MB"
		time : select_first([time_hr,24])
		disks : select_first([disks,"local-disk 100 HDD"])
	}
}

task bowtie2 {
	# parameters from workflow
	File idx_tar 		# reference bowtie2 index tar
	Array[File] fastqs 	# [read_end_id]
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
		File read_len_log = glob("*.read_length.txt")[0] # read_len
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
	String? dup_marker 				# picard.jar MarkDuplicates (picard) or 
									# sambamba markdup (sambamba)
	Int? mapq_thresh				# threshold for low MAPQ reads removal
	Boolean? no_dup_removal 		# no dupe reads removal when filtering BAM
									# dup.qc and pbc.qc will be empty files
									# and nodup_bam in the output is 
	# resource						# filtered bam with dupes	
	Int? cpu
	Int? mem_mb
	Int? time_hr
	String? disks

	command {
		python $(which encode_filter.py) \
			${bam} \
			${if paired_end then "--paired-end" else ""} \
			${"--multimapping " + multimapping} \
			${"--dup-marker " + select_first([dup_marker,'picard'])} \
			${"--mapq-thresh " + select_first([mapq_thresh,30])} \
			${if select_first([no_dup_removal,false]) then "--no-dup-removal" else ""} \
			${"--nth " + cpu}
		# ugly part to deal with optional outputs with Google JES backend
		${if select_first([no_dup_removal,false]) then "touch null.dup.qc null.pbc.qc null.mito_dup.txt; " else ""}
		touch null
	}
	output {
		File nodup_bam = glob("*.bam")[0]
		File nodup_bai = glob("*.bai")[0]
		File flagstat_qc = glob("*.flagstat.qc")[0]
		File dup_qc = if select_first([no_dup_removal,false]) then glob("null")[0] else glob("*.dup.qc")[0]
		File pbc_qc = if select_first([no_dup_removal,false]) then glob("null")[0] else glob("*.pbc.qc")[0]
		File mito_dup_log = if select_first([no_dup_removal,false]) then glob("null")[0] else glob("*.mito_dup.txt")[0] # mito_dups, fract_dups_from_mito
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
	String? regex_grep_v_ta   	# Perl-style regular expression pattern 
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
			${"--regex-grep-v-ta " +"'"+select_first([regex_grep_v_ta,'chrM'])+"'"} \
			${"--subsample " + select_first([subsample,0])} \
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
		cpu : 1
		memory : "${select_first([mem_mb,'16000'])} MB"
		time : 1
		disks : "local-disk 50 HDD"
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
	runtime {
		cpu : 1
		memory : "4000 MB"
		time : 1
		disks : "local-disk 50 HDD"
	}
}

task xcor {
	# parameters from workflow
	File ta
	Boolean paired_end
	# optional
	Int? subsample  # number of reads to subsample TAGALIGN
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
		memory : "${select_first([mem_mb,'16000'])} MB"
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
	Float? pval_thresh  # p.value threshold
	Int? smooth_win 	# size of smoothing window
	Boolean? make_signal
	File blacklist 		# blacklist BED to filter raw peaks
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
			${"--pval-thresh "+ select_first([pval_thresh,'0.01'])} \
			${"--smooth-win "+ select_first([smooth_win,150])} \
			${if select_first([make_signal,false]) then "--make-signal" else ""} \
			${"--blacklist "+ blacklist}
		
		# ugly part to deal with optional outputs with Google JES backend
		${if select_first([make_signal,false]) then "" 
			else "touch null.pval.signal.bigwig null.fc.signal.bigwig"}
		touch null 
	}
	output {
		File npeak = glob("*[!.][!b][!f][!i][!l][!t].narrowPeak.gz")[0]
		File bfilt_npeak = glob("*.bfilt.narrowPeak.gz")[0]
		File bfilt_npeak_bb = glob("*.bfilt.narrowPeak.bb")[0]
		File sig_pval = if select_first([make_signal,false]) then glob("*.pval.signal.bigwig")[0] else glob("null")[0]
		File sig_fc = if select_first([make_signal,false]) then glob("*.fc.signal.bigwig")[0] else glob("null")[0]
		File frip_qc = glob("*.frip.qc")[0]
	}
	runtime {
		cpu : 1
		memory : "${select_first([mem_mb,'16000'])} MB"
		time : select_first([time_hr,24])
		disks : select_first([disks,"local-disk 100 HDD"])
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
	File? ta		# to calculate FRiP
	File chrsz			# 2-col chromosome sizes file
	String peak_type
	String rank

	command {
		python $(which encode_idr.py) \
			${peak1} ${peak2} ${peak_pooled} \
			${"--prefix " + prefix} \
			${"--idr-thresh " + select_first([idr_thresh,'0.1'])} \
			${"--peak-type " + peak_type} \
			--idr-rank ${rank} \
			${"--chrsz " + chrsz} \
			${"--blacklist "+ blacklist} \
			${"--ta " + ta}

		# ugly part to deal with optional outputs with Google backend
		${if defined(ta) then "" else "touch null.frip.qc"}
		touch null 
	}
	output {
		File idr_peak = glob("*[!.][!b][!f][!i][!l][!t]."+peak_type+".gz")[0]
		File bfilt_idr_peak = glob("*.bfilt."+peak_type+".gz")[0]
		File bfilt_idr_peak_bb = glob("*.bfilt."+peak_type+".bb")[0]
		File idr_plot = glob("*.txt.png")[0]
		File idr_unthresholded_peak = glob("*.txt.gz")[0]
		File idr_log = glob("*.log")[0]
		File frip_qc = if defined(ta) then glob("*.frip.qc")[0] else glob("null")[0]
	}
	runtime {
		cpu : 1
		memory : "8000 MB"
		time : 1
		disks : "local-disk 50 HDD"
	}	
}

task overlap {
	# parameters from workflow
	String prefix 		# prefix for IDR output file
	File peak1
	File peak2
	File peak_pooled
	File blacklist 	# blacklist BED to filter raw peaks
	File? ta		# to calculate FRiP
	File chrsz			# 2-col chromosome sizes file
	String peak_type

	command {
		python $(which encode_naive_overlap.py) \
			${peak1} ${peak2} ${peak_pooled} \
			${"--prefix " + prefix} \
			${"--peak-type " + peak_type} \
			${"--chrsz " + chrsz} \
			${"--blacklist "+ blacklist} \
			${"--ta " + ta}

		# ugly part to deal with optional outputs with Google backend
		${if defined(ta) then "" else "touch null.frip.qc"}
		touch null 
	}
	output {
		File overlap_peak = glob("*[!.][!b][!f][!i][!l][!t]."+peak_type+".gz")[0]
		File bfilt_overlap_peak = glob("*.bfilt."+peak_type+".gz")[0]
		File bfilt_overlap_peak_bb = glob("*.bfilt."+peak_type+".bb")[0]
		File frip_qc = if defined(ta) then glob("*.frip.qc")[0] else glob("null")[0]
	}
	runtime {
		cpu : 1
		memory : "4000 MB"
		time : 1
		disks : "local-disk 50 HDD"
	}
}

task reproducibility {
	# parameters from workflow
	String prefix
	Array[File]? peaks # peak files from pair of true replicates
						# in a sorted order. for example of 4 replicates,
						# 1,2 1,3 1,4 2,3 2,4 3,4.
                        # x,y means peak file from rep-x vs rep-y
	Array[File]? peaks_pr	# peak files from pseudo replicates
	File? peak_ppr			# Peak file from pooled pseudo replicate.
	String peak_type
	File chrsz			# 2-col chromosome sizes file

	command {
		python $(which encode_reproducibility_qc.py) \
			${sep=' ' peaks} \
			--peaks-pr ${sep=' ' peaks_pr} \
			${"--peak-ppr "+ peak_ppr} \
			--prefix ${prefix} \
			${"--peak-type " + peak_type} \
			${"--chrsz " + chrsz}
	}
	output {
		File optimal_peak = glob("optimal_peak.gz")[0]
		File conservative_peak = glob("conservative_peak.gz")[0]
		File optimal_peak_bb = glob("optimal_peak.*.bb")[0]
		File conservative_peak_bb = glob("conservative_peak.*.bb")[0]
		File reproducibility_qc = glob("*reproducibility.qc")[0]
	}
	runtime {
		cpu : 1
		memory : "4000 MB"
		time : 1
		disks : "local-disk 50 HDD"
	}
}

task ataqc { # generate ATAQC report
	Boolean paired_end
	File read_len_log
	File flagstat_log
	File bowtie2_log
	File bam
	File nodup_flagstat_log
	File mito_dup_log
	File dup_log
	File pbc_log
	File nodup_bam
	File ta
	File? peak
	File? idr_peak 
	File? overlap_peak
	File bigwig
	# from genome database
	File ref_fa
	File chrsz
	File tss_enrich
	File blacklist
	File dnase
	File prom
	File enh
	File reg2map_bed
	File reg2map
	File roadmap_meta
	# resource
	Int? mem_mb
	Int? time_hr
	String? disks

	command {
		export _JAVA_OPTIONS="-Xms256M -Xmx${select_first([mem_mb,'16000'])}M -XX:ParallelGCThreads=1"

		python $(which encode_ataqc.py) \
			${if paired_end then "--paired-end" else ""} \
			--read-len-log ${read_len_log} \
			--flagstat-log ${flagstat_log} \
			--bowtie2-log ${bowtie2_log} \
			--bam ${bam} \
			--nodup-flagstat-log ${nodup_flagstat_log} \
			--mito-dup-log ${mito_dup_log} \
			--dup-log ${dup_log} \
			--pbc-log ${pbc_log} \
			--nodup-bam ${nodup_bam} \
			--ta ${ta} \
			--bigwig ${bigwig} \
			${"--peak " + peak} \
			${"--idr-peak " + idr_peak} \
			${"--overlap-peak " + overlap_peak} \
			--ref-fa ${ref_fa} \
			--blacklist ${blacklist} \
			--chrsz ${chrsz} \
			--dnase ${dnase} \
			--tss-enrich ${tss_enrich} \
			--prom ${prom} \
			--enh ${enh} \
			--reg2map-bed ${reg2map_bed} \
			--reg2map ${reg2map} \
			--roadmap-meta ${roadmap_meta}
	}
	output {
	    #File raw_peak_summ = glob("*.raw_peak_summ.tsv")[0]
	    #File raw_peak_dist = glob("*.raw_peak_dist.png")[0]
	    #File naive_peak_summ =  if defined(overlap_peak) then glob("*.naive_peak_summ.tsv")[0] else glob("null")[0]
	    #File naive_peak_dist = if defined(overlap_peak) then glob("*.naive_peak_dist.png")[0] else glob("null")[0]
	    #File idr_peak_summ = if defined(idr_peak) then glob("*.idr_peak_summ.tsv")[0] else glob("null")[0]
	    #File idr_peak_dist = if defined(idr_peak) then glob("*.idr_peak_dist.png")[0] else glob("null")[0]
		#File preseq_data = glob("*.preseq.data")[0]
		#File preseq_log = glob("*.preseq.log")[0]
		#File gc_plot = glob("*_gcPlot.pdf")[0]
		#File gc_out = glob("*_gc.txt")[0]
		#File gc_summary = glob("*_gcSummary.txt")[0]
		# optional
		#File tss_plot_file = if defined(tss_enrich) then glob("*_tss-enrich.png")[0] else glob("null")[0]
		#File tss_plot_large_file = if defined(tss_enrich) then glob("*_large_tss-enrich.png")[0] else glob("null")[0]
		#File roadmap_compare_plot = glob("*.signal")[0]
		#File tar = glob("*.tar.gz")[0]
		File html = glob("*_qc.html")[0]
		File txt = glob("*_qc.txt")[0]
	}
	runtime {
		cpu : 1
		memory : "${select_first([mem_mb,'16000'])} MB"
		time : select_first([time_hr,24])
		disks : select_first([disks,"local-disk 100 HDD"])
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
	Array[File]? flagstat_qcs
	Array[File]? nodup_flagstat_qcs
	Array[File]? dup_qcs
	Array[File]? pbc_qcs
	Array[File]? xcor_plots
	Array[File]? xcor_scores
	Array[File]? idr_plots
	Array[File]? idr_plots_pr
	File? idr_plot_ppr
	Array[File]? frip_macs2_qcs
	Array[File]? frip_macs2_qcs_pr1
	Array[File]? frip_macs2_qcs_pr2
	File? frip_macs2_qc_pooled
	File? frip_macs2_qc_ppr1 
	File? frip_macs2_qc_ppr2 
	Array[File]? frip_idr_qcs
	Array[File]? frip_idr_qcs_pr
	File? frip_idr_qc_ppr 
	Array[File]? frip_overlap_qcs
	Array[File]? frip_overlap_qcs_pr
	File? frip_overlap_qc_ppr
	File? idr_reproducibility_qc
	File? overlap_reproducibility_qc

	File? qc_json_ref

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
			--frip-macs2-qcs ${sep=' ' frip_macs2_qcs} \
			--frip-macs2-qcs-pr1 ${sep=' ' frip_macs2_qcs_pr1} \
			--frip-macs2-qcs-pr2 ${sep=' ' frip_macs2_qcs_pr2} \
			${"--frip-macs2-qc-pooled " + frip_macs2_qc_pooled} \
			${"--frip-macs2-qc-ppr1 " + frip_macs2_qc_ppr1} \
			${"--frip-macs2-qc-ppr2 " + frip_macs2_qc_ppr2} \
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
		
		diff qc.json ${if defined(qc_json_ref) then qc_json_ref else "/dev/null"} | wc -l > qc_json_match.txt
	}
	output {
		File report = glob('*qc.html')[0]
		File qc_json = glob('*qc.json')[0]
		Boolean qc_json_match = read_int("qc_json_match.txt")==0
	}
	runtime {
		cpu : 1
		memory : "4000 MB"
		time : 1
		disks : "local-disk 50 HDD"		
	}
}

task read_genome_tsv {
	File genome_tsv
	command {
		cat ${genome_tsv}
	}
	output {
		Map[String,String] genome = read_map(stdout())
	}
	runtime {
		cpu : 1
		memory : "4000 MB"
		time : 1
		disks : "local-disk 50 HDD"		
	}
}

task compare_md5sum {
	Array[String] labels
	Array[File] files
	Array[File] ref_files

	command <<<
		python <<CODE	
		from collections import OrderedDict
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

		result = OrderedDict()
		match = OrderedDict()
		match_overall = True

		result['tasks'] = []
		result['failed_task_labels'] = []
		result['succeeded_task_labels'] = []
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
			result['tasks'].append(OrderedDict([
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
			if matched:
				result['succeeded_task_labels'].append(label)
			else:
				result['failed_task_labels'].append(label)		
		result['match_overall'] = match_overall

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
		String json_str = read_string('result.json') # details (string)
	}
	runtime {
		cpu : 1
		memory : "4000 MB"
		time : 1
		disks : "local-disk 50 HDD"		
	}
}
