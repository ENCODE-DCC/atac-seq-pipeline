# ENCODE DCC ATAC-Seq/DNase-Seq pipeline
# Author: Jin Lee (leepc12@gmail.com)

workflow atac {
	# pipeline version
	String pipeline_ver = 'v1.2.0'

	# general sample information
	String title = 'Untitled'
	String description = 'No description'

	# endedness for input data
	Boolean? paired_end				# to define endedness for all replciates
									#	if defined, this will override individual endedness below
	Array[Boolean] paired_ends = []	# to define endedness for individual replicate

	# genome TSV
	# 	you can define any genome parameters either in this TSV
	#	or individually in an input JSON file
	# 	individually defined parameters will override those defined in this TSV	
	File? genome_tsv 				# reference genome data TSV file including
									# all genome-specific file paths and parameters
	# individual genome parameters
	File? ref_fa					# reference fasta (*.fa.gz)
	File? bowtie2_idx_tar 			# bowtie2 index tar (uncompressed)
	File? chrsz 					# 2-col chromosome sizes file
	File? blacklist 				# blacklist BED (peaks overlapping will be filtered out)
	String? gensz 					# genome sizes (hs for human, mm for mouse or sum of 2nd col in chrsz)
	# individual genome parameters for ATAqC
	File? tss 						# TSS BED file
	File? dnase 					# open chromatin region BED file
	File? prom 						# promoter region BED file
	File? enh 						# enhancer region BED file
	File? reg2map 					# file with cell type signals
	File? reg2map_bed 				# file of regions used to generate reg2map signals
	File? roadmap_meta 				# roadmap metedata

	# parameters for pipeline
	String pipeline_type = 'atac'	# atac (default), dnase
									# tn5 shiting will be enabled for atac only
	Boolean align_only = false		# disable all post-align analyses (peak-calling, overlap, idr, ...)
	Boolean true_rep_only = false 	# disable all analyses for pseudo replicates
									# 	if activated, overlap and idr will be disabled

	# parameters for trim_adapter
	Boolean auto_detect_adapter = false
									# automatically detect/trim adapters
									# 	can detect three adapters only
									# 	see /src/detect_adapter.py for details
	String cutadapt_param = '-e 0.1 -m 5'
									# cutadapt parameters (err_rate=0.1, min_trim_len=5)

	# parameters for align (align FASTQs and create raw BAM)
	#String aligner = 'bowtie2' 		# bowtie2, custom
	Int multimapping = 0			# for samples with multimapping reads
	String bowtie2_param_se = '--local'
									# params for bowtie2 (single-ended samples)
	String bowtie2_param_pe = '-X2000 --mm --local' 
									# params for bowtie2 (paired-ended samples)

	# parameters for filter (filter/dedup raw BAM)
	String dup_marker = 'picard'	# picard, sambamba
	Int mapq_thresh = 30			# threshold for low MAPQ reads removal
	Boolean no_dup_removal = false 	# keep all dupes in final BAM
									
	# parameters for bam2ta (convert filtered/deduped BAM to TAG-ALIGN)
	String mito_chr_name = 'chrM' 	# name of mito chromosome
									# 	THIS IS NOT A REG-EX!
									#	you can define only one name for mito chrom
	String regex_filter_reads = 'chrM' 	
									# Perl-style regular expression pattern 
									#	for chr name to filter out reads
									# 	THIS IS A REG-EX!
									# 	you can define multiple chrom names e.g. \(chrM\|chr21\|chr19\)
									# 	make sure that you escape (, ) and | with \
	Int subsample_reads = 0			# number of reads to subsample TAGALIGN
									# 	0 for no subsampling

	# parameters for cross-correlation analysis
	Boolean enable_xcor = false 	# enable cross-corr analysis
	Int xcor_subsample_reads = 25000000	
									# number of reads to subsample TAG-ALIGN
									# 	this will be used for cross-corr only
									# 	will not affect any downstream analyses

	# parameters for blacklist filtering peaks
	Boolean keep_irregular_chr_in_bfilt_peak = false 
									# peaks with irregular chr name will not be filtered out
									# 	in bfilt_peak (blacklist filtered peak) file
									# 	(e.g. chr1_AABBCC, AABR07024382.1, ...)
									# 	reg-ex pattern for "regular" chr name is chr[\dXY]+\b

	# parameters for peak calling
	String peak_type = 'narrowPeak'
	Int cap_num_peak = 300000		# cap number of raw peaks called
	Float pval_thresh = 0.01		# p.value threshold for peak caller
	Int smooth_win = 150			# size of smoothing window for peak caller

	# parameters for signal tracks
	Boolean enable_count_signal_track = false # generate count signal track

	# parameters for IDR
	Boolean enable_idr = false 		# enable IDR analysis on raw peaks
	Float idr_thresh = 0.1			# IDR threshold
	String idr_rank = 'p.value' 	# IDR ranking method (p.value, q.value, score)

	# parameters for ATAqC
	Boolean disable_ataqc = false 	# disable ATAqC (extra annotation-based analysis)

	# resources 	
	#	these variables will be automatically ignored if they are not supported by platform
	# 	"disks" is for cloud platforms (Google Cloud Platform, DNAnexus) only
	Int trim_adapter_cpu = 2
	Int trim_adapter_mem_mb = 12000
	Int trim_adapter_time_hr = 24
	String trim_adapter_disks = "local-disk 100 HDD"

	Int bowtie2_cpu = 4
	Int bowtie2_mem_mb = 20000
	Int bowtie2_time_hr = 48
	String bowtie2_disks = "local-disk 100 HDD"

	Int filter_cpu = 2
	Int filter_mem_mb = 20000
	Int filter_time_hr = 24
	String filter_disks = "local-disk 100 HDD"

	Int bam2ta_cpu = 2
	Int bam2ta_mem_mb = 10000
	Int bam2ta_time_hr = 6
	String bam2ta_disks = "local-disk 100 HDD"

	Int spr_mem_mb = 16000

	Int xcor_cpu = 2
	Int xcor_mem_mb = 16000
	Int xcor_time_hr = 6
	String xcor_disks = "local-disk 100 HDD"

	Int macs2_mem_mb = 16000
	Int macs2_time_hr = 24
	String macs2_disks = "local-disk 100 HDD"

	Int ataqc_mem_mb = 16000
	Int ataqc_mem_java_mb = 15000
	Int ataqc_time_hr = 24
	String ataqc_disks = "local-disk 400 HDD"

	# input file definition
	# supported types: fastq, bam, nodup_bam (or filtered bam), ta (tagAlign), peak
	# 	pipeline can start from any type of inputs
	# 	leave all other types undefined
	# 	you can define up to 10 replicates

 	# fastqs and adapters  	
	# 	if auto_detect_adapter == true, undefined adapters will be detected/trimmed 
	# 	otherwise, only defined adapters will be trimmed
	# 	so you can selectively detect/trim adapters for a specific fastq
	Array[File] fastqs_rep1_R1 = []		# FASTQs to be merged for rep1 R1
	Array[File] fastqs_rep1_R2 = [] 	# do not define if single-ended
	Array[File] fastqs_rep2_R1 = [] 	# do not define if unreplicated
	Array[File] fastqs_rep2_R2 = []		# ...
	Array[File] fastqs_rep3_R1 = []
	Array[File] fastqs_rep3_R2 = []
	Array[File] fastqs_rep4_R1 = []
	Array[File] fastqs_rep4_R2 = []
	Array[File] fastqs_rep5_R1 = []
	Array[File] fastqs_rep5_R2 = []
	Array[File] fastqs_rep6_R1 = []
	Array[File] fastqs_rep6_R2 = []
	Array[File] fastqs_rep7_R1 = []
	Array[File] fastqs_rep7_R2 = []
	Array[File] fastqs_rep8_R1 = []
	Array[File] fastqs_rep8_R2 = []
	Array[File] fastqs_rep9_R1 = []
	Array[File] fastqs_rep9_R2 = []
	Array[File] fastqs_rep10_R1 = []
	Array[File] fastqs_rep10_R2 = []

	String? adapter
	Array[String] adapters_rep1_R1 = []
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
	Array[String] adapters_rep7_R1 = []
	Array[String] adapters_rep7_R2 = []
	Array[String] adapters_rep8_R1 = []
	Array[String] adapters_rep8_R2 = []
	Array[String] adapters_rep9_R1 = []
	Array[String] adapters_rep9_R2 = []
	Array[String] adapters_rep10_R1 = []
	Array[String] adapters_rep10_R2 = []

	# other input types (bam, nodup_bam, ta). they are per replicate
	Array[File?] trim_merged_fastqs_R1 = []
	Array[File?] trim_merged_fastqs_R2 = [] 
	Array[File?] bams = []
	Array[File?] nodup_bams = []
	Array[File?] tas = []

	# other input types (peak)
	Array[File?] peaks = []			# per replicate
	Array[File?] peaks_pr1 = []		# per replicate. do not define if true_rep_only==true
	Array[File?] peaks_pr2 = []		# per replicate. do not define if true_rep_only==true
	File? peak_ppr1					# do not define if unreplicated or true_rep_only==true
	File? peak_ppr2					# do not define if unreplicated or true_rep_only==true
	File? peak_pooled				# do not define if unreplicated or true_rep_only==true

	####################### pipeline starts here #######################
	# DO NOT DEFINE ANY VARIABLES DECLARED BELOW IN AN INPUT JSON FILE #
	# THEY ARE TEMPORARY/INTERMEDIATE SYSTEM VARIABLES                 #
	####################### pipeline starts here #######################
	
	# read genome data and paths
	if ( defined(genome_tsv) ) {
		call read_genome_tsv { input: genome_tsv = genome_tsv }
	}
	File? ref_fa_ = if defined(ref_fa) then ref_fa
		else read_genome_tsv.ref_fa
	File? bowtie2_idx_tar_ = if defined(bowtie2_idx_tar) then bowtie2_idx_tar
		else read_genome_tsv.bowtie2_idx_tar
	File? chrsz_ = if defined(chrsz) then chrsz
		else read_genome_tsv.chrsz
	String? gensz_ = if defined(gensz) then gensz
		else read_genome_tsv.gensz
	File? blacklist_ = if defined(blacklist) then blacklist
		else read_genome_tsv.blacklist
	File? tss_ = if defined(tss) then tss
		else read_genome_tsv.tss
	File? dnase_ = if defined(dnase) then dnase
		else read_genome_tsv.dnase
	File? prom_ = if defined(prom) then prom
		else read_genome_tsv.prom
	File? enh_ = if defined(enh) then enh
		else read_genome_tsv.enh
	File? reg2map_ = if defined(reg2map) then reg2map
		else read_genome_tsv.reg2map
	File? reg2map_bed_ = if defined(reg2map_bed) then reg2map_bed
		else read_genome_tsv.reg2map_bed
	File? roadmap_meta_ = if defined(roadmap_meta) then roadmap_meta
		else read_genome_tsv.roadmap_meta

	# temporary files used for resumer, output from task align
	Array[File?] bowtie2__align_log = []
	Array[File?] bowtie2__flagstat_qc = []
	Array[File?] bowtie2__read_len_log = []

	Array[File?] filter__flagstat_qc = []
	Array[File?] filter__dup_qc = []
	Array[File?] filter__pbc_qc = []
	Array[File?] filter__mito_dup_log = []

	Array[File?] spr__ta_pr1 = []
	Array[File?] spr__ta_pr2 = []

	Array[File?] xcor__plot_png = []
	Array[File?] xcor__score = []

	Array[File?] count_signal_track__pos_bw = []
	Array[File?] count_signal_track__neg_bw = []

	Array[File?] macs2_signal_track__pval_bw = []
	Array[File?] macs2_signal_track__fc_bw = []

	Array[File?] macs2__frip_qc = []
	Array[File?] macs2_pr1__frip_qc = []
	Array[File?] macs2_pr2__frip_qc = []

	Array[File?] idr__idr_peak = [] # per pair of replicate
	Array[File?] idr__bfilt_idr_peak = []
	Array[File?] idr__idr_plot = []
	Array[File?] idr__frip_qc = []

	Array[File?] idr_pr__idr_peak = []	# per replicate
	Array[File?] idr_pr__bfilt_idr_peak = []
	Array[File?] idr_pr__idr_plot = []
	Array[File?] idr_pr__frip_qc = []

	Array[File?] overlap__overlap_peak = []	# per pair of replicate
	Array[File?] overlap__bfilt_overlap_peak = []
	Array[File?] overlap__frip_qc = []

	Array[File?] overlap_pr__overlap_peak = [] # per replicate
	Array[File?] overlap_pr__bfilt_overlap_peak = []
	Array[File?] overlap_pr__frip_qc = []

	Array[File?] ataqc__html = []
	Array[File?] ataqc__txt = []

	File? pool_ta__ta_pooled
	File? pool_ta_pr1__ta_pooled
	File? pool_ta_pr2__ta_pooled

	File? count_signal_track_pooled__pos_bw
	File? count_signal_track_pooled__neg_bw

	File? macs2_signal_track_pooled__pval_bw
	File? macs2_signal_track_pooled__fc_bw

	File? macs2_pooled__frip_qc
	File? macs2_ppr1__frip_qc
	File? macs2_ppr2__frip_qc

	File? idr_ppr__idr_peak
	File? idr_ppr__bfilt_idr_peak
	File? idr_ppr__idr_plot
	File? idr_ppr__frip_qc
	
	File? overlap_ppr__overlap_peak
	File? overlap_ppr__bfilt_overlap_peak
	File? overlap_ppr__frip_qc

	File? reproducibility_idr__optimal_peak
	File? reproducibility_idr__reproducibility_qc

	File? reproducibility_overlap__optimal_peak
	File? reproducibility_overlap__reproducibility_qc

	# temporary 2-dim fastqs array [rep_id][merge_id]
	Array[Array[File]] fastqs_R1 = 
		if length(fastqs_rep10_R1)>0 then
			[fastqs_rep1_R1, fastqs_rep2_R1, fastqs_rep3_R1, fastqs_rep4_R1, fastqs_rep5_R1,
			fastqs_rep6_R1, fastqs_rep7_R1, fastqs_rep8_R1, fastqs_rep9_R1, fastqs_rep10_R1]
		else if length(fastqs_rep9_R1)>0 then
			[fastqs_rep1_R1, fastqs_rep2_R1, fastqs_rep3_R1, fastqs_rep4_R1, fastqs_rep5_R1,
			fastqs_rep6_R1, fastqs_rep7_R1, fastqs_rep8_R1, fastqs_rep9_R1]
		else if length(fastqs_rep8_R1)>0 then
			[fastqs_rep1_R1, fastqs_rep2_R1, fastqs_rep3_R1, fastqs_rep4_R1, fastqs_rep5_R1,
			fastqs_rep6_R1, fastqs_rep7_R1, fastqs_rep8_R1]
		else if length(fastqs_rep7_R1)>0 then
			[fastqs_rep1_R1, fastqs_rep2_R1, fastqs_rep3_R1, fastqs_rep4_R1, fastqs_rep5_R1,
			fastqs_rep6_R1, fastqs_rep7_R1]
		else if length(fastqs_rep6_R1)>0 then
			[fastqs_rep1_R1, fastqs_rep2_R1, fastqs_rep3_R1, fastqs_rep4_R1, fastqs_rep5_R1,
			fastqs_rep6_R1]
		else if length(fastqs_rep5_R1)>0 then
			[fastqs_rep1_R1, fastqs_rep2_R1, fastqs_rep3_R1, fastqs_rep4_R1, fastqs_rep5_R1]
		else if length(fastqs_rep4_R1)>0 then
			[fastqs_rep1_R1, fastqs_rep2_R1, fastqs_rep3_R1, fastqs_rep4_R1]
		else if length(fastqs_rep3_R1)>0 then
			[fastqs_rep1_R1, fastqs_rep2_R1, fastqs_rep3_R1]
		else if length(fastqs_rep2_R1)>0 then
			[fastqs_rep1_R1, fastqs_rep2_R1]
		else if length(fastqs_rep1_R1)>0 then
			[fastqs_rep1_R1]
		else []
	# no need to do that for R2 (R1 array will be used to determine presense of fastq for each rep)
	Array[Array[File]] fastqs_R2 = 
		[fastqs_rep1_R2, fastqs_rep2_R2, fastqs_rep3_R2, fastqs_rep4_R2, fastqs_rep5_R2,
		fastqs_rep6_R2, fastqs_rep7_R2, fastqs_rep8_R2, fastqs_rep9_R2, fastqs_rep10_R2]

	# temporary 2-dim adapters array [rep_id][merge_id]
	Array[Array[String]] adapters_R1 = 
		[adapters_rep1_R1, adapters_rep2_R1, adapters_rep3_R1, adapters_rep4_R1, adapters_rep5_R1,
		adapters_rep6_R1, adapters_rep7_R1, adapters_rep8_R1, adapters_rep9_R1, adapters_rep10_R1]
	Array[Array[String]] adapters_R2 = 
		[adapters_rep1_R2, adapters_rep2_R2, adapters_rep3_R2, adapters_rep4_R2, adapters_rep5_R2,
		adapters_rep6_R2, adapters_rep7_R2, adapters_rep8_R2, adapters_rep9_R2, adapters_rep10_R2]

	# temporary variables to get number of replicates
	# 	WDLic implementation of max(A,B,C,...)
	Int num_rep_fastq = length(fastqs_R1)
	Int num_rep_trim_merged_fastq = if length(trim_merged_fastqs_R1)<num_rep_fastq then num_rep_fastq
		else length(trim_merged_fastqs_R1)
	Int num_rep_bam = if length(bams)<num_rep_trim_merged_fastq then num_rep_trim_merged_fastq
		else length(bams)
	Int num_rep_nodup_bam = if length(nodup_bams)<num_rep_bam then num_rep_bam
		else length(nodup_bams)
	Int num_rep_ta = if length(tas)<num_rep_nodup_bam then num_rep_nodup_bam
		else length(tas)
	Int num_rep_peak = if length(peaks)<num_rep_ta then num_rep_ta
		else length(peaks)
	Int num_rep = num_rep_peak

	# align each replicate
	scatter(i in range(num_rep)) {
		# to override endedness definition for individual replicate
		# 	paired_end will override paired_ends[i]
		Boolean? paired_end_ = if !defined(paired_end) && i<length(paired_ends) then paired_ends[i]
			else paired_end

		Boolean has_input_of_trim_adapter = i<length(fastqs_R1) && length(fastqs_R1[i])>0		
		Boolean has_output_of_trim_adapter = i<length(trim_merged_fastqs_R1) &&
			defined(trim_merged_fastqs_R1[i])

		# skip if we already have output of this step
		if ( has_input_of_trim_adapter && !has_output_of_trim_adapter ) {
			call trim_adapter { input :
				fastqs_R1 = fastqs_R1[i],
				fastqs_R2 = fastqs_R2[i],
				adapter = adapter,
				adapters_R1 = adapters_R1[i],
				adapters_R2 = adapters_R2[i],
				paired_end = paired_end_,
				auto_detect_adapter = auto_detect_adapter,
				cutadapt_param = cutadapt_param,
				# resource
				cpu = trim_adapter_cpu,
				mem_mb = trim_adapter_mem_mb,
				time_hr = trim_adapter_time_hr,
				disks = trim_adapter_disks,
			}
		}

		File? trim_merged_fastq_R1_ = if has_output_of_trim_adapter
			then trim_merged_fastqs_R1[i]
			else trim_adapter.trim_merged_fastq_R1
		File? trim_merged_fastq_R2_ = if i<length(trim_merged_fastqs_R2) &&
			defined(trim_merged_fastqs_R2[i])
			then trim_merged_fastqs_R2[i]
			else trim_adapter.trim_merged_fastq_R2

		Boolean has_input_of_bowtie2 = has_output_of_trim_adapter ||
			defined(trim_adapter.trim_merged_fastq_R1)
		Boolean has_output_of_bowtie2 = i<length(bams) && defined(bams[i])

		if ( has_input_of_bowtie2 && !has_output_of_bowtie2 ) {
			call bowtie2 { input :
				fastq_R1 = trim_merged_fastq_R1_,
				fastq_R2 = trim_merged_fastq_R2_,
				paired_end = paired_end_,
				#aligner = aligner,
				multimapping = multimapping,
				idx_tar = bowtie2_idx_tar_,
				bowtie2_param_se = bowtie2_param_se,
				bowtie2_param_pe = bowtie2_param_pe,
				# resource
				cpu = bowtie2_cpu,
				mem_mb = bowtie2_mem_mb,
				time_hr = bowtie2_time_hr,
				disks = bowtie2_disks,
			}
		}

		File? bam_ = if has_output_of_bowtie2 then bams[i]
			else bowtie2.bam
		File? flagstat_qc_ = if i<length(bowtie2__flagstat_qc) &&
			defined(bowtie2__flagstat_qc[i]) then bowtie2__flagstat_qc[i]
			else bowtie2.flagstat_qc
		File? align_log_ = if i<length(bowtie2__align_log) &&
			defined(bowtie2__align_log[i]) then bowtie2__align_log[i]
			else bowtie2.align_log
		File? read_len_log_ = if i<length(bowtie2__read_len_log) &&
			defined(bowtie2__read_len_log[i]) then bowtie2__read_len_log[i]
			else bowtie2.read_len_log

		Boolean has_input_of_filter = has_output_of_bowtie2 || defined(bowtie2.bam)
		Boolean has_output_of_filter = i<length(nodup_bams) && defined(nodup_bams[i])
		
		# skip if we already have output of this step
		if ( has_input_of_filter && !has_output_of_filter ) {
			call filter { input :
				bam = bam_,
				paired_end = paired_end_,
				dup_marker = dup_marker,
				mapq_thresh = mapq_thresh,
				no_dup_removal = no_dup_removal,
				multimapping = multimapping,
				mito_chr_name = mito_chr_name,

				cpu = filter_cpu,
				mem_mb = filter_mem_mb,
				time_hr = filter_time_hr,
				disks = filter_disks,
			}
		}

		File? nodup_bam_ = if has_output_of_filter then nodup_bams[i]
			else filter.nodup_bam
		File? nodup_flagstat_qc_ = if i<length(filter__flagstat_qc) &&
			defined(filter__flagstat_qc[i]) then filter__flagstat_qc[i]
			else filter.flagstat_qc
		File? dup_qc_ = if i<length(filter__dup_qc) &&
			defined(filter__dup_qc[i]) then filter__dup_qc[i]
			else filter.dup_qc
		File? pbc_qc_ = if i<length(filter__pbc_qc) &&
			defined(filter__pbc_qc[i]) then filter__pbc_qc[i]
			else filter.pbc_qc
		File? mito_dup_log_ = if i<length(filter__mito_dup_log) &&
			defined(filter__mito_dup_log[i]) then filter__mito_dup_log[i]
			else filter.mito_dup_log

		Boolean has_input_of_bam2ta = has_output_of_filter || defined(filter.nodup_bam)
		Boolean has_output_of_bam2ta = i<length(tas) && defined(tas[i])

		if ( has_input_of_bam2ta && !has_output_of_bam2ta ) {
			call bam2ta { input :
				bam = nodup_bam_,
				disable_tn5_shift = if pipeline_type=='atac' then false else true,
				regex_grep_v_ta = regex_filter_reads,
				subsample = subsample_reads,
				paired_end = paired_end_,
				mito_chr_name = mito_chr_name,

				cpu = bam2ta_cpu,
				mem_mb = bam2ta_mem_mb,
				time_hr = bam2ta_time_hr,
				disks = bam2ta_disks,
			}
		}

		File? ta_ = if has_output_of_bam2ta then tas[i] else bam2ta.ta

		Boolean has_input_of_xcor = has_output_of_bam2ta || defined(bam2ta.ta)
		Boolean has_output_of_xcor = i<length(xcor__score) && defined(xcor__score[i])

		if ( has_input_of_xcor && !has_output_of_xcor && 
			enable_xcor ) {
			# subsample tagalign (non-mito) and cross-correlation analysis
			call xcor { input :
				ta = ta_,
				subsample = xcor_subsample_reads,
				paired_end = paired_end_,
				mito_chr_name = mito_chr_name,

				cpu = xcor_cpu,
				mem_mb = xcor_mem_mb,
				time_hr = xcor_time_hr,
				disks = xcor_disks,
			}
		}

		File? xcor_score_ = if has_output_of_xcor then xcor__score[i] else xcor.score
		File? xcor_plot_ = if i<length(xcor__plot_png) &&
			defined(xcor__plot_png[i]) then xcor__plot_png[i] else xcor.plot_png

		Boolean has_input_of_macs2_signal_track = has_output_of_bam2ta || defined(bam2ta.ta)
		Boolean has_output_of_macs2_signal_track = i<length(macs2_signal_track__pval_bw) && 
			defined(macs2_signal_track__pval_bw[i])

		if ( has_input_of_macs2_signal_track && !has_output_of_macs2_signal_track ) {
			# generate count signal track
			call macs2_signal_track { input :
				ta = ta_,
				gensz = gensz_,
				chrsz = chrsz_,
				pval_thresh = pval_thresh,
				smooth_win = smooth_win,

				mem_mb = macs2_mem_mb,
				disks = macs2_disks,
				time_hr = macs2_time_hr,
			}
		}

		File? pval_bw_ = if has_output_of_macs2_signal_track then macs2_signal_track__pval_bw[i]
			else macs2_signal_track.pval_bw
		File? fc_bw_ = if i<length(macs2_signal_track__fc_bw) &&
			defined(macs2_signal_track__fc_bw[i]) then macs2_signal_track__fc_bw[i]
			else macs2_signal_track.fc_bw

		Boolean has_input_of_macs2 = has_output_of_bam2ta || defined(bam2ta.ta)
		Boolean has_output_of_macs2 = i<length(peaks) && defined(peaks[i])

		if ( has_input_of_macs2 && !has_output_of_macs2 &&
			!align_only ) {
			# call peaks on tagalign
			call macs2 { input :
				ta = ta_,
				gensz = gensz_,
				chrsz = chrsz_,
				cap_num_peak = cap_num_peak,
				pval_thresh = pval_thresh,
				smooth_win = smooth_win,
				blacklist = blacklist_,
				keep_irregular_chr_in_bfilt_peak = keep_irregular_chr_in_bfilt_peak,

				mem_mb = macs2_mem_mb,
				disks = macs2_disks,
				time_hr = macs2_time_hr,
			}
		}

		File? peak_ = if has_output_of_macs2 then peaks[i] else macs2.npeak
		File? macs2_frip_qc_ = if i<length(macs2__frip_qc) &&
			defined(macs2__frip_qc[i]) then macs2__frip_qc[i] else macs2.frip_qc

		Boolean has_input_of_spr = has_output_of_bam2ta || defined(bam2ta.ta)
		Boolean has_output_of_spr = i<length(spr__ta_pr1) && defined(spr__ta_pr1[i])

		if ( has_input_of_spr && !has_output_of_spr &&
			!align_only && !true_rep_only ) {
			call spr { input :
				ta = ta_,
				paired_end = paired_end_,
				mem_mb = spr_mem_mb,
			}
		}

		File? ta_pr1_ = if has_output_of_spr then spr__ta_pr1[i] else spr.ta_pr1
		File? ta_pr2_ = if i<length(spr__ta_pr2) &&
			defined(spr__ta_pr2[i]) then spr__ta_pr2[i]
			else spr.ta_pr2

		Boolean has_input_of_macs2_pr1 = has_output_of_spr || defined(spr.ta_pr1)
		Boolean has_output_of_macs2_pr1 = i<length(peaks_pr1) && defined(peaks_pr1[i])

		if ( has_input_of_macs2_pr1 && !has_output_of_macs2_pr1 &&
			!align_only && !true_rep_only ) {
			# call peaks on 1st pseudo replicated tagalign 
			call macs2 as macs2_pr1 { input :
				ta = ta_pr1_,
				gensz = gensz_,
				chrsz = chrsz_,
				cap_num_peak = cap_num_peak,
				pval_thresh = pval_thresh,
				smooth_win = smooth_win,
				blacklist = blacklist_,
				keep_irregular_chr_in_bfilt_peak = keep_irregular_chr_in_bfilt_peak,

				mem_mb = macs2_mem_mb,
				disks = macs2_disks,
				time_hr = macs2_time_hr,
			}
		}

		File? peak_pr1_ = if has_output_of_macs2_pr1 then peaks_pr1[i]
			else macs2_pr1.npeak
		File? macs2_pr1_frip_qc_ = if i<length(macs2_pr1__frip_qc) &&
			defined(macs2_pr1__frip_qc[i]) then macs2_pr1__frip_qc[i]
			else macs2_pr1.frip_qc

		Boolean has_input_of_macs2_pr2 = has_output_of_spr || defined(spr.ta_pr2)
		Boolean has_output_of_macs2_pr2 = i<length(peaks_pr2) && defined(peaks_pr2[i])

		if ( has_input_of_macs2_pr2 && !has_output_of_macs2_pr2 &&
			!align_only && !true_rep_only ) {
			# call peaks on 2nd pseudo replicated tagalign 
			call macs2 as macs2_pr2 { input :
				ta = ta_pr2_,
				gensz = gensz_,
				chrsz = chrsz_,
				cap_num_peak = cap_num_peak,
				pval_thresh = pval_thresh,
				smooth_win = smooth_win,
				blacklist = blacklist_,
				keep_irregular_chr_in_bfilt_peak = keep_irregular_chr_in_bfilt_peak,

				mem_mb = macs2_mem_mb,
				disks = macs2_disks,
				time_hr = macs2_time_hr,
			}
		}

		File? peak_pr2_ = if has_output_of_macs2_pr2 then peaks_pr2[i]
			else macs2_pr2.npeak
		File? macs2_pr2_frip_qc_ = if i<length(macs2_pr2__frip_qc) &&
			defined(macs2_pr2__frip_qc[i]) then macs2_pr2__frip_qc[i]
			else macs2_pr2.frip_qc

		Boolean has_input_of_count_signal_track = has_output_of_bam2ta || defined(bam2ta.ta)
		Boolean has_output_of_count_signal_track = i<length(count_signal_track__pos_bw) && 
			defined(count_signal_track__pos_bw[i])

		if ( has_input_of_count_signal_track && !has_output_of_count_signal_track &&
			enable_count_signal_track ) {
			# generate count signal track
			call count_signal_track { input :
				ta = ta_,
				chrsz = chrsz_,
			}
		}

		File? pos_bw_ = if has_output_of_count_signal_track then count_signal_track__pos_bw[i]
			else count_signal_track.pos_bw
		File? neg_bw_ = if i<length(count_signal_track__neg_bw) &&
			defined(count_signal_track__neg_bw[i]) then count_signal_track__neg_bw[i]
			else count_signal_track.neg_bw
	}

	# if there are TAs for ALL replicates then pool them
	Boolean has_all_inputs_of_pool_ta = length(select_all(ta_))==num_rep
	Boolean has_output_of_pool_ta = defined(pool_ta__ta_pooled)

	if ( has_all_inputs_of_pool_ta && !has_output_of_pool_ta &&
		num_rep>1 ) {
		# pool tagaligns from true replicates
		call pool_ta { input :
			tas = ta_,
		}
	}

	File? ta_pooled_ = if has_output_of_pool_ta then pool_ta__ta_pooled
		else pool_ta.ta_pooled	

	# if there are pr1 TAs for ALL replicates then pool them
	Boolean has_all_inputs_of_pool_ta_pr1 = length(select_all(ta_pr1_))==num_rep
	Boolean has_output_of_pool_ta_pr1 = defined(pool_ta_pr1__ta_pooled)

	if ( has_all_inputs_of_pool_ta_pr1 && !has_output_of_pool_ta_pr1 &&
		!align_only && !true_rep_only && num_rep>1 ) {
		# pool tagaligns from pseudo replicate 1
		call pool_ta as pool_ta_pr1 { input :
			tas = ta_pr1_,
		}
	}

	File? ta_ppr1_ = if has_output_of_pool_ta_pr1 then pool_ta_pr1__ta_pooled
		else pool_ta_pr1.ta_pooled	

	# if there are pr2 TAs for ALL replicates then pool them
	Boolean has_all_inputs_of_pool_ta_pr2 = length(select_all(ta_pr2_))==num_rep
	Boolean has_output_of_pool_ta_pr2 = defined(pool_ta_pr2__ta_pooled)

	if ( has_all_inputs_of_pool_ta_pr1 && !has_output_of_pool_ta_pr1 &&
		!align_only && !true_rep_only && num_rep>1 ) {
		# pool tagaligns from pseudo replicate 2
		call pool_ta as pool_ta_pr2 { input :
			tas = ta_pr2_,
		}
	}

	File? ta_ppr2_ = if has_output_of_pool_ta_pr2 then pool_ta_pr2__ta_pooled
		else pool_ta_pr2.ta_pooled	

	Boolean has_input_of_macs2_pooled = has_output_of_pool_ta || defined(pool_ta.ta_pooled)
	Boolean has_output_of_macs2_pooled = defined(peak_pooled)

	if ( has_input_of_macs2_pooled && !has_output_of_macs2_pooled &&
		!align_only && num_rep>1 ) {
		# call peaks on pooled replicate
		call macs2 as macs2_pooled { input :
			ta = ta_pooled_,
			gensz = gensz_,
			chrsz = chrsz_,
			cap_num_peak = cap_num_peak,
			pval_thresh = pval_thresh,
			smooth_win = smooth_win,
			blacklist = blacklist_,
			keep_irregular_chr_in_bfilt_peak = keep_irregular_chr_in_bfilt_peak,

			mem_mb = macs2_mem_mb,
			disks = macs2_disks,
			time_hr = macs2_time_hr,
		}
	}

	File? peak_pooled_ = if has_output_of_macs2_pooled then peak_pooled
		else macs2_pooled.npeak
	File? frip_macs2_qc_pooled_ = if defined(macs2_pooled__frip_qc) then macs2_pooled__frip_qc
		else macs2_pooled.frip_qc

	Boolean has_input_of_count_signal_track_pooled = has_output_of_pool_ta || defined(pool_ta.ta_pooled)
	Boolean has_output_of_count_signal_track_pooled = defined(count_signal_track_pooled__pos_bw)

	if ( has_input_of_count_signal_track_pooled && !has_output_of_count_signal_track_pooled &&
		enable_count_signal_track && num_rep>1 ) {
		call count_signal_track as count_signal_track_pooled { input :
			ta = ta_pooled_,
			chrsz = chrsz_,
		}
	}

	Boolean has_input_of_macs2_signal_track_pooled = has_output_of_pool_ta || defined(pool_ta.ta_pooled)
	Boolean has_output_of_macs2_signal_track_pooled = defined(macs2_signal_track_pooled__pval_bw)

	if ( has_input_of_macs2_signal_track_pooled && !has_output_of_macs2_signal_track_pooled &&
		num_rep>1 ) {
		call macs2_signal_track as macs2_signal_track_pooled { input :
			ta = ta_pooled_,
			gensz = gensz_,
			chrsz = chrsz_,
			pval_thresh = pval_thresh,
			smooth_win = smooth_win,

			mem_mb = macs2_mem_mb,
			disks = macs2_disks,
			time_hr = macs2_time_hr,
		}
	}

	Boolean has_input_of_macs2_ppr1 = has_output_of_pool_ta_pr1 || defined(pool_ta_pr1.ta_pooled)
	Boolean has_output_of_macs2_ppr1 = defined(peak_ppr1)

	if ( has_input_of_macs2_ppr1 && !has_output_of_macs2_ppr1 &&
		!align_only && !true_rep_only && num_rep>1 ) {
		# call peaks on 1st pooled pseudo replicates
		call macs2 as macs2_ppr1 { input :
			ta = ta_ppr1_,
			gensz = gensz_,
			chrsz = chrsz_,
			cap_num_peak = cap_num_peak,
			pval_thresh = pval_thresh,
			smooth_win = smooth_win,
			blacklist = blacklist_,
			keep_irregular_chr_in_bfilt_peak = keep_irregular_chr_in_bfilt_peak,

			mem_mb = macs2_mem_mb,
			disks = macs2_disks,
			time_hr = macs2_time_hr,
		}
	}

	File? peak_ppr1_ = if has_output_of_macs2_ppr1 then peak_ppr1
		else macs2_ppr1.npeak
	File? frip_macs2_qc_ppr1_ = if defined(macs2_ppr1__frip_qc) then macs2_ppr1__frip_qc
		else macs2_ppr1.frip_qc

	Boolean has_input_of_macs2_ppr2 = has_output_of_pool_ta_pr2 || defined(pool_ta_pr2.ta_pooled)
	Boolean has_output_of_macs2_ppr2 = defined(peak_ppr2)

	if ( has_input_of_macs2_ppr2 && !has_output_of_macs2_ppr2 &&
		!align_only && !true_rep_only && num_rep>1 ) {
		# call peaks on 2nd pooled pseudo replicates
		call macs2 as macs2_ppr2 { input :
			ta = ta_ppr2_,
			gensz = gensz_,
			chrsz = chrsz_,
			cap_num_peak = cap_num_peak,
			pval_thresh = pval_thresh,
			smooth_win = smooth_win,
			blacklist = blacklist_,
			keep_irregular_chr_in_bfilt_peak = keep_irregular_chr_in_bfilt_peak,

			mem_mb = macs2_mem_mb,
			disks = macs2_disks,
			time_hr = macs2_time_hr,
		}
	}

	File? peak_ppr2_ = if has_output_of_macs2_ppr2 then peak_ppr2
		else macs2_ppr2.npeak
	File? frip_macs2_qc_ppr2_ = if defined(macs2_ppr2__frip_qc) then macs2_ppr2__frip_qc
		else macs2_ppr2.frip_qc

	# IDR/overlap (post-peak-call analyses)

	# pipeline cannot resume from here
	# 	code will get too complicated for that

	# do IDR/overlap on all pairs of two replicates (i,j)
	# 	where i and j are zero-based indices and 0 <= i < j < num_rep
	Array[Pair[Int, Int]] pairs_ = cross(range(num_rep),range(num_rep))
	scatter( pair in pairs_ ) {
		Pair[Int, Int]? null_pair
		Pair[Int, Int]? pairs__ = if pair.left<pair.right then pair else null_pair
	}
	Array[Pair[Int, Int]] pairs = select_all(pairs__)

	scatter( pair in pairs ) {
		# pair.left = 0-based index of 1st replicate
		# pair.right = 0-based index of 2nd replicate
		Boolean has_input_of_overlap = 
			defined(peak_[pair.left]) && defined(peak_[pair.right]) && defined(peak_pooled_)
		Boolean has_output_of_overlap = pair.left<length(overlap__overlap_peak) &&
			defined(overlap__overlap_peak[pair.left])
		# for the case without frip_qc 
		#	where pipeline starts from peaks so tag-align is missing
 		Boolean has_output_of_overlap_frip_qc = pair.left<length(overlap__frip_qc) &&
			defined(overlap__frip_qc[pair.left])

		if ( has_input_of_overlap && !has_output_of_overlap &&
			pair.left<pair.right && !align_only ) { # only for repX-repY where X<Y
			# Naive overlap on every pair of true replicates
			call overlap { input :
				prefix = 'rep'+(pair.left+1)+"_rep"+(pair.right+1),
				peak1 = peak_[pair.left],
				peak2 = peak_[pair.right],
				peak_pooled = peak_pooled_,
				peak_type = peak_type,
				blacklist = blacklist_,
				chrsz = chrsz_,
				keep_irregular_chr_in_bfilt_peak = keep_irregular_chr_in_bfilt_peak,
				ta = ta_pooled_,
			}
		}
	}

	scatter( pair in pairs ) {
		# pair.left = 0-based index of 1st replicate
		# pair.right = 0-based index of 2nd replicate
		Boolean has_input_of_idr =
			defined(peak_[pair.left]) && defined(peak_[pair.right]) && defined(peak_pooled_)
		Boolean has_output_of_idr = pair.left<length(idr__idr_peak) &&
			defined(idr__idr_peak[pair.left])
		Boolean has_output_of_idr_frip_qc = pair.left<length(idr__frip_qc) &&
			defined(idr__frip_qc[pair.left])

		if ( has_input_of_idr && !has_output_of_idr &&
			pair.left<pair.right && !align_only && enable_idr) {
			# IDR on every pair of true replicates
			call idr { input :
				prefix = 'rep'+(pair.left+1)+"_rep"+(pair.right+1),
				peak1 = peak_[pair.left],
				peak2 = peak_[pair.right],
				peak_pooled = peak_pooled_,
				idr_thresh = idr_thresh,
				peak_type = peak_type,
				rank = idr_rank,
				blacklist = blacklist_,
				chrsz = chrsz_,
				keep_irregular_chr_in_bfilt_peak = keep_irregular_chr_in_bfilt_peak,
				ta = ta_pooled_,
			}
		}
	}

	scatter( i in range(length(pairs)) ) {
		# take raw peak if blacklist doesn't exist
		File? overlap_overlap_peak_ = 
			if has_output_of_overlap[i] && defined(blacklist_) then overlap__bfilt_overlap_peak[i]
			else if has_output_of_overlap[i] && !defined(blacklist_) then overlap__overlap_peak[i]
			else if defined(blacklist_) then overlap.bfilt_overlap_peak[i]
			else overlap.overlap_peak[i]
		File? overlap_frip_qc_ = if has_output_of_overlap_frip_qc[i] then overlap__frip_qc[i]
			else overlap.frip_qc[i]

		# take raw peak if blacklist doesn't exist
		File? idr_idr_peak_ = 
			if has_output_of_idr[i] && defined(blacklist_) then idr__bfilt_idr_peak[i]
			else if has_output_of_idr[i] && !defined(blacklist_) then idr__idr_peak[i]
			else if defined(blacklist_) then idr.bfilt_idr_peak[i]
			else idr.idr_peak[i]
		File? idr_idr_plot_ = if has_output_of_idr[i] then idr__idr_plot[i]
			else idr.idr_plot[i]
		File? idr_frip_qc_ = if has_output_of_idr_frip_qc[i] then idr__frip_qc[i]
			else idr.frip_qc[i]
	}

	# overlap on pseudo-replicates (pr1, pr2) for each true replicate
	scatter( i in range(num_rep) ) {
		Boolean has_input_of_overlap_pr =
			defined(peak_pr1_[i]) && defined(peak_pr2_[i]) && defined(peak_[i])
		Boolean has_output_of_overlap_pr = i<length(overlap_pr__overlap_peak) &&
			defined(overlap_pr__overlap_peak[i])
		Boolean has_output_of_overlap_pr_frip_qc = i<length(overlap_pr__frip_qc) &&
			defined(overlap_pr__frip_qc[i])

		if ( has_input_of_overlap_pr && !has_output_of_overlap_pr &&
			!align_only && !true_rep_only ) {
			call overlap as overlap_pr { input :
				prefix = "rep"+(i+1)+"-pr",
				peak1 = peak_pr1_[i],
				peak2 = peak_pr2_[i],
				peak_pooled = peak_[i],
				peak_type = peak_type,
				blacklist = blacklist_,
				chrsz = chrsz_,
				keep_irregular_chr_in_bfilt_peak = keep_irregular_chr_in_bfilt_peak,
				ta = ta_[i],
			}
		}
	}

	scatter( i in range(num_rep) ) {
		Boolean has_input_of_idr_pr = 
			defined(peak_pr1_[i]) && defined(peak_pr2_[i]) && defined(peak_[i])
		Boolean has_output_of_idr_pr = i<length(idr_pr__idr_peak) &&
			defined(idr_pr__idr_peak[i])
		Boolean has_output_of_idr_pr_frip_qc = i<length(idr_pr__frip_qc) &&
			defined(idr_pr__frip_qc[i])

		if ( has_input_of_idr_pr && !has_output_of_idr_pr &&
			!align_only && !true_rep_only && enable_idr ) {
			# IDR on pseduo replicates
			call idr as idr_pr { input :
				prefix = "rep"+(i+1)+"-pr",
				peak1 = peak_pr1_[i],
				peak2 = peak_pr2_[i],
				peak_pooled = peak_[i],
				idr_thresh = idr_thresh,
				peak_type = peak_type,
				rank = idr_rank,
				blacklist = blacklist_,
				chrsz = chrsz_,
				keep_irregular_chr_in_bfilt_peak = keep_irregular_chr_in_bfilt_peak,
				ta = ta_[i],
			}
		}
	}

	# gather outputs
	scatter( i in range(num_rep) ) {	
		File? overlap_pr_overlap_peak_ = 
			if has_output_of_overlap_pr[i] && defined(blacklist_) then overlap_pr__bfilt_overlap_peak[i]
			else if has_output_of_overlap_pr[i] && !defined(blacklist_) then overlap_pr__overlap_peak[i]
			else if defined(blacklist_) then overlap_pr.bfilt_overlap_peak[i]
			else overlap_pr.overlap_peak[i]
		File? overlap_pr_frip_qc_ = 
			if has_output_of_overlap_pr_frip_qc[i] then overlap_pr__frip_qc[i]
			else overlap_pr.frip_qc[i]

		File? idr_pr_idr_peak_ = 
			if has_output_of_idr_pr[i] && defined(blacklist_) then idr_pr__bfilt_idr_peak[i]
			else if has_output_of_idr_pr[i] && !defined(blacklist_) then idr_pr__idr_peak[i]
			else if defined(blacklist_) then idr_pr.bfilt_idr_peak[i]
			else idr_pr.idr_peak[i]
		File? idr_pr_idr_plot_ = 
			if has_output_of_idr_pr[i] then idr_pr__idr_plot[i]
			else idr_pr.idr_plot[i]
		File? idr_pr_frip_qc_ = 
			if has_output_of_idr_pr_frip_qc[i] then idr_pr__frip_qc[i]
			else idr_pr.frip_qc[i]
	}

	Boolean has_input_of_overlap_ppr =
		(has_output_of_macs2_ppr1 || defined(macs2_ppr1.npeak)) &&
		(has_output_of_macs2_ppr2 || defined(macs2_ppr2.npeak)) &&
		(has_output_of_macs2_pooled || defined(macs2_pooled.npeak))
	Boolean has_output_of_overlap_ppr = defined(overlap_ppr__overlap_peak)

	if ( has_input_of_overlap_ppr && !has_output_of_overlap_ppr &&
		!align_only && !true_rep_only && num_rep>1 ) {
		# Naive overlap on pooled pseudo replicates
		call overlap as overlap_ppr { input :
			prefix = "ppr",
			peak1 = peak_ppr1_,
			peak2 = peak_ppr2_,
			peak_pooled = peak_pooled_,
			peak_type = peak_type,
			blacklist = blacklist_,
			chrsz = chrsz_,
			keep_irregular_chr_in_bfilt_peak = keep_irregular_chr_in_bfilt_peak,
			ta = ta_pooled_,
		}
	}

	File? overlap_ppr_overlap_peak_ =
		if has_output_of_overlap_ppr && defined(blacklist_) then overlap_ppr__bfilt_overlap_peak
		else if has_output_of_overlap_ppr && !defined(blacklist_) then overlap_ppr__overlap_peak
		else if defined(blacklist_) then overlap_ppr.bfilt_overlap_peak
		else overlap_ppr.overlap_peak
	File? overlap_ppr_frip_qc_ =
		if defined(overlap_ppr__frip_qc) then overlap_ppr__frip_qc
		else overlap_ppr.frip_qc

	Boolean has_input_of_idr_ppr =
		(has_output_of_macs2_ppr1 || defined(macs2_ppr1.npeak)) &&
		(has_output_of_macs2_ppr2 || defined(macs2_ppr2.npeak)) &&
		(has_output_of_macs2_pooled || defined(macs2_pooled.npeak))
	Boolean has_output_of_idr_ppr = defined(idr_ppr__idr_peak)

	if ( has_input_of_idr_ppr && !has_output_of_idr_ppr &&
		!align_only && !true_rep_only && num_rep>1 ) {
		# IDR on pooled pseduo replicates
		call idr as idr_ppr { input :
			prefix = "ppr",
			peak1 = peak_ppr1_,
			peak2 = peak_ppr2_,
			peak_pooled = peak_pooled_,
			idr_thresh = idr_thresh,
			peak_type = peak_type,
			rank = idr_rank,
			blacklist = blacklist_,
			chrsz = chrsz_,
			keep_irregular_chr_in_bfilt_peak = keep_irregular_chr_in_bfilt_peak,
			ta = ta_pooled_,
		}
	}

	File? idr_ppr_idr_peak_ =
		if has_output_of_idr_ppr && defined(blacklist_) then idr_ppr__bfilt_idr_peak
		else if has_output_of_idr_ppr && !defined(blacklist_) then idr_ppr__idr_peak
		else if defined(blacklist_) then idr_ppr.bfilt_idr_peak
		else idr_ppr.idr_peak
	File? idr_ppr_idr_plot_ =
		if has_output_of_idr_ppr then idr_ppr__idr_plot
		else idr_ppr.idr_plot
	File? idr_ppr_frip_qc_ =
		if defined(idr_ppr__frip_qc) then idr_ppr__frip_qc
		else idr_ppr.frip_qc

	Boolean has_input_of_reproducibility_overlap =
		length(select_all(overlap_pr_overlap_peak_))==num_rep
	Boolean has_output_of_reproducibility_overlap =
		defined(reproducibility_overlap__reproducibility_qc)

	# reproducibility QC for overlap/IDR peaks
	if ( has_input_of_reproducibility_overlap && !has_output_of_reproducibility_overlap &&
		!align_only && !true_rep_only ) {
		# reproducibility QC for overlapping peaks
		call reproducibility as reproducibility_overlap { input :
			prefix = 'overlap',
			peaks = overlap_overlap_peak_,
			peaks_pr = overlap_pr_overlap_peak_,
			peak_ppr = overlap_ppr_overlap_peak_,
			peak_type = peak_type,
			chrsz = chrsz_,
			keep_irregular_chr_in_bfilt_peak = keep_irregular_chr_in_bfilt_peak,
		}
	}

	File? reproducibility_overlap_optimal_peak_ =
		if has_output_of_reproducibility_overlap then reproducibility_overlap__optimal_peak
		else reproducibility_overlap.optimal_peak
	File? reproducibility_overlap_reproducibility_qc_ =
		if has_output_of_reproducibility_overlap then reproducibility_overlap__reproducibility_qc
		else reproducibility_overlap.reproducibility_qc

	Boolean has_input_of_reproducibility_idr =
		length(select_all(idr_pr_idr_peak_))==num_rep
	Boolean has_output_of_reproducibility_idr =
		defined(reproducibility_idr__reproducibility_qc)

	if ( has_input_of_reproducibility_idr && !has_output_of_reproducibility_idr &&
		!align_only && !true_rep_only && enable_idr ) {
		# reproducibility QC for IDR peaks
		call reproducibility as reproducibility_idr { input :
			prefix = 'idr',
			peaks = idr_idr_peak_,
			peaks_pr = idr_pr_idr_peak_,
			peak_ppr = idr_ppr_idr_peak_,
			peak_type = peak_type,
			chrsz = chrsz_,
			keep_irregular_chr_in_bfilt_peak = keep_irregular_chr_in_bfilt_peak,
		}
	}

	File? reproducibility_idr_optimal_peak_ =
		if has_output_of_reproducibility_idr then reproducibility_idr__optimal_peak
		else reproducibility_idr.optimal_peak
	File? reproducibility_idr_reproducibility_qc_ =
		if has_output_of_reproducibility_idr then reproducibility_idr__reproducibility_qc
		else reproducibility_idr.reproducibility_qc

	# ATAqC
	scatter( i in range(num_rep) ) {

		Boolean has_input_of_ataqc = true # can run for any given inputs
		Boolean has_output_of_ataqc = i<length(ataqc__txt) &&
			defined(ataqc__txt[i])

		if ( has_input_of_ataqc && !has_output_of_ataqc &&
			!disable_ataqc ) {
			call ataqc { input :
				paired_end = paired_end_[i],
				read_len_log = read_len_log_[i],
				flagstat_qc = flagstat_qc_[i],
				bowtie2_log = align_log_[i],
				pbc_qc = pbc_qc_[i],
				dup_qc = dup_qc_[i],
				bam = bam_[i],
				nodup_flagstat_qc = nodup_flagstat_qc_[i],
				mito_dup_log = mito_dup_log_[i],
				nodup_bam = nodup_bam_[i],
				ta = ta_[i],
				peak = if defined(idr_pr_idr_peak_[i]) then idr_pr_idr_peak_[i]
					else reproducibility_overlap_optimal_peak_,
				idr_peak = reproducibility_idr_optimal_peak_, #idr_peaks_ataqc[i],
				overlap_peak= reproducibility_overlap_optimal_peak_, #overlap_peaks_ataqc[i],
				pval_bw = pval_bw_[i],
				ref_fa = ref_fa_,
				chrsz = chrsz_,
				tss = tss_,
				blacklist = blacklist_,
				dnase = dnase_,
				prom = prom_,
				enh = enh_,
				reg2map_bed = reg2map_bed_,
				reg2map = reg2map_,
				roadmap_meta = roadmap_meta_,
				mito_chr_name = mito_chr_name,

				mem_mb = ataqc_mem_mb,
				mem_java_mb = ataqc_mem_java_mb,
				time_hr = ataqc_time_hr,
				disks = ataqc_disks,
			}
		}
	}

	scatter( i in range(num_rep) ) {
		File? ataqc_txt_ = if has_output_of_ataqc[i] then ataqc__txt[i] else ataqc.txt[i]
		File? ataqc_html_ = if has_output_of_ataqc[i] then ataqc__html[i] else ataqc.html[i]
	}

	# Generate final QC report and JSON
	call qc_report { input :
		pipeline_ver = pipeline_ver,
		title = title,
		description = description,
		genome = basename(select_first([genome_tsv, ref_fa_, chrsz_, 'None'])),
		multimapping = multimapping,
		paired_ends = paired_end_,
		pipeline_type = pipeline_type,
		peak_caller = 'macs2',
		macs2_cap_num_peak = cap_num_peak,
		idr_thresh = idr_thresh,
		flagstat_qcs = flagstat_qc_,
		nodup_flagstat_qcs = nodup_flagstat_qc_,
		dup_qcs = dup_qc_,
		pbc_qcs = pbc_qc_,
		xcor_plots = xcor_plot_,
		xcor_scores = xcor_score_,

		frip_macs2_qcs = macs2_frip_qc_,
		frip_macs2_qcs_pr1 = macs2_pr1_frip_qc_,
		frip_macs2_qcs_pr2 = macs2_pr2_frip_qc_,

		frip_macs2_qc_pooled = frip_macs2_qc_pooled_,
		frip_macs2_qc_ppr1 = frip_macs2_qc_ppr1_,
		frip_macs2_qc_ppr2 = frip_macs2_qc_ppr2_,
		
		idr_plots = idr_idr_plot_,
		idr_plots_pr = idr_pr_idr_plot_,
		idr_plot_ppr = idr_ppr_idr_plot_,
		frip_idr_qcs = idr_frip_qc_,
		frip_idr_qcs_pr = idr_pr_frip_qc_,
		frip_idr_qc_ppr = idr_ppr_frip_qc_,
		frip_overlap_qcs = overlap_frip_qc_,
		frip_overlap_qcs_pr = overlap_pr_frip_qc_,
		frip_overlap_qc_ppr = overlap_ppr_frip_qc_,
		idr_reproducibility_qc = reproducibility_idr_reproducibility_qc_,
		overlap_reproducibility_qc = reproducibility_overlap_reproducibility_qc_,

		ataqc_txts = ataqc_txt_,
		ataqc_htmls = ataqc_html_,
	}

	output {
		File report = qc_report.report
		File qc_json = qc_report.qc_json
		Boolean qc_json_ref_match = qc_report.qc_json_ref_match
	}
}

# trim adapters and merge trimmed fastqs
task trim_adapter { 
	Array[File] fastqs_R1 		# [merge_id]
	Array[File] fastqs_R2

	String? adapter 	# adapter for all fastqs,
						#	this will override individual adapters in adapters_R1/R2
	Array[String] adapters_R1
	Array[String] adapters_R2
	Boolean paired_end
	Boolean auto_detect_adapter
	String cutadapt_param 
	# resource
	Int cpu
	Int mem_mb
	Int time_hr
	String disks

	# tmp vars
	File? null_f
	Array[Array[File]] tmp_fastqs = if paired_end then transpose([fastqs_R1, fastqs_R2])
				else transpose([fastqs_R1])
	Array[Array[String]] tmp_adapters = if paired_end then transpose([adapters_R1, adapters_R2])
				else transpose([adapters_R1])
	command {
		python $(which encode_trim_adapter.py) \
			${write_tsv(tmp_fastqs)} \
			${"--adapter " + adapter} \
			--adapters ${write_tsv(tmp_adapters)} \
			${if paired_end then "--paired-end" else ""} \
			${if auto_detect_adapter then "--auto-detect-adapter" else ""} \
			--cutadapt-param ' ${cutadapt_param}' \
			${"--nth " + cpu}
	}
	output {
		File trim_merged_fastq_R1 = glob("R1/*.fastq.gz")[0]
		File? trim_merged_fastq_R2 = if paired_end then glob("R2/*.fastq.gz")[0] else null_f
	}
	runtime {
		cpu : cpu
		memory : "${mem_mb} MB"
		time : time_hr
		disks : disks
	}
}

task bowtie2 {
	File idx_tar 		# reference bowtie2 index tar
	File? fastq_R1 		# [read_end_id]
	File? fastq_R2
	Boolean paired_end
	Int multimapping
	String bowtie2_param_se
	String bowtie2_param_pe
	Int cpu
	Int mem_mb
	Int time_hr
	String disks

	command {
		python $(which encode_bowtie2.py) \
			${idx_tar} \
			${fastq_R1} ${fastq_R2} \
			${if paired_end then "--paired-end" else ""} \
			${"--multimapping " + multimapping} \
			--bowtie2-param-se ' ${bowtie2_param_se}' \
			--bowtie2-param-pe ' ${bowtie2_param_pe}' \
			${"--nth " + cpu}
	}
	output {
		File bam = glob("*.bam")[0]
		File bai = glob("*.bai")[0]
		File align_log = glob("*.align.log")[0]
		File flagstat_qc = glob("*.flagstat.qc")[0]
		File read_len_log = glob("*.read_length.txt")[0] # read_len
	}
	runtime {
		cpu : cpu
		memory : "${mem_mb} MB"
		time : time_hr
		disks : disks
		preemptible: 0
	}
}

task filter {
	File bam
	Boolean paired_end
	Int multimapping
	String dup_marker 			# picard.jar MarkDuplicates (picard) or 
								# sambamba markdup (sambamba)
	Int mapq_thresh				# threshold for low MAPQ reads removal
	Boolean no_dup_removal 		# no dupe reads removal when filtering BAM
	String mito_chr_name
	Int cpu
	Int mem_mb
	Int time_hr
	String disks

	command {
		python $(which encode_filter.py) \
			${bam} \
			${if paired_end then "--paired-end" else ""} \
			${"--multimapping " + multimapping} \
			${"--dup-marker " + dup_marker} \
			${"--mapq-thresh " + mapq_thresh} \
			${if no_dup_removal then "--no-dup-removal" else ""} \
			${"--mito-chr-name " + mito_chr_name} \
			${"--nth " + cpu}
	}
	output {
		File nodup_bam = glob("*.bam")[0]
		File nodup_bai = glob("*.bai")[0]
		File flagstat_qc = glob("*.flagstat.qc")[0]
		File dup_qc = glob("*.dup.qc")[0]
		File pbc_qc = glob("*.pbc.qc")[0]
		File mito_dup_log = glob("*.mito_dup.txt")[0] # mito_dups, fract_dups_from_mito
	}
	runtime {
		cpu : cpu
		memory : "${mem_mb} MB"
		time : time_hr
		disks : disks
	}
}

task bam2ta {
	File bam
	Boolean paired_end
	Boolean disable_tn5_shift 	# no tn5 shifting (it's for dnase-seq)
	String regex_grep_v_ta   	# Perl-style regular expression pattern 
                        		# to remove matching reads from TAGALIGN
	String mito_chr_name 		# mito chromosome name
	Int subsample 				# number of reads to subsample TAGALIGN
								# this affects all downstream analysis
	Int cpu
	Int mem_mb
	Int time_hr
	String disks

	command {
		python $(which encode_bam2ta.py) \
			${bam} \
			${if paired_end then "--paired-end" else ""} \
			${if disable_tn5_shift then "--disable-tn5-shift" else ""} \
			${if regex_grep_v_ta!="" then "--regex-grep-v-ta '"+regex_grep_v_ta+"'" else ""} \
			${"--mito-chr-name " + mito_chr_name} \
			${"--subsample " + subsample} \
			${"--nth " + cpu}
	}
	output {
		File ta = glob("*.tagAlign.gz")[0]
	}
	runtime {
		cpu : cpu
		memory : "${mem_mb} MB"
		time : time_hr
		disks : disks
	}
}

task spr { # make two self pseudo replicates
	File ta
	Boolean paired_end

	Int mem_mb

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
		memory : "${mem_mb} MB"
		time : 1
		disks : "local-disk 50 HDD"
	}
}

task pool_ta {
	# input variables
	Array[File?] tas 	# TAG-ALIGNs to be merged

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
	File ta
	Boolean paired_end
	String mito_chr_name
	Int subsample  # number of reads to subsample TAGALIGN
				# this will be used for xcor only
				# will not affect any downstream analysis
	Int cpu
	Int mem_mb	
	Int time_hr
	String disks

	command {
		python $(which encode_xcor.py) \
			${ta} \
			${if paired_end then "--paired-end" else ""} \
			${"--mito-chr-name " + mito_chr_name} \
			${"--subsample " + subsample} \
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
		cpu : cpu
		memory : "${mem_mb} MB"
		time : time_hr
		disks : disks
	}
}

task count_signal_track {
	File ta 			# tag-align
	File chrsz			# 2-col chromosome sizes file

	command {
		python $(which encode_count_signal_track.py) \
			${ta} \
			${"--chrsz " + chrsz}
	}
	output {
		File pos_bw = glob("*.positive.bigwig")[0]
		File neg_bw = glob("*.negative.bigwig")[0]
	}
	runtime {
		cpu : 1
		memory : "8000 MB"
		time : 4
		disks : "local-disk 50 HDD"
	}
}

task macs2 {
	File ta
	String gensz		# Genome size (sum of entries in 2nd column of 
                        # chr. sizes file, or hs for human, ms for mouse)
	File chrsz			# 2-col chromosome sizes file
	Int cap_num_peak	# cap number of raw peaks called from MACS2
	Float pval_thresh  	# p.value threshold
	Int smooth_win 		# size of smoothing window
	File? blacklist 	# blacklist BED to filter raw peaks
	Boolean	keep_irregular_chr_in_bfilt_peak
	
	Int mem_mb
	Int time_hr
	String disks

	File? null_f

	command {
		python $(which encode_macs2_atac.py) \
			${ta} \
			${"--gensz "+ gensz} \
			${"--chrsz " + chrsz} \
			${"--cap-num-peak " + cap_num_peak} \
			${"--pval-thresh "+ pval_thresh} \
			${"--smooth-win "+ smooth_win} \
			${if keep_irregular_chr_in_bfilt_peak then "--keep-irregular-chr" else ""} \
			${"--blacklist "+ blacklist}
	}
	output {
		File npeak = glob("*[!.][!b][!f][!i][!l][!t].narrowPeak.gz")[0]
		File? bfilt_npeak = if defined(blacklist) then glob("*.bfilt.narrowPeak.gz")[0] else null_f
		File? bfilt_npeak_bb = if defined(blacklist) then glob("*.bfilt.narrowPeak.bb")[0] else null_f
		# hammock (.hammock.gz) and its index (.tbi)
		# Array[File] bfilt_npeak_hammock = glob("*.bfilt.narrowPeak.hammock.gz*")
		File? bfilt_npeak_hammock = if defined(blacklist) then glob("*.bfilt.narrowPeak.hammock.gz*")[0] else null_f
		File? bfilt_npeak_hammock_tbi = if defined(blacklist) then glob("*.bfilt.narrowPeak.hammock.gz*")[1] else null_f
		File frip_qc = glob("*.frip.qc")[0]
	}
	runtime {
		cpu : 1
		memory : "${mem_mb} MB"
		time : time_hr
		disks : disks
	}
}

task macs2_signal_track {
	File ta
	String gensz		# Genome size (sum of entries in 2nd column of 
                        # chr. sizes file, or hs for human, ms for mouse)
	File chrsz			# 2-col chromosome sizes file
	Float pval_thresh  	# p.value threshold
	Int smooth_win 		# size of smoothing window
	
	Int mem_mb
	Int time_hr
	String disks

	command {
		python $(which encode_macs2_signal_track_atac.py) \
			${ta} \
			${"--gensz "+ gensz} \
			${"--chrsz " + chrsz} \
			${"--pval-thresh "+ pval_thresh} \
			${"--smooth-win "+ smooth_win}
	}
	output {
		File pval_bw = glob("*.pval.signal.bigwig")[0]
		File fc_bw = glob("*.fc.signal.bigwig")[0]
	}
	runtime {
		cpu : 1
		memory : "${mem_mb} MB"
		time : time_hr
		disks : disks
	}
}

task idr {
	String prefix 		# prefix for IDR output file
	File peak1 			
	File peak2
	File peak_pooled
	Float idr_thresh
	File? blacklist 	# blacklist BED to filter raw peaks
	Boolean	keep_irregular_chr_in_bfilt_peak
	# parameters to compute FRiP
	File? ta			# to calculate FRiP
	File chrsz			# 2-col chromosome sizes file
	String peak_type
	String rank

	File? null_f

	command {
		python $(which encode_idr.py) \
			${peak1} ${peak2} ${peak_pooled} \
			${"--prefix " + prefix} \
			${"--idr-thresh " + idr_thresh} \
			${"--peak-type " + peak_type} \
			--idr-rank ${rank} \
			${"--chrsz " + chrsz} \
			${"--blacklist "+ blacklist} \
			${if keep_irregular_chr_in_bfilt_peak then "--keep-irregular-chr" else ""} \
			${"--ta " + ta}
	}
	output {
		File idr_peak = glob("*[!.][!b][!f][!i][!l][!t]."+peak_type+".gz")[0]
		File? bfilt_idr_peak = if defined(blacklist) then glob("*.bfilt."+peak_type+".gz")[0] else null_f
		File? bfilt_idr_peak_bb = if defined(blacklist) then glob("*.bfilt."+peak_type+".bb")[0] else null_f
		# Array[File] bfilt_idr_peak_hammock = glob("*.bfilt."+peak_type+".hammock.gz*")
		File? bfilt_idr_peak_hammock = if defined(blacklist) then glob("*.bfilt."+peak_type+".hammock.gz*")[0] else null_f
		File? bfilt_idr_peak_hammock_tbi = if defined(blacklist) then glob("*.bfilt."+peak_type+".hammock.gz*")[1] else null_f
		File idr_plot = glob("*.txt.png")[0]
		File idr_unthresholded_peak = glob("*.txt.gz")[0]
		File idr_log = glob("*.idr*.log")[0]
		File? frip_qc = if defined(ta) then glob("*.frip.qc")[0] else null_f
	}
	runtime {
		cpu : 1
		memory : "8000 MB"
		time : 1
		disks : "local-disk 50 HDD"
	}
}

task overlap {
	String prefix 		# prefix for IDR output file
	File peak1
	File peak2
	File peak_pooled
	File? blacklist 	# blacklist BED to filter raw peaks
	Boolean	keep_irregular_chr_in_bfilt_peak
	File? ta		# to calculate FRiP
	File chrsz			# 2-col chromosome sizes file
	String peak_type

	File? null_f

	command {
		python $(which encode_naive_overlap.py) \
			${peak1} ${peak2} ${peak_pooled} \
			${"--prefix " + prefix} \
			${"--peak-type " + peak_type} \
			${"--chrsz " + chrsz} \
			${"--blacklist "+ blacklist} \
			--nonamecheck \
			${if keep_irregular_chr_in_bfilt_peak then "--keep-irregular-chr" else ""} \
			${"--ta " + ta}
	}
	output {
		File overlap_peak = glob("*[!.][!b][!f][!i][!l][!t]."+peak_type+".gz")[0]
		File? bfilt_overlap_peak = if defined(blacklist) then glob("*.bfilt."+peak_type+".gz")[0] else null_f
		File? bfilt_overlap_peak_bb = if defined(blacklist) then glob("*.bfilt."+peak_type+".bb")[0] else null_f
		# Array[File] bfilt_overlap_peak_hammock = glob("*.bfilt."+peak_type+".hammock.gz*")
		File? bfilt_overlap_peak_hammock = if defined(blacklist) then glob("*.bfilt."+peak_type+".hammock.gz*")[0] else null_f
		File? bfilt_overlap_peak_hammock_tbi = if defined(blacklist) then glob("*.bfilt."+peak_type+".hammock.gz*")[1] else null_f
		File? frip_qc = if defined(ta) then glob("*.frip.qc")[0] else null_f
	}
	runtime {
		cpu : 1
		memory : "4000 MB"
		time : 1
		disks : "local-disk 50 HDD"
	}
}

task reproducibility {
	String prefix
	Array[File?] peaks # peak files from pair of true replicates
						# in a sorted order. for example of 4 replicates,
						# 1,2 1,3 1,4 2,3 2,4 3,4.
                        # x,y means peak file from rep-x vs rep-y
	Array[File?] peaks_pr	# peak files from pseudo replicates
	File? peak_ppr			# Peak file from pooled pseudo replicate.
	String peak_type
	File chrsz			# 2-col chromosome sizes file
	Boolean	keep_irregular_chr_in_bfilt_peak

	command {
		python $(which encode_reproducibility_qc.py) \
			${sep=' ' peaks} \
			--peaks-pr ${sep=' ' peaks_pr} \
			${"--peak-ppr "+ peak_ppr} \
			--prefix ${prefix} \
			${"--peak-type " + peak_type} \
			${if keep_irregular_chr_in_bfilt_peak then "--keep-irregular-chr" else ""} \
			${"--chrsz " + chrsz}
	}
	output {
		File optimal_peak = glob("optimal_peak.*.gz")[0]
		File conservative_peak = glob("conservative_peak.*.gz")[0]
		File optimal_peak_bb = glob("optimal_peak.*.bb")[0]
		File conservative_peak_bb = glob("conservative_peak.*.bb")[0]
		# Array[File] optimal_peak_hammock = glob("optimal_peak.*.hammock.gz*")
		# Array[File] conservative_peak_hammock = glob("conservative_peak.*.hammock.gz*")
		File optimal_peak_hammock = glob("optimal_peak.*.hammock.gz*")[0]
		File optimal_peak_hammock_tbi = glob("optimal_peak.*.hammock.gz*")[1]
		File conservative_peak_hammock = glob("conservative_peak.*.hammock.gz*")[0]
		File conservative_peak_hammock_tbi = glob("conservative_peak.*.hammock.gz*")[1]
		File reproducibility_qc = glob("*reproducibility.qc")[0]
	}
	runtime {
		cpu : 1
		memory : "4000 MB"
		time : 1
		disks : "local-disk 50 HDD"
	}
}

# annotation-based analysis
task ataqc {
	Boolean paired_end
	File? read_len_log
	File? flagstat_qc
	File? bowtie2_log
	File? bam
	File? nodup_flagstat_qc
	File? mito_dup_log
	File? dup_qc
	File? pbc_qc
	File? nodup_bam
	File? ta
	File? peak
	File? idr_peak 
	File? overlap_peak
	File? pval_bw
	# from genome database
	File? ref_fa
	File? chrsz
	File? tss
	File? blacklist
	File? dnase
	File? prom
	File? enh
	File? reg2map_bed
	File? reg2map
	File? roadmap_meta
	String mito_chr_name

	Int mem_mb
	Int mem_java_mb
	Int time_hr
	String disks

	command {
		export _JAVA_OPTIONS="-Xms256M -Xmx${mem_java_mb}M -XX:ParallelGCThreads=1 $_JAVA_OPTIONS"

		python $(which encode_ataqc.py) \
			${if paired_end then "--paired-end" else ""} \
			${"--read-len-log " + read_len_log} \
			${"--flagstat-log " + flagstat_qc} \
			${"--bowtie2-log " + bowtie2_log} \
			${"--bam " + bam} \
			${"--nodup-flagstat-log " + nodup_flagstat_qc} \
			${"--mito-dup-log " + mito_dup_log} \
			${"--dup-log " + dup_qc} \
			${"--pbc-log " + pbc_qc} \
			${"--nodup-bam " + nodup_bam} \
			${"--ta " + ta} \
			${"--bigwig " + pval_bw} \
			${"--peak " + peak} \
			${"--idr-peak " + idr_peak} \
			${"--overlap-peak " + overlap_peak} \
			${"--ref-fa " + ref_fa} \
			${"--blacklist " + blacklist} \
			${"--chrsz " + chrsz} \
			${"--dnase " + dnase} \
			${"--tss " + tss} \
			${"--prom " + prom} \
			${"--enh " + enh} \
			${"--reg2map-bed " + reg2map_bed} \
			${"--reg2map " + reg2map} \
			${"--roadmap-meta " + roadmap_meta} \
			${"--mito-chr-name " + mito_chr_name}

	}
	output {
		File html = glob("*_qc.html")[0]
		File txt = glob("*_qc.txt")[0]
	}
	runtime {
		cpu : 1
		memory : "${mem_mb} MB"
		time : time_hr
		disks : disks
	}
}

# gather all outputs and generate 
# - qc.html		: organized final HTML report
# - qc.json		: all QCs
task qc_report {
	String pipeline_ver
 	String title
	String description
	String? genome	
	# workflow params
	Int multimapping
	Array[Boolean?] paired_ends
	String pipeline_type
	String peak_caller
	Int? macs2_cap_num_peak
	Int? spp_cap_num_peak
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
	Array[File?] frip_macs2_qcs
	Array[File?] frip_macs2_qcs_pr1
	Array[File?] frip_macs2_qcs_pr2
	File? frip_macs2_qc_pooled
	File? frip_macs2_qc_ppr1 
	File? frip_macs2_qc_ppr2 
	Array[File?] frip_idr_qcs
	Array[File?] frip_idr_qcs_pr
	File? frip_idr_qc_ppr 
	Array[File?] frip_overlap_qcs
	Array[File?] frip_overlap_qcs_pr
	File? frip_overlap_qc_ppr
	File? idr_reproducibility_qc
	File? overlap_reproducibility_qc
	Array[File?] ataqc_txts
	Array[File?] ataqc_htmls

	File? qc_json_ref

	command {
		python $(which encode_qc_report.py) \
			${"--pipeline-ver " + pipeline_ver} \
			${"--title '" + sub(title,"'","_") + "'"} \
			${"--desc '" + sub(description,"'","_") + "'"} \
			${"--genome " + genome} \
			${"--multimapping " + multimapping} \
			--paired-ends ${sep=" " paired_ends} \
			--pipeline-type ${pipeline_type} \
			--peak-caller ${peak_caller} \
			${"--macs2-cap-num-peak " + macs2_cap_num_peak} \
			${"--spp-cap-num-peak " + spp_cap_num_peak} \
			--idr-thresh ${idr_thresh} \
			--flagstat-qcs ${sep="_:_" flagstat_qcs} \
			--nodup-flagstat-qcs ${sep="_:_" nodup_flagstat_qcs} \
			--dup-qcs ${sep="_:_" dup_qcs} \
			--pbc-qcs ${sep="_:_" pbc_qcs} \
			--xcor-plots ${sep="_:_" xcor_plots} \
			--xcor-scores ${sep="_:_" xcor_scores} \
			--idr-plots ${sep="_:_" idr_plots} \
			--idr-plots-pr ${sep="_:_" idr_plots_pr} \
			${"--idr-plot-ppr " + idr_plot_ppr} \
			--frip-macs2-qcs ${sep="_:_" frip_macs2_qcs} \
			--frip-macs2-qcs-pr1 ${sep="_:_" frip_macs2_qcs_pr1} \
			--frip-macs2-qcs-pr2 ${sep="_:_" frip_macs2_qcs_pr2} \
			${"--frip-macs2-qc-pooled " + frip_macs2_qc_pooled} \
			${"--frip-macs2-qc-ppr1 " + frip_macs2_qc_ppr1} \
			${"--frip-macs2-qc-ppr2 " + frip_macs2_qc_ppr2} \
			--frip-idr-qcs ${sep="_:_" frip_idr_qcs} \
			--frip-idr-qcs-pr ${sep="_:_" frip_idr_qcs_pr} \
			${"--frip-idr-qc-ppr " + frip_idr_qc_ppr} \
			--frip-overlap-qcs ${sep="_:_" frip_overlap_qcs} \
			--frip-overlap-qcs-pr ${sep="_:_" frip_overlap_qcs_pr} \
			${"--frip-overlap-qc-ppr " + frip_overlap_qc_ppr} \
			${"--idr-reproducibility-qc " + idr_reproducibility_qc} \
			${"--overlap-reproducibility-qc " + overlap_reproducibility_qc} \
			--ataqc-txts ${sep="_:_" ataqc_txts} \
			--ataqc-htmls ${sep="_:_" ataqc_htmls} \
			--out-qc-html qc.html \
			--out-qc-json qc.json \
			${"--qc-json-ref " + qc_json_ref}
	}
	output {
		File report = glob('*qc.html')[0]
		File qc_json = glob('*qc.json')[0]
		Boolean qc_json_ref_match = read_string("qc_json_ref_match.txt")=="True"
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

	String? null_s
	command <<<
		# create empty files for all entries
		touch ref_fa bowtie2_idx_tar chrsz gensz blacklist
		touch tss tss_enrich # for backward compatibility
		touch dnase prom enh reg2map reg2map_bed roadmap_meta

		python <<CODE
		import os
		with open("${genome_tsv}",'r') as fp:
			for line in fp:
				arr = line.strip('\n').split('\t')
				if arr:
					key, val = arr
					with open(key,'w') as fp2:
						fp2.write(val)
		CODE
	>>>
	output {
		String? ref_fa = if size('ref_fa')==0 then null_s else read_string('ref_fa')
		String? bowtie2_idx_tar = if size('bowtie2_idx_tar')==0 then null_s else read_string('bowtie2_idx_tar')
		String? chrsz = if size('chrsz')==0 then null_s else read_string('chrsz')
		String? gensz = if size('gensz')==0 then null_s else read_string('gensz')
		String? blacklist = if size('blacklist')==0 then null_s else read_string('blacklist')
		String? tss = if size('tss')!=0 then read_string('tss')
			else if size('tss_enrich')!=0 then read_string('tss_enrich') else null_s
		String? dnase = if size('dnase')==0 then null_s else read_string('dnase')
		String? prom = if size('prom')==0 then null_s else read_string('prom')
		String? enh = if size('enh')==0 then null_s else read_string('enh')
		String? reg2map = if size('reg2map')==0 then null_s else read_string('reg2map')
		String? reg2map_bed = if size('reg2map_bed')==0 then null_s else read_string('reg2map_bed')
		String? roadmap_meta = if size('roadmap_meta')==0 then null_s else read_string('roadmap_meta')
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
