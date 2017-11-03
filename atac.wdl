task detect_adapter {
	Array[Array[File]] fastqs # fastqs[tech_rep_id][pair_id]
	Boolean? auto_detect_adapter
}

task trim_adapter {
	Array[Array[File]] fastqs # fastqs[tech_rep_id][pair_id]
	Array[Array[String]] adapters # fastqs[tech_rep_id][pair_id]
}

task trim { # detect/trim adapter
	Array[Array[File]] fastqs # fastqs[tech_rep_id][pair_id]
	# optional
	Boolean? auto_detect_adapter
	Int? min_trim_len
	Float? err_rate
	# do not touch
	Array[Array[String]] adapters # adapters[tech_rep_id][pair_id]
	Boolean? paired_end
	# resource
	Int? cpu
	Int? mem_mb
	Int? time_hr
	String? queue

	command {
		python $(which encode_dcc_trim.py) \
			${write_tsv(fastqs)} \
			${if select_first([auto_detect_adapter,false])
				then "--auto-detect-adapter" else ""} \
			${"--min-trim-len " + min_trim_len} \
			${"--err-rate " + err_rate} \
			${"--adapters " + write_tsv(adapters)} \
			${if select_first([paired_end,false])
				then "--paired-end" else ""} \
			${"--nth " + cpu}
	}
	output {
		# trimmed_fastqs[tech_rep_id][pair_id]
		Array[Array[File]] trimmed_fastqs = read_tsv("out.tsv")
	}
	runtime {
		cpu : "${select_first([cpu,1])}"
		memory : "${select_first([mem_mb,'100'])} MB"
		time : "${select_first([time_hr,1])}"
		queue : queue
	}
}

task merge { # merge fastqs
	Array[Array[File]] fastqs # fastqs[tech_rep_id][pair_id]
	# resource
	Int? cpu
	Int? mem_mb
	Int? time_hr
	String? queue

	command {
		python $(which encode_dcc_merge_fastq.py) \
			${write_tsv(fastqs)}
	}
	output {
		# merged_fastqs[pair_id]
		Array[File] merged_fastqs = read_lines("out.lines")
	}
	runtime {
		cpu : "${select_first([cpu,1])}"
		memory : "${select_first([mem_mb,'100'])} MB"
		time : "${select_first([time_hr,1])}"
		queue : queue
	}
}

task bowtie2 {
	File idx_tar # reference bowtie2 index tar
	Array[File] fastqs # fastqs[pair_id]
	# optional
	String? score_min 	# min acceptable alignment score func
						# w.r.t read length
	# do not touch
	Boolean? paired_end
	Int? multimapping						
	# resource
	Int? cpu
	Int? mem_mb
	Int? time_hr
	String? queue

	command {
		python $(which encode_dcc_bowtie2.py) \
			${idx_tar} \
			${sep=' ' fastqs} \
			${"--score-min " + score_min} \
			${if select_first([paired_end,false])
				then "--paired-end" else ""} \
			${"--multimapping " + multimapping} \
			${"--nth " + cpu}
	}
	output {
		File bam = glob("*.bam")[0]
		File bai = glob("*.bai")[0]
		File align_log = glob("*.align.log")[0]
		File flagstat_qc = glob("*.flagstat.qc")[0]
	}
	runtime {
		cpu : "${select_first([cpu,1])}"
		memory : "${select_first([mem_mb,'100'])} MB"
		time : "${select_first([time_hr,1])}"
		queue : queue
	}
}

task filter {
	File bam
	# optional
	String? dup_marker
	Int? mapq_thresh
	Boolean? no_dup_removal
	# do not touch
	Boolean? paired_end
	Int? multimapping
	# resource
	Int? cpu
	Int? mem_mb
	Int? time_hr
	String? queue

	command {
		python $(which encode_dcc_filter.py) \
			${bam} \
			${"--dup-marker " + dup_marker} \
			${"--mapq-thresh " + mapq_thresh} \
			${if select_first([no_dup_removal,false])
				then "--no-dup-removal" else ""} \
			${if select_first([paired_end,false])
				then "--paired-end" else ""} \
			${"--multimapping " + multimapping} \
			${"--nth " + cpu}			
	}
	output {
		File nodup_bam = glob("*.bam")[0]
		File nodup_bai = glob("*.bai")[0]
		File flagstat_qc = glob("*.flagstat.qc")[0]
		File dup_qc = glob("*.dup.qc")[0]
		File pbc_qc = glob("*.pbc.qc")[0]
	}

	runtime {
		cpu : "${select_first([cpu,1])}"
		memory : "${select_first([mem_mb,'100'])} MB"
		time : "${select_first([time_hr,1])}"
		queue : queue
	}
}

task bam2ta {
	File bam
	# optional
	Boolean? disable_tn5_shift # for dnase-seq
	String? regex_grep_v_ta # 
	Int? subsample # number of reads to subsample
	# do not touch
	Boolean? paired_end
	# resource
	Int? cpu
	Int? mem_mb
	Int? time_hr
	String? queue

	command {
		python $(which encode_dcc_bam2ta.py) \
			${bam} \
			${if select_first([disable_tn5_shift,false])
				then "--disable-tn5-shift" else ""} \
			${"--regex-grep-v-ta " +"'"+regex_grep_v_ta+"'"} \
			${"--subsample " + subsample} \
			${if select_first([paired_end,false])
				then "--paired-end" else ""} \
			${"--nth " + cpu}
	}
	output {
		File ta = glob("*.tagAlign.gz")[0]
	}

	runtime {
		cpu : "${select_first([cpu,1])}"
		memory : "${select_first([mem_mb,'100'])} MB"
		time : "${select_first([time_hr,1])}"
		queue : queue
	}
}

task spr { # make two self pseudo replicates
	File ta
	# do not touch
	Boolean? paired_end
	# resource
	Int? cpu
	Int? mem_mb
	Int? time_hr
	String? queue

	command {
		python $(which encode_dcc_spr.py) \
			${ta} \
			${if select_first([paired_end,false])
				then "--paired-end" else ""}
	}
	output {
		File ta_pr1 = glob("*.pr1.tagAlign.gz")[0]
		File ta_pr2 = glob("*.pr2.tagAlign.gz")[0]
	}
	runtime {
		cpu : "${select_first([cpu,1])}"
		memory : "${select_first([mem_mb,'100'])} MB"
		time : "${select_first([time_hr,1])}"
		queue : queue
	}
}

task pool_ta {
	Array[File] tas
	# resource
	Int? cpu
	Int? mem_mb
	Int? time_hr
	String? queue

	command {
		python $(which encode_dcc_pool_ta.py) \
			${sep=' ' tas}
	}
	output {
		File ta_pooled = glob("*.tagAlign.gz")[0]
	}
	runtime {
		cpu : "${select_first([cpu,1])}"
		memory : "${select_first([mem_mb,'100'])} MB"
		time : "${select_first([time_hr,1])}"
		queue : queue
	}
}

task xcor {
	File ta
	# optional
	Int? subsample # number of reads to subsample
	# do not touch
	Boolean? paired_end
	# resource
	Int? cpu
	Int? mem_mb
	Int? time_hr
	String? queue

	command {
		python $(which encode_dcc_xcor.py) \
			${ta} \
			${"--subsample " + subsample} \
			--speak=0 \
			${if select_first([paired_end,false])
				then "--paired-end" else ""} \
			${"--nth " + cpu}
	}
	output {
		File plot = glob("*.cc.plot.pdf")[0]
		File score = glob("*.cc.qc")[0]
	}
	runtime {
		cpu : "${select_first([cpu,1])}"
		memory : "${select_first([mem_mb,'100'])} MB"
		time : "${select_first([time_hr,1])}"
		queue : queue
	}
}

task macs2 {
	File ta
	File chrsz
	# optional
	Int? cap_num_peak # cap number of peaks with top scores
	String? gensz
	Float? pval_thresh
    Int? smooth_win
    # do not touch
    Boolean? make_signal
	# resource
	Int? cpu
	Int? mem_mb
	Int? time_hr
	String? queue

	command {
		python $(which encode_dcc_macs2.py) \
			${ta} \
			${"--chrsz " + chrsz} \
			${"--cap-num-peak " + cap_num_peak} \
			${"--gensz "+ gensz} \
			${"--p-val-thresh "+ pval_thresh} \
			${"--smooth-win "+ smooth_win} \
			${if select_first([make_signal,false])
				then "--make-signal" else ""}
	}
	output {
		File npeak = glob("*.narrowPeak.gz")[0]
		File sig_pval = glob("*.pval.signal.bigwig")[0]
		File sig_fc = glob("*.fc.signal.bigwig")[0]
	}
	runtime {
		cpu : "${select_first([cpu,1])}"
		memory : "${select_first([mem_mb,'100'])} MB"
		time : "${select_first([time_hr,1])}"
		queue : queue
	}
}

task reproducibility { # IDR and naive overlap
	Array[File] peaks
	Array[File] peaks_pr1
	Array[File] peaks_pr2
	File peak_pooled
	File peak_ppr1
	File peak_ppr2
	String method # idr, overlap
	Float? idr_thresh
	String? idr_rank

	command {
		python $(which encode_dcc_macs2.py) \
			${ta} \
			--peaks ${sep=' ' peaks} \
			--peak-pr1 ${sep=' ' peaks_pr1} \
			--peak-pr2 ${sep=' ' peaks_pr2} \
			${"--peak-pooled "+ peak_pooled} \
			${"--peak-ppr1 "+ peak_ppr1} \
			${"--peak-ppr2 "+ peak_ppr2} \
			${"--method "+ method} \
			${"--idr-thresh "+ idr_tresh} \
			${"--idr-rank "+ idr_rank}
	}
	output {
		# IDR/overlapping peaks
		Map[String,File] out_peaks_tr # true 
		Array[File] out_peak_pr  # pseudo replicated
		File out_peak_ppr # pooled pseudo replicate
		File reproducibility_qc 
	}
	runtime {
		cpu : "${select_first([cpu,1])}"
		memory : "${select_first([mem_mb,'100'])} MB"
		time : "${select_first([time_hr,1])}"
		queue : queue
	}
}

task blacklist_filter {
	File peak
	File blacklist
	Boolean? keep_irregular_chr

	command {
		python $(which encode_dcc_blacklist_filter.py) \
			${peak} \
			${blacklist} \
			${if select_first([keep_irregular_chr,false])
				then "--keep-irregular-chr" else ""}
	}
	output {
		File filtered_peak = glob('*.gz')[0]
	}
	runtime {
		cpu : "${select_first([cpu,1])}"
		memory : "${select_first([mem_mb,'100'])} MB"
		time : "${select_first([time_hr,1])}"
		queue : queue
	}
}

task frip { # IDR and naive overlap
	File peak
	File ta

	command {
		python $(which encode_dcc_frip.py) \
			${peak} \
			${ta}
	}
	output {
		File frip_qc = glob('*.frip.qc')[0]
	}
	runtime {
		cpu : "${select_first([cpu,1])}"
		memory : "${select_first([mem_mb,'100'])} MB"
		time : "${select_first([time_hr,1])}"
		queue : queue
	}
}

workflow atac {
	# mandatory parameters
	Boolean paired_end 	# endedness of sample
	File genome_tsv 	# reference genome data TSV file
	Map[String,String] genome = read_map(genome_tsv)

	# optional but important
	Int? multimapping 	# used for bowtie2 and filter
	Int? cpu 			# total number of threads to process this sample
						# must be multiples of num_rep
	# optional
	String? name 		# name of sample
	String? description # description for sample
	String? accession_id # accession ID of sample on ENCODE portal

	# mandatory input files
	Int num_rep 		# number of replicates
	String input_type 	# choices = ['fastq','bam','nodup_bam','ta']
	Array[Array[Array[String]]] fastqs 
						# fastqs[bio_rep_id][tech_rep_id][pair_id]
	Array[Array[Array[String]]] adapters 
						# adapters[bio_rep_id][tech_rep_id][pair_id]
	Array[String] bams 			# if starting from bams
	Array[String] nodup_bams 	# if starting from filtered bams
	Array[String] tas 			# if starting from tag-aligns

	if (input_type=='fastq') {
		# trim adapters 
		scatter(i in range(num_rep)) {
			call trim {
				input:
					fastqs = fastqs[i],
					adapters = adapters[i],
					paired_end = paired_end,
					cpu = cpu/num_rep,
			}
		}
		# merge fastqs from technical replicates
		scatter(i in range(num_rep)) {
			call merge {
				input: 
					fastqs = trim.trimmed_fastqs[i]
			}
		}
		# align trimmed/merged fastqs with bowtie2
		scatter(i in range(num_rep)) {
			call bowtie2 {
				input:
					fastqs = merge.merged_fastqs[i],
					idx_tar = genome["bowtie2_idx_tar"],
					paired_end = paired_end,
					multimapping = multimapping,
					cpu = cpu/num_rep,
			}
		}
	}

	if (input_type=='fastq' || input_type=='bam') {
		# filter/dedup bam
		scatter(i in range(num_rep)) {
			call filter {
				input:
					bam = select_first([bowtie2.bam,bams])[i],
					paired_end = paired_end,
					multimapping = multimapping,
					cpu = cpu/num_rep,
			}
		}
	}

	if (input_type=='fastq' || input_type=='bam' ||
		input_type=='nodup_bam') {
		# convert bam to tagalign and subsample it if necessary
		scatter(i in range(num_rep)) {
			call bam2ta {
				input:
					bam = select_first([filter.nodup_bam,nodup_bams])[i],
					paired_end = select_first([paired_end,true]),
					cpu = cpu/num_rep,
			}
		}
	}

	# subsample tagalign (non-mito) and cross-correlation analysis
	scatter(i in range(num_rep)) {
		call xcor {
			input:
				ta = select_first([bam2ta.ta,tas])[i],
				paired_end = select_first([paired_end,true]),
		}
	}

	# make two self pseudo replicates per true replicate
	scatter(i in range(num_rep)) {
		call spr {
			input:
				ta = select_first([bam2ta.ta,tas])[i],
				paired_end = select_first([paired_end,true]),
		}
	}

	# pool tagaligns from true/pseudo replicates
	call pool_ta {
		input :
			tas = select_first([bam2ta.ta,tas]),
	}
	call pool_ta as pool_ta_pr1 {
		input :
			tas = spr.ta_pr1,
	}
	call pool_ta as pool_ta_pr2 {
		input :
			tas = spr.ta_pr2,
	}

	# call peaks on tagalign
	scatter(i in range(num_rep)) {
		call macs2 {
			input:
				ta = select_first([bam2ta.ta,tas])[i],
				chrsz = genome["chrsz"],
				make_signal = true,
		}
	}
	# call peaks on 1st pseudo replicated tagalign 
	scatter(ta in spr.ta_pr1) {
		call macs2 as macs2_pr1 {
			input:
				ta = ta,
				chrsz = genome["chrsz"],
		}
	}
	# call peaks on 2nd pseudo replicated tagalign 
	scatter(ta in spr.ta_pr2) {
		call macs2 as macs2_pr2 {
			input:
				ta = ta,
				chrsz = genome["chrsz"],
		}
	}
	# call peaks on pooled replicate
	call macs2 as macs2_pooled {
		input:
			ta = pool_ta.ta_pooled,
			chrsz = genome["chrsz"],
			make_signal = true,
	}
	# call peaks on 1st pooled pseudo replicates
	call macs2 as macs2_ppr1 {
		input:
			ta = pool_ta_pr1.ta_pooled,
			chrsz = genome["chrsz"],
	}
	# call peaks on 2nd pooled pseudo replicates
	call macs2 as macs2_ppr2 {
		input:
			ta = pool_ta_pr2.ta_pooled,
			chrsz = genome["chrsz"],
	}

	# reproducibility QC and FRiP on called peaks
	## IDR
	call reproducibility as idr {
		input:
			peaks = macs2.
	}
	## naive overlap
	call reproducibility as overlap {
	}

	# filter out peaks overlapping blacklist
	call blacklist_filter as blacklist_filt_idr {
		input:
			blacklist = genome["blacklist"],
	}
}