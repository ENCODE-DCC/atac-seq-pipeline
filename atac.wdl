#### genomic tasks
task trim_adapter { # detect/trim adapter
	# mandatory
	Boolean auto_detect_adapter
	# optional
	Int? min_trim_len
	Float? err_rate
	# do not touch
	Array[Array[File]] fastqs # fastqs[tech_rep_id][pair_id]
	Array[Array[String]] adapters # adapters[tech_rep_id][pair_id]
	Boolean? paired_end
	# resource
	Int? cpu
	Int? mem_mb
	Int? time_hr
	String? queue

	command {
		python $(which encode_dcc_trim_adapter.py) \
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

task merge_fastq { # merge fastqs
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
		Array[File] merged_fastqs = read_lines("out.txt")
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

task idr {
	String? prefix
	File peak1
	File peak2
	File peak_pooled
	Float? idr_thresh

	command {
		python $(which encode_dcc_idr.py) \
			${peak1} ${peak2} ${peak_pooled} \
			${"--prefix " + prefix} \
			${"--idr-thresh " + idr_thresh} \
			--idr-rank signal.value
	}
	output {
		File idr_peak = glob("*peak.gz")[0]
		File idr_plot = glob("*.txt.png")[0]
		File idr_unthresholded_peak = glob("*.txt.gz")[0]
		File idr_log = glob("*.log")[0]
	}
	runtime {
		cpu : "${select_first([cpu,1])}"
		memory : "${select_first([mem_mb,'100'])} MB"
		time : "${select_first([time_hr,1])}"
		queue : queue
	}
}

task overlap {
	String? prefix
	File peak1
	File peak2
	File peak_pooled

	command {
		python $(which encode_dcc_naive_overlap.py) \
			${peak1} ${peak2} ${peak_pooled} \
			${"--prefix " + prefix}
	}
	output {
		File overlap_peak = glob("*peak.gz")[0]
	}
	runtime {
		cpu : "${select_first([cpu,1])}"
		memory : "${select_first([mem_mb,'100'])} MB"
		time : "${select_first([time_hr,1])}"
		queue : queue
	}
}

task reproducibility {
	Array[File]? peaks
	Array[File] peaks_pr
	File? peak_ppr

	command {
		python $(which encode_dcc_reproducibility_qc.py) \			
			${sep=' ' peaks} \
			--peaks-pr ${sep=' ' peaks_pr} \
			${"--peak-ppr "+ peak_ppr}
	}
	output {
		File reproducibility_qc = 
			glob("*.reproducibility.qc")[0]
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

#### workflow system tasks
task inputs { # determine input type and number of replicates
	Array[Array[Array[String]]] fastqs 
	Array[String] bams
	Array[String] nodup_bams
	Array[String] tas

	command <<<
		python <<CODE
		name = ['fastq','bam','nodup_bam','ta']
		arr = [${length(fastqs)},${length(bams)},
		       ${length(nodup_bams)},${length(tas)}]
		num_rep = max(arr)
		type = name[arr.index(num_rep)]
		with open('num_rep.txt','w') as fp:
		    fp.write(str(num_rep)) 
		with open('type.txt','w') as fp:
		    fp.write(type)		    
		CODE
	>>>
	output {
		String type = read_string("type.txt")
		Int num_rep = read_int("num_rep.txt")
	}
}

task pair_gen { # returns every pair of true replicate
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

#### workflow body
workflow atac {
	# mandatory parameters
	Boolean paired_end 	# endedness of sample
	File genome_tsv 	# reference genome data TSV file
	Map[String,String] genome = read_map(genome_tsv)

	# optional but important
	Boolean? enable_idr	# enable IDR for called peaks
	Boolean? true_rep_only # disable all analysis for pseudo replicates
	Int? multimapping 	# used for bowtie2 and filter
	Int? cpu 			# total number of threads to process this sample
						# must be multiples of num_rep
	# optional
	String? name 		# name of sample
	String? description # description for sample
	String? accession_id # accession ID of sample on ENCODE portal

	# mandatory input files
	Array[Array[Array[String]]] fastqs 
						# fastqs[bio_rep_id][tech_rep_id][pair_id]
	Array[Array[Array[String]]] adapters 
						# adapters[bio_rep_id][tech_rep_id][pair_id]
	Array[String] bams 			# if starting from bams
	Array[String] nodup_bams 	# if starting from filtered bams
	Array[String] tas 			# if starting from tag-aligns

	# determin input type and num_rep
	call inputs {
		input : 
			fastqs = fastqs,
			bams = bams,
			nodup_bams = nodup_bams,
			tas = tas,
	}

	# pipeline starts here
	if (inputs.type=='fastq') {
		# trim adapters
		scatter(i in range(inputs.num_rep)) {
			call trim_adapter {
				input:
					fastqs = fastqs[i],
					adapters = if (length(adapters)>0) 
							then adapters[i] else [],
					paired_end = paired_end,
					cpu = cpu/inputs.num_rep,
			}
		}		
		# merge fastqs from technical replicates
		scatter(i in range(inputs.num_rep)) {
			call merge_fastq {
				input: 
					fastqs = trim_adapter.trimmed_fastqs[i]
			}
		}
		# align trimmed/merged fastqs with bowtie2
		scatter(i in range(inputs.num_rep)) {
			call bowtie2 {
				input:
					fastqs = merge_fastq.merged_fastqs[i],
					idx_tar = genome["bowtie2_idx_tar"],
					paired_end = paired_end,
					multimapping = multimapping,
					cpu = cpu/inputs.num_rep,
			}
		}
	}

	if (inputs.type=='fastq' || inputs.type=='bam') {
		# filter/dedup bam
		scatter(i in range(inputs.num_rep)) {
			call filter {
				input:
					bam = select_first([bowtie2.bam,bams])[i],
					paired_end = paired_end,
					multimapping = multimapping,
					cpu = cpu/inputs.num_rep,
			}
		}
	}

	if (inputs.type=='fastq' || inputs.type=='bam' ||
		inputs.type=='nodup_bam') {
		# convert bam to tagalign and subsample it if necessary
		scatter(i in range(inputs.num_rep)) {
			call bam2ta {
				input:
					bam = select_first([filter.nodup_bam,nodup_bams])[i],
					paired_end = select_first([paired_end,true]),
					cpu = cpu/inputs.num_rep,
			}
		}
	}

	# subsample tagalign (non-mito) and cross-correlation analysis
	scatter(i in range(inputs.num_rep)) {
		call xcor {
			input:
				ta = select_first([bam2ta.ta,tas])[i],
				paired_end = select_first([paired_end,true]),
		}
	}

	# make two self pseudo replicates per true replicate
	scatter(i in range(inputs.num_rep)) {
		call spr {
			input:
				ta = select_first([bam2ta.ta,tas])[i],
				paired_end = select_first([paired_end,true]),
		}
	}

	if ( inputs.num_rep>1 ) {
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
	}

	# call peaks on tagalign
	scatter(i in range(inputs.num_rep)) {
		call macs2 {
			input:
				ta = select_first([bam2ta.ta,tas])[i],
				chrsz = genome["chrsz"],
				make_signal = true,
		}
	}
	# filter out peaks with blacklist
	scatter(i in range(inputs.num_rep)) {
		call blacklist_filter as blacklist_filt_macs2 {
			input:
				peak = macs2.npeak[i],
				blacklist = genome["blacklist"],				
		}		
	}
	if ( inputs.num_rep>1 ) {
		# call peaks on pooled replicate
		call macs2 as macs2_pooled {
			input:
				ta = pool_ta.ta_pooled,
				chrsz = genome["chrsz"],
				make_signal = true,
		}
		# filter out peaks with blacklist
		call blacklist_filter as blacklist_filt_macs2_pooled {
			input:
				peak = macs2_pooled.npeak,
				blacklist = genome["blacklist"],
		}
	}
	if ( !select_first([true_rep_only,false]) ) {
		# call peaks on 1st pseudo replicated tagalign 
		scatter(ta in spr.ta_pr1) {
			call macs2 as macs2_pr1 {
				input:
					ta = ta,
					chrsz = genome["chrsz"],
			}
		}
		scatter(peak in macs2_pr1.npeak) {
			# filter out peaks with blacklist
			call blacklist_filter as blacklist_filt_macs2_pr1 {
				input:
					peak = peak,
					blacklist = genome["blacklist"],
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
		scatter(peak in macs2_pr2.npeak) {
			# filter out peaks with blacklist
			call blacklist_filter as blacklist_filt_macs2_pr2 {
				input:
					peak = peak,
					blacklist = genome["blacklist"],
			}		
		}
		if ( inputs.num_rep>1 ) {
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
			# filter out peaks with blacklist
			call blacklist_filter as blacklist_filt_macs2_ppr1 {
				input:
					peak = macs2_ppr1.npeak,
					blacklist = genome["blacklist"],
			}
			call blacklist_filter as blacklist_filt_macs2_ppr2 {
				input:
					peak = macs2_ppr2.npeak,
					blacklist = genome["blacklist"],
			}
		}
	}

	# generate every pair of true replicates
	if ( inputs.num_rep>1 ) {
		call pair_gen { 
			input : num_rep = inputs.num_rep
		}
	}

	# Naive overlap on every pair of true replicates
	if ( inputs.num_rep>1 ) {
		scatter( pair in pair_gen.pairs ) {
			call overlap {
				input : 
					prefix = "rep"+(pair[0]+1)+
							"-rep"+(pair[1]+1),
					peak1 = macs2.npeak[(pair[0])],
					peak2 = macs2.npeak[(pair[1])],
					peak_pooled = macs2_pooled.npeak,
			}
		}
		scatter( peak in overlap.overlap_peak ) {
			call blacklist_filter as blacklist_filt_overlap {
				input:
					peak = peak,
					blacklist = genome["blacklist"],
			}
		}
	}
	if ( !select_first([true_rep_only,false]) ) {
		# Naive overlap on pseduo replicates
		scatter( i in range(inputs.num_rep) ) {
			call overlap as overlap_pr {
				input : 
					prefix = "rep"+(i+1)+"pr",
					peak1 = macs2_pr1.npeak[i],
					peak2 = macs2_pr2.npeak[i],
					peak_pooled = macs2.npeak[i],
			}
		}
		scatter( peak in overlap_pr.overlap_peak ) {
			call blacklist_filter as blacklist_filt_overlap_pr {
				input:
					peak = peak,
					blacklist = genome["blacklist"],
			}
		}
		if ( inputs.num_rep>1 ) {
			# Naive overlap on pooled pseudo replicates
			call overlap as overlap_ppr {
				input : 
					prefix = "ppr",
					peak1 = macs2_ppr1.npeak,
					peak2 = macs2_ppr2.npeak,
					peak_pooled = macs2_pooled.npeak,
			}
			call blacklist_filter as blacklist_filt_overlap_ppr {
				input:
					peak = overlap_ppr.overlap_peak,
					blacklist = genome["blacklist"],
			}
		}
		# reproducibility QC for overlapping peaks
		call reproducibility as reproducibility_overlap {
			input:
				peaks = select_first(
					blacklist_filt_overlap.filtered_peak,[]),
				peaks_pr = blacklist_filt_overlap_pr.filtered_peak,
				peak_ppr = blacklist_filt_overlap_ppr.filtered_peak,
		}
	}

	if ( select_first([enable_idr,false]) ) {
		if ( inputs.num_rep>1 ) {
			scatter( pair in pair_gen.pairs ) {
				# IDR on every pair of true replicates
				call idr {
					input : 
						prefix = "rep"+(pair[0]+1)
								+"-rep"+(pair[1]+1),
						peak1 = macs2.npeak[(pair[0])],
						peak2 = macs2.npeak[(pair[1])],
						peak_pooled = macs2_pooled.npeak,
				}
			}
			scatter( peak in idr.idr_peak ) {
				# filter out peaks with blacklist
				call blacklist_filter as blacklist_filt_idr {
					input:
						peak = peak,
						blacklist = genome["blacklist"],
				}
			}
		}
		if ( !select_first([true_rep_only,false]) ) {
			# IDR on pseduo replicates
			scatter( i in range(inputs.num_rep) ) {
				call idr as idr_pr {
					input : 
						prefix = "rep"+(i+1)+"pr",
						peak1 = macs2_pr1.npeak[i],
						peak2 = macs2_pr2.npeak[i],
						peak_pooled = macs2.npeak[i],
				}
			}
			scatter( peak in idr_pr ) {			
				call blacklist_filter as blacklist_filt_idr_pr {
					input:
						peak = peak,
						blacklist = genome["blacklist"],
				}		
			}
			if ( inputs.num_rep>1 ) {
				call idr as idr_ppr {
					input : 
						prefix = "ppr",
						peak1 = macs2_ppr1.npeak,
						peak2 = macs2_ppr2.npeak,
						peak_pooled = macs2_pooled.npeak,
				}
				call blacklist_filter as blacklist_filt_idr_ppr {
					input:
						peak = idr_ppr.idr_peak,
						blacklist = genome["blacklist"],
				}
			}
			# reproducibility QC for IDR peaks
			call reproducibility as reproducibility_idr {
				input:
					peaks = select_first(
						blacklist_filt_idr.filtered_peak,[]),
					peaks_pr = blacklist_filt_idr_pr.filtered_peak,
					peak_ppr = blacklist_filt_idr_ppr.filtered_peak,
			}
		}
	}	
}