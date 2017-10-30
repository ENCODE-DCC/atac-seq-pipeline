task trim_merge { # detect/trim adapter
	Array[Array[File]] fastqs # fastqs[tech_rep_id][pair_id]
	Array[Array[String]]? adapters # adapters[tech_rep_id][pair_id]
	Boolean? paired_end 
	Boolean? auto_detect_adapter
	Int? min_trim_len
	Float? err_rate
	# resource
	Int? cpu
	Int? mem_mb
	Int? time_hr
	String? queue

	command {
		python $(which encode_dcc_trim_merge.py) \
			${write_tsv(fastqs)} \
			${if select_first([paired_end,false])
				then "--paired-end" else ""} \
			${"--adapters " + write_tsv(adapters)} \
			${if select_first([auto_detect_adapter,false])
				then "--auto-detect-adapter" else ""} \
			${"--min-trim-len " + min_trim_len} \
			${"--err-rate " + err_rate} \
			${"--nth " + cpu} \
			--out-dir . \
			--out-meta-tsv meta.tmp
	}
	output {
		Array[File] trimmed_fastqs = read_tsv("meta.tmp")[0]
	}
	runtime {
        cpu : "${select_first([cpu,1])}"
        memory : "${select_first([mem_mb,'100'])} MB"
        time : "${select_first([time_hr,1])}"
        queue : queue
	}
}

task bowtie2 {
	Array[File] fastqs # fastqs[pair_id]
	Boolean? paired_end
	Int? multimapping
	File idx_tar # reference index tar
	# resource
	Int? cpu
	Int? mem_mb
	Int? time_hr
	String? queue

	command {
		python $(which encode_dcc_bowtie2.py) \
			${idx_tar} \
			${sep=' ' fastqs} \
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
	# input from workflow
	File bam
	Boolean? paired_end
	Int? multimapping
	String? dup_marker
	Int? mapq_thresh
	Boolean? no_dup_removal
	# resource
	Int? cpu
	Int? mem_mb
	Int? time_hr
	String? queue

	command {
		python $(which encode_dcc_filter.py) \
			${bam} \
			${if select_first([paired_end,false])
				then "--paired-end" else ""} \
			${"--multimapping " + multimapping} \
			${"--dup-marker " + dup_marker} \
			${"--mapq-thresh " + mapq_thresh} \
			${if select_first([no_dup_removal,false])
				then "--no-dup-removal" else ""} \
			${"--nth " + cpu}			
	}
	output {
		File dedup_bam = glob("*.nodup.bam")[0]
		File dedup_bai = glob("*.bai")[0]
		File flagstat_qc = glob("*.flagstat.qc")[0]
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
	Boolean? paired_end
	Boolean? disable_tn5_shift # for dnase-seq
	Int? subsample # number of reads to subsample
	# resource
	Int? cpu
	Int? mem_mb
	Int? time_hr
	String? queue

	command {
		python $(which encode_dcc_bam2ta.py) \
			${bam} \
			${if select_first([paired_end,false])
				then "--paired-end" else ""} \
			${if select_first([disable_tn5_shift,false])
				then "--disable-tn5-shift" else ""} \
			${"--nth " + cpu}			
			${"--subsample " + subsample}
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
	Boolean? paired_end
	# resource
	Int? cpu
	Int? mem_mb
	Int? time_hr
	String? queue

	command {
		python $(which encode_dcc_spr.py) \
			${if select_first([paired_end,false])
				then "--paired-end" else ""} \
			${ta}
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
		File ta_pooled = glob("*.pooled.tagAlign.gz")[0]
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
	Int? cap_num_peak # cap number of peaks with top scores
	String? rank # ranking column (pval, qval, ...)
	# resource
	Int? cpu
	Int? mem_mb
	Int? time_hr
	String? queue

	command {
		python $(which encode_dcc_macs2.py) \
			${ta} \
			${"--cap-num-peak " + cap_num_peak}
			${"--rank " + rank}
	}
	output {
		File npeak = glob("*.narrowPeak.gz")[0]
	}
	runtime {
        cpu : "${select_first([cpu,1])}"
        memory : "${select_first([mem_mb,'100'])} MB"
        time : "${select_first([time_hr,1])}"
        queue : queue
	}
}

#task naive_overlap { # including FRiP calculation
#}

#task idr { # including FRiP calculation
#}

workflow atac {
	# general
	File genome_tsv # reference genome data TSV file
	Map[String,String] genome = read_map(genome_tsv)
	Array[Array[Array[File]]] fastqs # fastqs[bio_rep_id][tech_rep_id][pair_id]
	Array[Array[Array[String]]]? adapters # adapters[bio_rep_id][tech_rep_id][pair_id]
	Boolean? paired_end
	Int? multimapping # used for bowtie2 and filter
	Int? cpu # must be multiples of number of replicates

	# trim adapters and merge fastqs from technical replicates
	scatter(i in range(length(fastqs))) {
		call trim_merge {
			input:
				fastqs = fastqs[i],
				adapters = adapters[i],
				paired_end = paired_end,
				cpu = cpu/length(fastqs),
		}
	}

	# align trimmed/merged fastqs with bowtie2
	scatter(fastqs in trim_merge.trimmed_fastqs) {
		call bowtie2 {
			input:
				fastqs = fastqs,
				idx_tar = sub(genome["bowtie2_idx_tar"],"[\r\n]+",""),
				paired_end = paired_end,
				multimapping = multimapping,
				cpu = cpu/length(fastqs),
		}
	}
	
	# filter/dedup bam
	scatter(bam in bowtie2.bam) {
		call filter {
			input:
				bam = bam,
				paired_end = paired_end,
				multimapping = multimapping,
				cpu = cpu/length(trim_merge.trimmed_fastqs),
		}
	}

	# convert bam to tagalign
	scatter(bam in filter.dedup_bam) {
		call bam2ta {
			input:
				bam = bam,
				paired_end = select_first([paired_end,true]),
				cpu = cpu/length(filter.dedup_bam),
		}
	}

	# make self-pseudo replicates
	scatter(ta in bam2ta.ta) {
		call spr {
			input:
				ta = ta,
				paired_end = select_first([paired_end,true]),
		}
	}

	call pool_ta {
		input :
			tas = bam2ta.ta,
	}

	call pool_ta as pool_ta_pr1 {
		input :
			tas = spr.ta_pr1,
	}

	call pool_ta as pool_ta_pr2 {
		input :
			tas = spr.ta_pr2,
	}

#	# cross-correlation analysis on subsampled tagalign
#	scatter(ta in bam2ta.ta) {
#		call xcor {
#			input:
#				ta = ta,
#				paired_end = select_first([paired_end,true]),
#		}
#	}

#	# call peaks on tagalign with macs2
#	scatter(ta in bam2ta.ta) {
#		call macs2 as macs2_tr {
#			input:
#				ta = ta,
#				chrsz = sub(genome["chrsz"],"[\r\n]+",""),
#		}
#	}

#	# call peaks on pseudo replicated tagalign with macs2
#	scatter(ta in bam2ta.ta) {
#		call macs2 as macs2_pr {
#			input:
#				ta = ta,
#				chrsz = sub(genome["chrsz"],"[\r\n]+",""),
#		}
#	}

#	# call peaks on pooled pseudo replicates
#	call macs2 as macs2_pooled {
#		input:
#			ta = pool_ta_tr.pooled_ta,
#			chrsz = sub(genome["chrsz"],"[\r\n]+",""),
#	}

#	# call peaks on pseudo replicated pooled tagalign with macs2
#	call macs2 as macs2_ppr {
#		input:
#			ta = pool_ta_pr.pooled_ta,
#			chrsz = sub(genome["chrsz"],"[\r\n]+",""),
#	}

#	# naive overlap
#	call naive_overlap {
#		tas = bam2ta.ta,
#		npeaks_tr = macs2_tr.npeak
#		npeaks_pr = macs2_pr.npeak
#		ppr_ta = pool_ta_pr.pooled_ta,
#		blacklist = sub(genome["blacklist"],"[\r\n]+",""),
#	}

#	# naive overlap and idr
#	call idr {
#		tas = bam2ta.ta,
#		pr_tas = spr.pr_ta,
#		pooled_ta = pool_ta_tr.pooled_ta,
#		ppr_ta = pool_ta_pr.pooled_ta,
#		blacklist = sub(genome["blacklist"],"[\r\n]+",""),
#	}

}
