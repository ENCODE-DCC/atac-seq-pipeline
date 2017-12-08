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

	# mandaory adapters
	Array[Array[Array[String]]] adapters 
								# [rep_id][merge_id][end_id]
	# mandatory genome param
	String genome_tsv 		# reference genome data TSV file including
							# all important genome specific data file paths
							# and parameters
	Boolean paired_end 		# endedness of sample

	# optional resource (only for SGE and SLURM)
							# name of SGE queue or SLURM partition
							# all sub-tasks of pipeline will be sumitted to two queues
	String? queue_hard 		# queue for hard/long multi-threaded tasks
							# (trim_adapter, bowtie2, filter, bam2ta, macs2)
	String? queue_short 	# queue for easy/short tasks (all others)

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
		if (inputs.type=='fastq') {
			# trim adapters
			call trim_adapter {
				input:
					fastqs = fastqs[i],
					adapters = if length(adapters)>0 
							then adapters[i] else [],
					paired_end = paired_end,
					queue = queue_hard,
			}
			# merge fastqs from technical replicates
			call merge_fastq {
				input: 
					fastqs_R1 = trim_adapter.trimmed_fastqs_R1,
					fastqs_R2 = if paired_end then 
						trim_adapter.trimmed_fastqs_R2 else [],
					paired_end = paired_end,
					queue = queue_short,
			}
			# align trimmed/merged fastqs with bowtie2
			call bowtie2 {
				input:
					idx_tar = inputs.bowtie2_idx_tar,
					fastqs = merge_fastq.merged_fastqs, #[R1,R2]
					paired_end = paired_end,
					multimapping = multimapping,
					queue = queue_hard,
			}
		}
	}

}

# genomic tasks
task trim_adapter { # detect/trim adapter
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
	String? queue

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

	command {
		python $(which encode_dcc_trim_adapter.py) \
			${write_tsv(fastqs)} \
			${"--adapters " + write_tsv(adapters)} \
			${if paired_end then "--paired-end" else ""} \
			${if select_first([auto_detect_adapter,false])
				then "--auto-detect-adapter" else ""} \
			${"--min-trim-len " + min_trim_len} \
			${"--err-rate " + err_rate} \
			${"--nth " + cpu}
	}
	output {
		# Google Cloud does not support 2-dim array in output
		# read_tsv() for JES always fails here (Cromwell bug?)
		Array[File] trimmed_fastqs_R1 = read_lines("out_R1.txt")
		Array[File] trimmed_fastqs_R2 = if paired_end then 
			 read_lines("out_R2.txt") else []
	}
	runtime {
		cpu : "${select_first([cpu,2])}"
		disks: "local-disk 50 SSD"
		queue : queue
	}
}

task merge_fastq { # merge fastqs
	# parameters from workflow
	# Google Cloud does not support 2-dim array in output
	#Array[Array[File]] fastqs # [merge_id][end_id]
	Array[File] fastqs_R1 # [merge_id]
	Array[File] fastqs_R2 # [merge_id]
	Boolean paired_end
	# resource
	Int? cpu
	String? queue

	command {
		python $(which encode_dcc_merge_fastq.py) \
			--fastqs-R1 ${sep=' ' fastqs_R1} \
			${if paired_end then "--paired-end" else ""} \
			--fastqs-R2 ${sep=' ' fastqs_R2}
	}
	output {
		# merged_fastqs[end_id]
		Array[File] merged_fastqs = read_lines("out.txt")
	}
	runtime {
		cpu : "${select_first([cpu,2])}"
		disks: "local-disk 50 SSD"
		queue : queue
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
	String? queue

	command {
		python $(which encode_dcc_bowtie2.py) \
			${idx_tar} \
			${sep=' ' fastqs} \
			${if paired_end then "--paired-end" else ""} \
			${"--multimapping " + multimapping} \
			${"--score-min " + score_min} \
			${"--nth " + cpu}
	}
	output {
		File bam = glob("*.bam")[0]
		File bai = glob("*.bai")[0]
		File align_log = glob("*.align.log")[0]
		File flagstat_qc = glob("*.flagstat.qc")[0]
	}
	runtime {
		cpu : "${select_first([cpu,4])}"
		memory : "${select_first([mem_mb,'20000'])} MB"
		time : "${select_first([time_hr,48])}"
		disks: "local-disk 100 SSD"
		queue : queue
	}
}

# workflow system tasks
# to reduce overhead of provisioning extra nodes
# we have only one task to
# 	1) determine input type and number of replicates	
# 	2) generate pair (rep-x_vs_rep-y) of all true replicate
# 	3) read genome_tsv
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
		Array[Array[Int]] pairs = if num_rep>1 then 
									read_tsv(stdout()) else [[]]
		File ref_fa = read_map(genome_tsv)['ref_fa']
		File bowtie2_idx_tar = read_map(genome_tsv)['bowtie2_idx_tar']
		#File bwa_idx_tar = read_map(genome_tsv)['bwa_idx_tar']
		String blacklist = read_map(genome_tsv)['blacklist']
		String chrsz = read_map(genome_tsv)['chrsz']
		String gensz = read_map(genome_tsv)['gensz']
	}
	runtime {
		disks: "local-disk 20 HDD"		
	}
}
