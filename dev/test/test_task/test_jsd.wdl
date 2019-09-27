# ENCODE DCC atac-Seq pipeline tester for task jsd
# Author: Jin Lee (leepc12@gmail.com)
import '../../../atac.wdl' as atac
import 'compare_md5sum.wdl' as compare_md5sum

workflow test_jsd {
	Array[File] se_nodup_bams
	File se_blacklist
	File se_fake_blacklist
	Array[File] ref_se_jsd_logs
	Array[File] ref_se_jsd_fake_blacklist_logs
	# task level test data (BAM) is generated from BWA
	# so we keep using 30 here, this should be 255 for bowtie2 BAMs
	Int mapq_thresh = 30

	Int jsd_cpu = 1
	Int jsd_mem_mb = 12000
	Int jsd_time_hr = 6
	String jsd_disks = 'local-disk 100 HDD'

	call atac.jsd as se_jsd { input :
		nodup_bams = se_nodup_bams,
		blacklist = se_blacklist,
		mapq_thresh = mapq_thresh,

		cpu = jsd_cpu,
		mem_mb = jsd_mem_mb,
		time_hr = jsd_time_hr,
		disks = jsd_disks,
	}

	call atac.jsd as se_jsd_fake_blacklist { input :
		nodup_bams = se_nodup_bams,
		blacklist = se_fake_blacklist,
		mapq_thresh = mapq_thresh,

		cpu = jsd_cpu,
		mem_mb = jsd_mem_mb,
		time_hr = jsd_time_hr,
		disks = jsd_disks,
	}

	# take first 8 columns (vaule in other columns are random)
	#scatter(i in range(2)){
	#	call take_8_cols { input :
	#		f = se_jsd.jsd_qcs[i],
	#	}
	#	call take_8_cols as ref_take_8_cols { input :
	#		f = ref_se_jsd_logs[i],
	#	}
	#}

	call compare_md5sum.compare_md5sum { input :
		labels = [
			'se_jsd_rep1',
			'se_jsd_fake_blacklist_rep1',
		],
		files = [
			#take_8_cols.out[0],
			se_jsd.jsd_qcs[0],
			se_jsd_fake_blacklist.jsd_qcs[0],
		],
		ref_files = [
			#ref_take_8_cols.out[0],
			ref_se_jsd_logs[0],
			ref_se_jsd_fake_blacklist_logs[0],
		],
	}
}

task take_8_cols {
	File f
	command {
		cut -f 1-8 ${f} > out.txt
	}
	output {
		File out = 'out.txt'
	}
}
