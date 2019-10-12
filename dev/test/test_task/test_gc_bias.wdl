# ENCODE DCC ATAC-Seq/DNase-Seq pipeline tester
# Author: Jin Lee (leepc12@gmail.com)
import '../../../atac.wdl' as atac
import 'compare_md5sum.wdl' as compare_md5sum

workflow test_gc_bias {
	File nodup_bam

	File ref_fa

	File ref_gc_log

	call atac.gc_bias { input : 
		nodup_bam = nodup_bam,

		ref_fa = ref_fa,
	}

	call remove_comments_from_gc_log { input :
		gc_log = gc_bias.gc_log
	}

	call remove_comments_from_gc_log as remove_comments_from_gc_log_ref { input :
		gc_log = ref_gc_log
	}

	call compare_md5sum.compare_md5sum { input :
		labels = [
			'test_gc_log',
		],
		files = [
			remove_comments_from_gc_log.filt_gc_log,
		],
		ref_files = [
			remove_comments_from_gc_log_ref.filt_gc_log,
		],
	}
}

task remove_comments_from_gc_log {
	File gc_log
	command {
		zcat -f ${gc_log} | grep -v '# ' \
			> ${basename(gc_log) + '.date_filt_out'}
	}
	output {
		File filt_gc_log = glob('*.date_filt_out')[0]
	}
}
