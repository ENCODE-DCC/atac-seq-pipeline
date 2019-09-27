# ENCODE DCC ATAC-Seq/DNase-Seq pipeline tester
# Author: Jin Lee (leepc12@gmail.com)
import '../../../atac.wdl' as atac
import 'compare_md5sum.wdl' as compare_md5sum

workflow test_tss_enrich {
	File read_len_log
	File nodup_bam
	File tss
	File chrsz

	File ref_tss_enrich_qc

	call atac.tss_enrich { input : 
		read_len_log = read_len_log,
		nodup_bam = nodup_bam,
		chrsz = chrsz,
		tss = tss,
	}

	call compare_md5sum.compare_md5sum { input :
		labels = [
			'test_tss_enrich_qc',
		],
		files = [
			tss_enrich.tss_enrich_qc,
		],
		ref_files = [
			ref_tss_enrich_qc,
		],
	}
}
