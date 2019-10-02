# ENCODE DCC ATAC-Seq/DNase-Seq pipeline tester
# Author: Jin Lee (leepc12@gmail.com)
import '../../../atac.wdl' as atac
import 'compare_md5sum.wdl' as compare_md5sum

workflow test_annot_enrich {
	File ta
	File blacklist
	File dnase
	File prom
	File enh

	File ref_annot_enrich_qc

	call atac.annot_enrich { input : 
		ta = ta,
		blacklist = blacklist,
		dnase = dnase,
		prom = prom,
		enh = enh,
	}

	call compare_md5sum.compare_md5sum { input :
		labels = [
			'test_annot_enrich_qc',
		],
		files = [
			annot_enrich.annot_enrich_qc,
		],
		ref_files = [
			ref_annot_enrich_qc,
		],
	}
}
