# ENCODE DCC ATAC-Seq/DNase-Seq pipeline tester
# Author: Jin Lee (leepc12@gmail.com)
import "../../atac.wdl" as atac

workflow test_ataqc {
	# subsampled paired end dataset aligned with full chromosome bowtie2 index
	# only full chromosome index works. chr19,chrM-only index fails in metaseq...
	File read_len_log
	File flagstat_log
	File bowtie2_log
	File pbc_log
	File dup_log
	File bam
	File nodup_flagstat_log
	File mito_dup_log
	File nodup_bam
	File ta
	File peak
	File idr_peak
	File overlap_peak
	File bigwig
	Boolean paired_end

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

	#File ref_qc_html
	File ref_qc_txt

	Int ataqc_mem_mb = 16000
	Int ataqc_mem_java_mb = 16000
	Int ataqc_time_hr = 24
	String ataqc_disks = "local-disk 100 HDD"

	call atac.ataqc { input : 
		paired_end = paired_end,
		read_len_log = read_len_log,
		flagstat_log = flagstat_log,
		bowtie2_log = bowtie2_log,
		pbc_log = pbc_log,
		dup_log = dup_log,
		bam = bam,
		nodup_flagstat_log = nodup_flagstat_log,
		mito_dup_log = mito_dup_log,
		nodup_bam = nodup_bam,
		ta = ta,
		peak = peak,
		idr_peak = idr_peak,
		overlap_peak = overlap_peak,
		bigwig = bigwig,
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

		mem_mb = ataqc_mem_mb,
		mem_java_mb = ataqc_mem_java_mb,
		time_hr = ataqc_time_hr,
		disks = ataqc_disks,
	}

	call atac.compare_md5sum { input :
		labels = [
			#'ataqc_html',
			'ataqc_txt',
		],
		files = [
			#ataqc.html,
			ataqc.txt,
		],
		ref_files = [
			#ref_qc_html,
			ref_qc_txt,
		],
	}
}
