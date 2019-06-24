# ENCODE DCC ATAC-Seq/DNase-Seq pipeline tester
# Author: Jin Lee (leepc12@gmail.com)
import "../../atac.wdl" as atac
import "compare_md5sum.wdl" as compare_md5sum

workflow test_ataqc {
	File read_len_log
	File flagstat_qc
	File bowtie2_log
	File pbc_qc
	File dup_qc
	File bam
	File nodup_flagstat_qc
	File mito_dup_log
	File nodup_bam
	File ta
	File peak
	File idr_peak
	File overlap_peak
	File pval_bw
	Boolean paired_end

	File ref_fa
	File chrsz
	File tss
	File blacklist
	File dnase
	File prom
	File enh
	File reg2map_bed
	File reg2map
	File roadmap_meta
	String mito_chr_name = 'chrM'

	#File ref_qc_html
	File ref_qc_txt
	File ref_qc_align_only_1_txt
	File ref_qc_align_only_2_txt
	File ref_qc_peak_only_txt
	File ref_qc_tss_enrich_txt

	Int ataqc_mem_mb = 7000
	Int ataqc_mem_java_mb = 6000
	Int ataqc_time_hr = 24
	String ataqc_disks = "local-disk 100 HDD"

	call atac.ataqc { input : 
		paired_end = paired_end,
		read_len_log = read_len_log,
		flagstat_qc = flagstat_qc,
		bowtie2_log = bowtie2_log,
		pbc_qc = pbc_qc,
		dup_qc = dup_qc,
		bam = bam,
		nodup_flagstat_qc = nodup_flagstat_qc,
		mito_dup_log = mito_dup_log,
		nodup_bam = nodup_bam,
		ta = ta,
		peak = peak,
		idr_peak = idr_peak,
		overlap_peak = overlap_peak,
		pval_bw = pval_bw,

		ref_fa = ref_fa,
		chrsz = chrsz,
		tss = tss,
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
		mito_chr_name = mito_chr_name,
	}

	call compare_md5sum.compare_md5sum { input :
		labels = [
			'ataqc_txt',
			'ataqc_align_tss_enrich_txt',
		],
		files = [
			ataqc.txt,
			ataqc_tss_enrich.txt,
		],
		ref_files = [
			ref_qc_txt,
			ref_qc_tss_enrich_txt,
		],
	}
}
