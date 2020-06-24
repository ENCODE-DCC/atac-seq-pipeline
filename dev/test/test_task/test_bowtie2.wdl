version 1.0
import '../../../atac.wdl' as atac
import 'compare_md5sum.wdl' as compare_md5sum

workflow test_bowtie2 {
    input {
        Int multimapping

        Array[String] pe_fastqs_R1
        Array[String] pe_fastqs_R2
        Array[String] se_fastqs_R1

        String se_chrsz
        String pe_chrsz
        String cutadapt_param = '-e 0.1 -m 5'

        # we don't compare BAM because BAM's header includes date
        # hence md5sums don't match all the time
        String ref_pe_flagstat
        String ref_pe_flagstat_no_multimapping

        String ref_se_flagstat
        String ref_se_flagstat_no_multimapping

        String pe_bowtie2_idx_tar
        String se_bowtie2_idx_tar

        Int bowtie2_cpu = 1
        Int bowtie2_mem_mb = 20000
        Int bowtie2_time_hr = 48
        String bowtie2_disks = 'local-disk 100 HDD'
    }

    call atac.align as pe_bowtie2 { input :
        aligner = 'bowtie2',
        idx_tar = pe_bowtie2_idx_tar,
        mito_chr_name = 'chrM',
        fastqs_R1 = pe_fastqs_R1,
        fastqs_R2 = pe_fastqs_R2,
        adapters_R1 = [],
        adapters_R2 = [],
        cutadapt_param = cutadapt_param,
        multimapping = multimapping,
        paired_end = true,
        chrsz = pe_chrsz,
        auto_detect_adapter = true,

        cpu = bowtie2_cpu,
        mem_mb = bowtie2_mem_mb,
        time_hr = bowtie2_time_hr,
        disks = bowtie2_disks,
    }
    call atac.align as pe_bowtie2_no_multimapping { input :
        aligner = 'bowtie2',
        idx_tar = pe_bowtie2_idx_tar,
        mito_chr_name = 'chrM',
        fastqs_R1 = pe_fastqs_R1,
        fastqs_R2 = pe_fastqs_R2,
        adapters_R1 = [],
        adapters_R2 = [],
        cutadapt_param = cutadapt_param,
        multimapping = 0,
        paired_end = true,
        chrsz = pe_chrsz,
        auto_detect_adapter = true,

        cpu = bowtie2_cpu,
        mem_mb = bowtie2_mem_mb,
        time_hr = bowtie2_time_hr,
        disks = bowtie2_disks,
    }
    call atac.align as se_bowtie2 { input :
        aligner = 'bowtie2',
        idx_tar = se_bowtie2_idx_tar,
        mito_chr_name = 'chrM',
        fastqs_R1 = se_fastqs_R1,
        fastqs_R2 = [],
        adapters_R1 = [],
        adapters_R2 = [],
        cutadapt_param = cutadapt_param,
        multimapping = multimapping,
        paired_end = false,
        chrsz = se_chrsz,
        auto_detect_adapter = true,

        cpu = bowtie2_cpu,
        mem_mb = bowtie2_mem_mb,
        time_hr = bowtie2_time_hr,
        disks = bowtie2_disks,
    }
    call atac.align as se_bowtie2_no_multimapping { input :
        aligner = 'bowtie2',
        idx_tar = se_bowtie2_idx_tar,
        mito_chr_name = 'chrM',
        fastqs_R1 = se_fastqs_R1,
        fastqs_R2 = [],
        adapters_R1 = [],
        adapters_R2 = [],
        cutadapt_param = cutadapt_param,
        multimapping = 0,
        paired_end = false,
        chrsz = se_chrsz,
        auto_detect_adapter = true,

        cpu = bowtie2_cpu,
        mem_mb = bowtie2_mem_mb,
        time_hr = bowtie2_time_hr,
        disks = bowtie2_disks,
    }

    call compare_md5sum.compare_md5sum { input :
        labels = [
            'pe_bowtie2',
            'pe_bowtie2_no_multimapping',
            'se_bowtie2',
            'se_bowtie2_no_multimapping',
        ],
        files = [
            pe_bowtie2.samstat_qc,
            pe_bowtie2_no_multimapping.samstat_qc,
            se_bowtie2.samstat_qc,
            se_bowtie2_no_multimapping.samstat_qc,
        ],
        ref_files = [
            ref_pe_flagstat,
            ref_pe_flagstat_no_multimapping,
            ref_se_flagstat,
            ref_se_flagstat_no_multimapping,
        ],
    }
}
