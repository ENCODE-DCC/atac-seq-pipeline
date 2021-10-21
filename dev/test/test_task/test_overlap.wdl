version 1.0
import '../../../atac.wdl' as atac
import 'compare_md5sum.wdl' as compare_md5sum

workflow test_overlap {
    input {
        String se_peak_rep1 # test overlap,idr for SE set only
        String se_peak_rep2
        String se_peak_pooled
        String se_ta_pooled

        String ref_se_overlap_peak
        String ref_se_overlap_bfilt_peak
        String ref_se_overlap_frip_qc

        String se_blacklist
        String se_chrsz

        String regex_bfilt_peak_chr_name = 'chr[\\dXY]+'
        String docker
    }
    RuntimeEnvironment runtime_environment = {
        "docker": docker,
        "singularity": "",
        "conda": ""
    }

    call atac.overlap as se_overlap { input :
        prefix = 'rep1-rep2',
        peak1 = se_peak_rep1,
        peak2 = se_peak_rep2,
        peak_pooled = se_peak_pooled,
        peak_type = 'narrowPeak',
        blacklist = se_blacklist,
        regex_bfilt_peak_chr_name = regex_bfilt_peak_chr_name,
        chrsz = se_chrsz,
        ta = se_ta_pooled,
        runtime_environment = runtime_environment,
    }

    call compare_md5sum.compare_md5sum { input :
        labels = [
            'se_overlap_peak',
            'se_overlap_bfilt_peak',
            'se_overlap_frip_qc',
        ],
        files = [
            se_overlap.overlap_peak,
            se_overlap.bfilt_overlap_peak,
            se_overlap.frip_qc,
        ],
        ref_files = [
            ref_se_overlap_peak,
            ref_se_overlap_bfilt_peak,
            ref_se_overlap_frip_qc,
        ],
    }
}
