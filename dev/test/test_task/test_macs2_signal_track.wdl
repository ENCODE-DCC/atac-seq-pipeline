version 1.0
import '../../../atac.wdl' as atac
import 'compare_md5sum.wdl' as compare_md5sum

workflow test_macs2_signal_track {
    input {
        Float pval_thresh
        Int smooth_win

        # test macs2 for SE set only
        String se_ta

        String ref_se_macs2_pval_bw # p-val signal

        String se_chrsz
        String se_gensz

        Float macs2_mem_factor = 0.0
        Int macs2_time_hr = 24
        Float macs2_disk_factor = 40.0
        String docker
    }
    RuntimeEnvironment runtime_environment = {
        "docker": docker,
        "singularity": "",
        "conda": ""
    }

    call atac.macs2_signal_track as se_macs2_signal_track { input :
        ta = se_ta,
        gensz = se_gensz,
        chrsz = se_chrsz,
        pval_thresh = pval_thresh,
        smooth_win = smooth_win,

        mem_factor = macs2_mem_factor,
        time_hr = macs2_time_hr,
        disk_factor = macs2_disk_factor,
        runtime_environment = runtime_environment,
    }

    call compare_md5sum.compare_md5sum { input :
        labels = [
            'se_macs2_pval_bw',
        ],
        files = [
            se_macs2_signal_track.pval_bw,
        ],
        ref_files = [
            ref_se_macs2_pval_bw,
        ],
    }
}
