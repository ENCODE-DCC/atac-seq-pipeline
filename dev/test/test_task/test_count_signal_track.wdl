version 1.0
import '../../../atac.wdl' as atac
import 'compare_md5sum.wdl' as compare_md5sum

workflow test_count_signal_track {
    input {    
        String se_ta

        String ref_se_count_signal_track_pos_bw
        String ref_se_count_signal_track_neg_bw

        String se_chrsz
        String docker
    }
    RuntimeEnvironment runtime_environment = {
        "docker": docker,
        "singularity": "",
        "conda": ""
    }

    call atac.count_signal_track as se_count_signal_track { input :
        ta = se_ta,
        chrsz = se_chrsz,
        runtime_environment = runtime_environment,
    }

    call compare_md5sum.compare_md5sum { input :
        labels = [
            'se_count_signal_track_pos_bw',
            'se_count_signal_track_neg_bw',
        ],
        files = [
            se_count_signal_track.pos_bw,
            se_count_signal_track.neg_bw,
        ],
        ref_files = [
            ref_se_count_signal_track_pos_bw,
            ref_se_count_signal_track_neg_bw,
        ],
    }
}
