version 1.0
import '../../../atac.wdl' as atac
import 'compare_md5sum.wdl' as compare_md5sum

workflow test_tss_enrich {
    input {
        File read_len_log
        File nodup_bam
        File tss
        File chrsz

        File ref_tss_enrich_qc
        String docker
    }
    RuntimeEnvironment runtime_environment = {
        "docker": docker,
        "singularity": "",
        "conda": ""
    }

    Int? read_len_ = read_int(read_len_log)

    call atac.tss_enrich { input : 
        read_len = read_len_,
        nodup_bam = nodup_bam,
        chrsz = chrsz,
        tss = tss,
        runtime_environment = runtime_environment,
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
