version 1.0
import '../../../atac.wdl' as atac
import 'compare_md5sum.wdl' as compare_md5sum

workflow test_fraglen_stat_pe {
    input {
        File nodup_bam

        File ref_nucleosomal_qc
        String docker
    }
    RuntimeEnvironment runtime_environment = {
        "docker": docker,
        "singularity": "",
        "conda": ""
    }

    call atac.fraglen_stat_pe { input : 
        nodup_bam = nodup_bam,
        picard_java_heap = '4G',
        runtime_environment = runtime_environment,
    }

    call compare_md5sum.compare_md5sum { input :
        labels = [
            'test_nucleosomal_qc',
        ],
        files = [
            fraglen_stat_pe.nucleosomal_qc,
        ],
        ref_files = [
            ref_nucleosomal_qc,
        ],
    }
}
