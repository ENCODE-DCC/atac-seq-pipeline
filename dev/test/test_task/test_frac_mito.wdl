version 1.0
import '../../../atac.wdl' as atac
import 'compare_md5sum.wdl' as compare_md5sum

workflow test_frac_mito {
    input {
        File non_mito_samstat
        File mito_samstat

        File ref_frac_mito_qc
    }

    call atac.frac_mito as frac_mito { input:
        non_mito_samstat = non_mito_samstat,
        mito_samstat = mito_samstat,
    }

    call compare_md5sum.compare_md5sum { input :
        labels = [
            'frac_mito', 
        ],
        files = [
            frac_mito.frac_mito_qc,
        ],
        ref_files = [
            ref_frac_mito_qc,
        ],
    }
}
