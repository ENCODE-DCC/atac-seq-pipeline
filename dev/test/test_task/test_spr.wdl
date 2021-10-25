version 1.0
import '../../../atac.wdl' as atac
import 'compare_md5sum.wdl' as compare_md5sum

workflow test_spr {
    input {
        File pe_ta
        File se_ta

        File ref_pe_ta_pr1
        File ref_pe_ta_pr2
        File ref_se_ta_pr1
        File ref_se_ta_pr2
        File ref_pe_seed_10_ta_pr1
        File ref_pe_seed_10_ta_pr2
        File ref_se_seed_10_ta_pr1
        File ref_se_seed_10_ta_pr2

        Float spr_mem_factor = 0.0
        Float spr_disk_factor = 6.0
        String docker
    }
    RuntimeEnvironment runtime_environment = {
        "docker": docker,
        "singularity": "",
        "conda": ""
    }

    call atac.spr as pe_spr { input :
        ta = pe_ta,
        paired_end = true,
        pseudoreplication_random_seed = 0,
        mem_factor = spr_mem_factor,
        disk_factor = spr_disk_factor,
        runtime_environment = runtime_environment,
    }    
    call atac.spr as se_spr { input :
        ta = se_ta,
        paired_end = false,
        pseudoreplication_random_seed = 0,
        mem_factor = spr_mem_factor,
        disk_factor = spr_disk_factor,
        runtime_environment = runtime_environment,
    }
    call atac.spr as pe_spr_seed_10 { input :
        ta = pe_ta,
        paired_end = true,
        pseudoreplication_random_seed = 10,
        mem_factor = spr_mem_factor,
        disk_factor = spr_disk_factor,
        runtime_environment = runtime_environment,
    }
    call atac.spr as se_spr_seed_10 { input :
        ta = se_ta,
        paired_end = false,
        pseudoreplication_random_seed = 10,
        mem_factor = spr_mem_factor,
        disk_factor = spr_disk_factor,
        runtime_environment = runtime_environment,
    }

    call compare_md5sum.compare_md5sum { input :
        labels = [
            'pe_spr_pr1',
            'pe_spr_pr2',
            'se_spr_pr1',
            'se_spr_pr2',
            'pe_spr_seed_10_pr1',
            'pe_spr_seed_10_pr2',
            'se_spr_seed_10_pr1',
            'se_spr_seed_10_pr2',
        ],
        files = [
            pe_spr.ta_pr1,
            pe_spr.ta_pr2,
            se_spr.ta_pr1,
            se_spr.ta_pr2,
            pe_spr_seed_10.ta_pr1,
            pe_spr_seed_10.ta_pr2,
            se_spr_seed_10.ta_pr1,
            se_spr_seed_10.ta_pr2,
        ],
        ref_files = [
            ref_pe_ta_pr1,
            ref_pe_ta_pr2,
            ref_se_ta_pr1,
            ref_se_ta_pr2,
            ref_pe_seed_10_ta_pr1,
            ref_pe_seed_10_ta_pr2,
            ref_se_seed_10_ta_pr1,
            ref_se_seed_10_ta_pr2,
        ],
    }
}
