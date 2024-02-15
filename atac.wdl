version 1.0

struct RuntimeEnvironment {
    String docker
    String singularity
    String conda
}

workflow atac {
    String pipeline_ver = 'v2.2.3'

    meta {
        version: 'v2.2.3'

        author: 'Jin wook Lee'
        email: 'leepc12@gmail.com'
        description: 'ENCODE ATAC-Seq pipeline. See https://github.com/ENCODE-DCC/atac-seq-pipeline for more details. e.g. example input JSON for Terra/Anvil.'
        organization: 'ENCODE DCC'

        specification_document: 'https://docs.google.com/document/d/1f0Cm4vRyDQDu0bMehHD7P7KOMxTOP-HiNoIvL1VcBt8/edit?usp=sharing'

        default_docker: 'encodedcc/atac-seq-pipeline:v2.2.3'
        default_singularity: 'https://encode-pipeline-singularity-image.s3.us-west-2.amazonaws.com/atac-seq-pipeline_v2.2.3.sif'
        default_conda: 'encd-atac'
        croo_out_def: 'https://storage.googleapis.com/encode-pipeline-output-definition/atac.croo.v5.json'

        parameter_group: {
            runtime_environment: {
                title: 'Runtime environment',
                description: 'Runtime environment such as container URIs (Docker, Singularity) and Conda environment name.'
            },
            pipeline_metadata: {
                title: 'Pipeline metadata',
                description: 'Metadata for a pipeline (e.g. title and description).'
            },
            reference_genome: {
                title: 'Reference genome',
                description: 'Genome specific files. e.g. reference FASTA, bowtie2 index, chromosome sizes file.',
                help: 'Choose one atac.genome_tsv file that defines all genome specific parameters in it or define each genome specific parameter in input JSON to override those defined in genome TSV file. If you use Caper then use https://storage.googleapis.com/encode-pipeline-genome-data/genome_tsv/v1/[GENOME]_caper.tsv. Caper will automatically download/install all files defined in such TSV. Otherwise download genome TSV file by using a shell script (scripts/download_genome_data.sh [GENOME] [DEST_DIR]). Supported genomes are hg38, hg19, mm10 and mm9. See pipeline documentation if you want to build genome database from your own FASTA file. If some genome data are missing then analyses using such data will be skipped.'
            },
            input_genomic_data: {
                title: 'Input genomic data',
                description: 'Genomic input files for experiment.',
                help: 'Pipeline can start with any types of experiment data (e.g. FASTQ, BAM, NODUP_BAM, TAG-ALIGN, PEAK). Choose one type and leave others empty. FASTQs have a variable for each biological replicate. e.g. atac.fastqs_rep1_R1 and atac.fastqs_rep2_R1. You can define up to 10 experiment replicates. For other types, there is an array to define file for each biological replicate. e.g. atac.bams: ["rep1.bam", "rep1.bam"]. Define sequential endedness with atac.paired_end, if you have mixed SE and PE replicates then define atac.paired_ends instead for each replicate. e.g. atac.paired_ends: [false, true].'
            },
            adapter_trimming: {
                title: 'Adapter trimming',
                description: 'Parameters for adapter trimming.',
                help: 'Use atac.auto_detect_adapter to automatically detect/trim 3 adapters (Illumina: AGATCGGAAGAGC, Nextera: CTGTCTCTTATA, smallRNA: TGGAATTCTCGG) or manually define adapter sequence to be trimmed (atac.adapter or atac.adapters_repX_RY). Leave all parameters undefined/empty if your FASTQs are already trimmed.'
            },
            pipeline_parameter: {
                title: 'Pipeline parameter',
                description: 'Pipeline type and flags to turn on/off analyses.',
                help: 'Pipeline can run as DNAse-seq mode. The only difference is TN5-shifting of read in ATAC-seq mode. Use atac.align_only to align FASTQs without peak calling.'
            },
            alignment: {
                title: 'Alignment',
                description: 'Parameters for alignment.',
                help: 'Pipeline calculates mitochondrial fraction of reads in raw BAM. But after that it filters out mitochondrial reads (e.g. chrM, MT) from NODUP_BAMs (filtered/deduped). It is controlled by atac.filter_chrs array. If you want to keep mitochondrial reads then make this array empty.'
            },
            peak_calling: {
                title: 'Peak calling',
                description: 'Parameters for peak calling.',
                help: 'This group includes statistical thresholds for peak-calling or post-peak-calling analyses: p-val, FDR, IDR.'
            },
            resource_parameter: {
                title: 'Resource parameter',
                description: 'Number of CPUs (threads), max. memory and walltime for tasks.',
                help: 'Resource settings are used for determining an instance type on cloud backends (e.g. GCP, AWS) and used for submitting tasks to a cluster engine (e.g. SLURM, SGE, ...). Walltime (atac.*_time_hr) is only used for cluster engines. Other tasks default to use 1 CPU and 4GB of memory.'
            }
        }
    }
    input {
        # group: runtime_environment
        String docker = 'encodedcc/atac-seq-pipeline:v2.2.3'
        String singularity = 'https://encode-pipeline-singularity-image.s3.us-west-2.amazonaws.com/atac-seq-pipeline_v2.2.3.sif'
        String conda = 'encd-atac'
        String conda_macs2 = 'encd-atac-macs2'
        String conda_spp = 'encd-atac-spp'
        String conda_python2 = 'encd-atac-py2'

        # group: pipeline_metadata
        String title = 'Untitled'
        String description = 'No description'

        # group: reference_genome
        File? genome_tsv
        String? genome_name
        File? ref_fa
        File? ref_mito_fa
        File? bowtie2_idx_tar
        File? bowtie2_mito_idx_tar
        File? chrsz
        File? blacklist
        File? blacklist2
        String? mito_chr_name
        String? regex_bfilt_peak_chr_name
        String? gensz
        File? tss
        File? dnase
        File? prom
        File? enh
        File? reg2map
        File? reg2map_bed
        File? roadmap_meta

        # group: input_genomic_data
        Boolean? paired_end
        Array[Boolean] paired_ends = []
        Array[File] fastqs_rep1_R1 = []
        Array[File] fastqs_rep1_R2 = []
        Array[File] fastqs_rep2_R1 = []
        Array[File] fastqs_rep2_R2 = []
        Array[File] fastqs_rep3_R1 = []
        Array[File] fastqs_rep3_R2 = []
        Array[File] fastqs_rep4_R1 = []
        Array[File] fastqs_rep4_R2 = []
        Array[File] fastqs_rep5_R1 = []
        Array[File] fastqs_rep5_R2 = []
        Array[File] fastqs_rep6_R1 = []
        Array[File] fastqs_rep6_R2 = []
        Array[File] fastqs_rep7_R1 = []
        Array[File] fastqs_rep7_R2 = []
        Array[File] fastqs_rep8_R1 = []
        Array[File] fastqs_rep8_R2 = []
        Array[File] fastqs_rep9_R1 = []
        Array[File] fastqs_rep9_R2 = []
        Array[File] fastqs_rep10_R1 = []
        Array[File] fastqs_rep10_R2 = []
        Array[File] bams = []
        Array[File] nodup_bams = []
        Array[File] tas = []
        Array[File] peaks = []
        Array[File] peaks_pr1 = []
        Array[File] peaks_pr2 = []
        File? peak_pooled
        File? peak_ppr1
        File? peak_ppr2

        # group: pipeline_parameter
        String pipeline_type = 'atac'
        Boolean align_only = false
        Boolean true_rep_only = false
        Boolean enable_xcor = false
        Boolean enable_count_signal_track = false
        Boolean enable_idr = true
        Boolean enable_preseq = false
        Boolean enable_fraglen_stat = true
        Boolean enable_tss_enrich = true
        Boolean enable_annot_enrich = true
        Boolean enable_jsd = true
        Boolean enable_compare_to_roadmap = false
        Boolean enable_gc_bias = true

        # group: adapter_trimming
        String cutadapt_param = '-e 0.1 -m 5'
        Boolean auto_detect_adapter = false
        String? adapter
        Array[String] adapters_rep1_R1 = []
        Array[String] adapters_rep1_R2 = []
        Array[String] adapters_rep2_R1 = []
        Array[String] adapters_rep2_R2 = []
        Array[String] adapters_rep3_R1 = []
        Array[String] adapters_rep3_R2 = []
        Array[String] adapters_rep4_R1 = []
        Array[String] adapters_rep4_R2 = []
        Array[String] adapters_rep5_R1 = []
        Array[String] adapters_rep5_R2 = []
        Array[String] adapters_rep6_R1 = []
        Array[String] adapters_rep6_R2 = []
        Array[String] adapters_rep7_R1 = []
        Array[String] adapters_rep7_R2 = []
        Array[String] adapters_rep8_R1 = []
        Array[String] adapters_rep8_R2 = []
        Array[String] adapters_rep9_R1 = []
        Array[String] adapters_rep9_R2 = []
        Array[String] adapters_rep10_R1 = []
        Array[String] adapters_rep10_R2 = []

        # group: alignment
        Int multimapping = 4
        String dup_marker = 'picard'
        Boolean no_dup_removal = false
        Int mapq_thresh = 30
        Array[String] filter_chrs = ['chrM', 'MT']
        Int subsample_reads = 0
        Int xcor_subsample_reads = 25000000
        Array[Int?] read_len = []
        Int pseudoreplication_random_seed = 0

        # group: peak_calling
        Int cap_num_peak = 300000
        Float pval_thresh = 0.01
        Int smooth_win = 150
        Float idr_thresh = 0.05

        # group: resource_parameter
        Int align_cpu = 6
        Float align_mem_factor = 0.15
        Int align_time_hr = 48
        Float align_disk_factor = 8.0

        Int filter_cpu = 4
        Float filter_mem_factor = 0.4
        Int filter_time_hr = 24
        Float filter_disk_factor = 8.0

        Int bam2ta_cpu = 2
        Float bam2ta_mem_factor = 0.3
        Int bam2ta_time_hr = 12
        Float bam2ta_disk_factor = 4.0

        Float spr_mem_factor = 20.0
        Float spr_disk_factor = 30.0

        Int jsd_cpu = 4
        Float jsd_mem_factor = 0.1
        Int jsd_time_hr = 12
        Float jsd_disk_factor = 2.0

        Int xcor_cpu = 2
        Float xcor_mem_factor = 1.0
        Int xcor_time_hr = 6
        Float xcor_disk_factor = 4.5

        Int call_peak_cpu = 2
        Float call_peak_mem_factor = 4.0
        Int call_peak_time_hr = 24
        Float call_peak_disk_factor = 30.0

        Float macs2_signal_track_mem_factor = 12.0
        Int macs2_signal_track_time_hr = 24
        Float macs2_signal_track_disk_factor = 80.0

        Float preseq_mem_factor = 0.5
        Float preseq_disk_factor = 5.0

        String? filter_picard_java_heap
        String? preseq_picard_java_heap
        String? fraglen_stat_picard_java_heap
        String? gc_bias_picard_java_heap
    }

    parameter_meta {
        docker: {
            description: 'Default Docker image URI to run WDL tasks.',
            group: 'runtime_environment',
            example: 'ubuntu:20.04'
        }
        singularity: {
            description: 'Default Singularity image URI to run WDL tasks. For Singularity users only.',
            group: 'runtime_environment',
            example: 'docker://ubuntu:20.04'
        }
        conda: {
            description: 'Default Conda environment name to run WDL tasks. For Conda users only.',
            group: 'runtime_environment',
            example: 'encd-atac'
        }
        conda_macs2: {
            description: 'Conda environment name for task macs2. For Conda users only.',
            group: 'runtime_environment',
            example: 'encd-atac-macs2'
        }
        conda_spp: {
            description: 'Conda environment name for tasks spp/xcor. For Conda users only.',
            group: 'runtime_environment',
            example: 'encd-atac-spp'
        }
        conda_python2: {
            description: 'Conda environment name for tasks with python2 wrappers (tss_enrich). For Conda users only.',
            group: 'runtime_environment',
            example: 'encd-atac-py2'
        }
        title: {
            description: 'Experiment title.',
            group: 'pipeline_metadata',
            example: 'ENCSR356KRQ (subsampled 1/400)'
        }
        description: {
            description: 'Experiment description.',
            group: 'pipeline_metadata',
            example: 'ATAC-seq on primary keratinocytes in day 0.0 of differentiation (subsampled 1/400)'
        }
        genome_tsv: {
            description: 'Reference genome database TSV.',
            group: 'reference_genome',
            help: 'This TSV files includes all genome specific parameters (e.g. reference FASTA, bowtie2 index). You can still invidiaully define any parameters in it. Parameters defined in input JSON will override those defined in genome TSV.',
            example: 'https://storage.googleapis.com/encode-pipeline-genome-data/genome_tsv/v1/hg38_caper.tsv'
        }
        genome_name: {
            description: 'Genome name.',
            group: 'reference_genome'
        }
        ref_fa: {
            description: 'Reference FASTA file.',
            group: 'reference_genome'
        }
        ref_mito_fa: {
            description: 'Reference FASTA file (mitochondrial reads only).',
            group: 'reference_genome'
        }
        bowtie2_idx_tar: {
            description: 'BWA index TAR file.',
            group: 'reference_genome'
        }
        bowtie2_mito_idx_tar: {
            description: 'BWA index TAR file (mitochondrial reads only).',
            group: 'reference_genome'
        }
        chrsz: {
            description: '2-col chromosome sizes file.',
            group: 'reference_genome'
        }
        blacklist: {
            description: 'Blacklist file in BED format.',
            group: 'reference_genome',
            help: 'Peaks will be filtered with this file.'
        }
        blacklist2: {
            description: 'Secondary blacklist file in BED format.',
            group: 'reference_genome',
            help: 'If it is defined, it will be merged with atac.blacklist. Peaks will be filtered with merged blacklist.'
        }
        mito_chr_name: {
            description: 'Mitochondrial chromosome name.',
            group: 'reference_genome',
            help: 'e.g. chrM, MT. Mitochondrial reads defined here will be filtered out during filtering BAMs in "filter" task.'
        }
        regex_bfilt_peak_chr_name: {
            description: 'Reg-ex for chromosomes to keep while filtering peaks.',
            group: 'reference_genome',
            help: 'Chromosomes defined here will be kept. All other chromosomes will be filtered out in .bfilt. peak file. This is done along with blacklist filtering peak file.'
        }
        gensz: {
            description: 'Genome sizes. "hs" for human, "mm" for mouse or sum of 2nd columnin chromosome sizes file.',
            group: 'reference_genome'
        }
        tss: {
            description: 'TSS file in BED format.',
            group: 'reference_genome'
        }
        dnase: {
            description: 'Open chromatin regions file in BED format.',
            group: 'reference_genome'
        }
        prom: {
            description: 'Promoter regions file in BED format.',
            group: 'reference_genome'
        }
        enh: {
            description: 'Enhancer regions file in BED format.',
            group: 'reference_genome'
        }
        reg2map: {
            description: 'Cell type signals file.',
            group: 'reference_genome'
        }
        reg2map_bed: {
            description: 'File of regions used to generate reg2map signals.',
            group: 'reference_genome'
        }
        roadmap_meta: {
            description: 'Roadmap metadata.',
            group: 'reference_genome'
        }
        paired_end: {
            description: 'Sequencing endedness.',
            group: 'input_genomic_data',
            help: 'Setting this on means that all replicates are paired ended. For mixed samples, use atac.paired_ends array instead.',
            example: true
        }
        paired_ends: {
            description: 'Sequencing endedness array (for mixed SE/PE datasets).',
            group: 'input_genomic_data',
            help: 'Whether each biological replicate is paired ended or not.'
        }
        fastqs_rep1_R1: {
            description: 'Read1 FASTQs to be merged for a biological replicate 1.',
            group: 'input_genomic_data',
            help: 'Define if you want to start pipeline from FASTQ files. Pipeline can start from any type of inputs (e.g. FASTQs, BAMs, ...). Choose one type and fill paramters for that type and leave other undefined. Especially for FASTQs, we have individual variable for each biological replicate to allow FASTQs of technical replicates can be merged. Make sure that they are consistent with read2 FASTQs (atac.fastqs_rep1_R2). These FASTQs are usually technical replicates to be merged.',
            example: [
                "https://storage.googleapis.com/encode-pipeline-test-samples/encode-atac-seq-pipeline/ENCSR356KRQ/fastq_subsampled/rep1/pair1/ENCFF341MYG.subsampled.400.fastq.gz",
                "https://storage.googleapis.com/encode-pipeline-test-samples/encode-atac-seq-pipeline/ENCSR356KRQ/fastq_subsampled/rep1/pair1/ENCFF106QGY.subsampled.400.fastq.gz"
            ]
        }
        fastqs_rep1_R2: {
            description: 'Read2 FASTQs to be merged for a biological replicate 1.',
            group: 'input_genomic_data',
            help: 'Make sure that they are consistent with read1 FASTQs (atac.fastqs_rep1_R1). These FASTQs are usually technical replicates to be merged.',
            example: [
                "https://storage.googleapis.com/encode-pipeline-test-samples/encode-atac-seq-pipeline/ENCSR356KRQ/fastq_subsampled/rep1/pair2/ENCFF248EJF.subsampled.400.fastq.gz",
                "https://storage.googleapis.com/encode-pipeline-test-samples/encode-atac-seq-pipeline/ENCSR356KRQ/fastq_subsampled/rep1/pair2/ENCFF368TYI.subsampled.400.fastq.gz"
            ]
        }
        fastqs_rep2_R1: {
            description: 'Read1 FASTQs to be merged for a biological replicate 2.',
            group: 'input_genomic_data',
            help: 'Make sure that they are consistent with read2 FASTQs (atac.fastqs_rep2_R2). These FASTQs are usually technical replicates to be merged.',
            example: [
                "https://storage.googleapis.com/encode-pipeline-test-samples/encode-atac-seq-pipeline/ENCSR356KRQ/fastq_subsampled/rep2/pair1/ENCFF641SFZ.subsampled.400.fastq.gz",
                "https://storage.googleapis.com/encode-pipeline-test-samples/encode-atac-seq-pipeline/ENCSR356KRQ/fastq_subsampled/rep2/pair1/ENCFF751XTV.subsampled.400.fastq.gz",
                "https://storage.googleapis.com/encode-pipeline-test-samples/encode-atac-seq-pipeline/ENCSR356KRQ/fastq_subsampled/rep2/pair1/ENCFF927LSG.subsampled.400.fastq.gz",
                "https://storage.googleapis.com/encode-pipeline-test-samples/encode-atac-seq-pipeline/ENCSR356KRQ/fastq_subsampled/rep2/pair1/ENCFF859BDM.subsampled.400.fastq.gz",
                "https://storage.googleapis.com/encode-pipeline-test-samples/encode-atac-seq-pipeline/ENCSR356KRQ/fastq_subsampled/rep2/pair1/ENCFF193RRC.subsampled.400.fastq.gz",
                "https://storage.googleapis.com/encode-pipeline-test-samples/encode-atac-seq-pipeline/ENCSR356KRQ/fastq_subsampled/rep2/pair1/ENCFF366DFI.subsampled.400.fastq.gz"
            ]
        }
        fastqs_rep2_R2: {
            description: 'Read2 FASTQs to be merged for a biological replicate 2.',
            group: 'input_genomic_data',
            help: 'Make sure that they are consistent with read1 FASTQs (atac.fastqs_rep2_R1). These FASTQs are usually technical replicates to be merged.',
            example: [
                "https://storage.googleapis.com/encode-pipeline-test-samples/encode-atac-seq-pipeline/ENCSR356KRQ/fastq_subsampled/rep2/pair2/ENCFF031ARQ.subsampled.400.fastq.gz",
                "https://storage.googleapis.com/encode-pipeline-test-samples/encode-atac-seq-pipeline/ENCSR356KRQ/fastq_subsampled/rep2/pair2/ENCFF590SYZ.subsampled.400.fastq.gz",
                "https://storage.googleapis.com/encode-pipeline-test-samples/encode-atac-seq-pipeline/ENCSR356KRQ/fastq_subsampled/rep2/pair2/ENCFF734PEQ.subsampled.400.fastq.gz",
                "https://storage.googleapis.com/encode-pipeline-test-samples/encode-atac-seq-pipeline/ENCSR356KRQ/fastq_subsampled/rep2/pair2/ENCFF007USV.subsampled.400.fastq.gz",
                "https://storage.googleapis.com/encode-pipeline-test-samples/encode-atac-seq-pipeline/ENCSR356KRQ/fastq_subsampled/rep2/pair2/ENCFF886FSC.subsampled.400.fastq.gz",
                "https://storage.googleapis.com/encode-pipeline-test-samples/encode-atac-seq-pipeline/ENCSR356KRQ/fastq_subsampled/rep2/pair2/ENCFF573UXK.subsampled.400.fastq.gz"
            ]
        }
        fastqs_rep3_R1: {
            description: 'Read1 FASTQs to be merged for a biological replicate 3.',
            group: 'input_genomic_data',
            help: 'Make sure that they are consistent with read2 FASTQs (atac.fastqs_rep3_R2). These FASTQs are usually technical replicates to be merged.'
        }
        fastqs_rep3_R2: {
            description: 'Read2 FASTQs to be merged for a biological replicate 3.',
            group: 'input_genomic_data',
            help: 'Make sure that they are consistent with read1 FASTQs (atac.fastqs_rep3_R1). These FASTQs are usually technical replicates to be merged.'
        }
        fastqs_rep4_R1: {
            description: 'Read1 FASTQs to be merged for a biological replicate 4.',
            group: 'input_genomic_data',
            help: 'Make sure that they are consistent with read2 FASTQs (atac.fastqs_rep4_R2). These FASTQs are usually technical replicates to be merged.'
        }
        fastqs_rep4_R2: {
            description: 'Read2 FASTQs to be merged for a biological replicate 4.',
            group: 'input_genomic_data',
            help: 'Make sure that they are consistent with read1 FASTQs (atac.fastqs_rep4_R1). These FASTQs are usually technical replicates to be merged.'
        }
        fastqs_rep5_R1: {
            description: 'Read1 FASTQs to be merged for a biological replicate 5.',
            group: 'input_genomic_data',
            help: 'Make sure that they are consistent with read2 FASTQs (atac.fastqs_rep5_R2). These FASTQs are usually technical replicates to be merged.'
        }
        fastqs_rep5_R2: {
            description: 'Read2 FASTQs to be merged for a biological replicate 5.',
            group: 'input_genomic_data',
            help: 'Make sure that they are consistent with read1 FASTQs (atac.fastqs_rep5_R1). These FASTQs are usually technical replicates to be merged.'
        }
        fastqs_rep6_R1: {
            description: 'Read1 FASTQs to be merged for a biological replicate 6.',
            group: 'input_genomic_data',
            help: 'Make sure that they are consistent with read2 FASTQs (atac.fastqs_rep6_R2). These FASTQs are usually technical replicates to be merged.'
        }
        fastqs_rep6_R2: {
            description: 'Read2 FASTQs to be merged for a biological replicate 6.',
            group: 'input_genomic_data',
            help: 'Make sure that they are consistent with read1 FASTQs (atac.fastqs_rep6_R1). These FASTQs are usually technical replicates to be merged.'
        }
        fastqs_rep7_R1: {
            description: 'Read1 FASTQs to be merged for a biological replicate 7.',
            group: 'input_genomic_data',
            help: 'Make sure that they are consistent with read2 FASTQs (atac.fastqs_rep7_R2). These FASTQs are usually technical replicates to be merged.'
        }
        fastqs_rep7_R2: {
            description: 'Read2 FASTQs to be merged for a biological replicate 7.',
            group: 'input_genomic_data',
            help: 'Make sure that they are consistent with read1 FASTQs (atac.fastqs_rep7_R1). These FASTQs are usually technical replicates to be merged.'
        }
        fastqs_rep8_R1: {
            description: 'Read1 FASTQs to be merged for a biological replicate 8.',
            group: 'input_genomic_data',
            help: 'Make sure that they are consistent with read2 FASTQs (atac.fastqs_rep8_R2). These FASTQs are usually technical replicates to be merged.'
        }
        fastqs_rep8_R2: {
            description: 'Read2 FASTQs to be merged for a biological replicate 8.',
            group: 'input_genomic_data',
            help: 'Make sure that they are consistent with read1 FASTQs (atac.fastqs_rep8_R1). These FASTQs are usually technical replicates to be merged.'
        }
        fastqs_rep9_R1: {
            description: 'Read1 FASTQs to be merged for a biological replicate 9.',
            group: 'input_genomic_data',
            help: 'Make sure that they are consistent with read2 FASTQs (atac.fastqs_rep9_R2). These FASTQs are usually technical replicates to be merged.'
        }
        fastqs_rep9_R2: {
            description: 'Read2 FASTQs to be merged for a biological replicate 9.',
            group: 'input_genomic_data',
            help: 'Make sure that they are consistent with read1 FASTQs (atac.fastqs_rep9_R1). These FASTQs are usually technical replicates to be merged.'
        }
        fastqs_rep10_R1: {
            description: 'Read1 FASTQs to be merged for a biological replicate 10.',
            group: 'input_genomic_data',
            help: 'Make sure that they are consistent with read2 FASTQs (atac.fastqs_rep10_R2). These FASTQs are usually technical replicates to be merged.'
        }
        fastqs_rep10_R2: {
            description: 'Read2 FASTQs to be merged for a biological replicate 10.',
            group: 'input_genomic_data',
            help: 'Make sure that they are consistent with read1 FASTQs (atac.fastqs_rep10_R1). These FASTQs are usually technical replicates to be merged.'
        }
        bams: {
            description: 'List of unfiltered/raw BAM files for each biological replicate.',
            group: 'input_genomic_data',
            help: 'Define if you want to start pipeline from BAM files. Unfiltered/raw BAM file generated from aligner (e.g. bowtie2). Each entry for each biological replicate. e.g. [rep1.bam, rep2.bam, rep3.bam, ...].'
        }
        nodup_bams: {
            description: 'List of filtered/deduped BAM files for each biological replicate',
            group: 'input_genomic_data',
            help: 'Define if you want to start pipeline from filtered BAM files. Filtered/deduped BAM file. Each entry for each biological replicate. e.g. [rep1.nodup.bam, rep2.nodup.bam, rep3.nodup.bam, ...].'
        }
        tas: {
            description: 'List of TAG-ALIGN files for each biological replicate.',
            group: 'input_genomic_data',
            help: 'Define if you want to start pipeline from TAG-ALIGN files. TAG-ALIGN is in a 6-col BED format. It is a simplified version of BAM. Each entry for each biological replicate. e.g. [rep1.tagAlign.gz, rep2.tagAlign.gz, ...].'
        }
        peaks: {
            description: 'List of NARROWPEAK files (not blacklist filtered) for each biological replicate.',
            group: 'input_genomic_data',
            help: 'Define if you want to start pipeline from PEAK files. Each entry for each biological replicate. e.g. [rep1.narrowPeak.gz, rep2.narrowPeak.gz, ...]. Define other PEAK parameters (e.g. atac.peaks_pr1, atac.peak_pooled) according to your flag settings (e.g. atac.true_rep_only) and number of replicates. If you have more than one replicate then define atac.peak_pooled, atac.peak_ppr1 and atac.peak_ppr2. If atac.true_rep_only flag is on then do not define any parameters (atac.peaks_pr1, atac.peaks_pr2, atac.peak_ppr1 and atac.peak_ppr2) related to pseudo replicates.'
        }
        peaks_pr1: {
            description: 'List of NARROWPEAK files (not blacklist filtered) for pseudo-replicate 1 of each biological replicate.',
            group: 'input_genomic_data',
            help: 'Define if you want to start pipeline from PEAK files. Define if atac.true_rep_only flag is off.'
        }
        peaks_pr2: {
            description: 'List of NARROWPEAK files (not blacklist filtered) for pseudo-replicate 2 of each biological replicate.',
            group: 'input_genomic_data',
            help: 'Define if you want to start pipeline from PEAK files. Define if atac.true_rep_only flag is off.'
        }
        peak_pooled: {
            description: 'NARROWPEAK file for pooled true replicate.',
            group: 'input_genomic_data',
            help: 'Define if you want to start pipeline from PEAK files. Define if you have multiple biological replicates. Pooled true replicate means analysis on pooled biological replicates.'
        }
        peak_ppr1: {
            description: 'NARROWPEAK file for pooled pseudo replicate 1.',
            group: 'input_genomic_data',
            help: 'Define if you want to start pipeline from PEAK files. Define if you have multiple biological replicates and atac.true_rep_only flag is off. PPR1 means analysis on pooled 1st pseudo replicates. Each biological replicate is shuf/split into two pseudos. This is a pooling of each replicate\'s 1st pseudos.'
        }
        peak_ppr2: {
            description: 'NARROWPEAK file for pooled pseudo replicate 2.',
            group: 'input_genomic_data',
            help: 'Define if you want to start pipeline from PEAK files. Define if you have multiple biological replicates and atac.true_rep_only flag is off. PPR1 means analysis on pooled 2nd pseudo replicates. Each biological replicate is shuf/split into two pseudos. This is a pooling of each replicate\'s 2nd pseudos.'
        }

        pipeline_type: {
            description: 'Pipeline type. atac for ATAC-Seq or dnase for DNase-Seq.',
            group: 'pipeline_parameter',
            help: 'The only difference of two types is that TN5 shifting of TAG-ALIGN is done for atac. TAG-ALIGN is in 6-col BED format. It is a simplified version of BAM.',
            choices: ['atac', 'dnase'],
            example: 'atac'
        }
        align_only: {
            description: 'Align only mode.',
            group: 'pipeline_parameter',
            help: 'Reads will be aligned but there will be no peak-calling on them.'
        }
        true_rep_only: {
            description: 'Disables all analyses related to pseudo-replicates.',
            group: 'pipeline_parameter',
            help: 'Pipeline generates 2 pseudo-replicate from one biological replicate. This flag turns off all analyses related to pseudos (with prefix/suffix pr, ppr).'
        }
        enable_xcor: {
            description: 'Enables cross-correlation analysis.',
            group: 'pipeline_parameter',
            help: 'Generates cross-correlation plot.'
        }
        enable_count_signal_track: {
            description: 'Enables generation of count signal tracks.',
            group: 'pipeline_parameter'
        }
        enable_idr: {
            description: 'Enables IDR on MACS2 NARROWPEAKs.',
            group: 'pipeline_parameter'
        }
        enable_preseq: {
            description: 'Enables preseq analysis.',
            group: 'pipeline_parameter'
        }
        enable_fraglen_stat: {
            description: 'Enables calculation of fragment length distribution/statistics.',
            group: 'pipeline_parameter'
        }
        enable_tss_enrich: {
            description: 'Enables TSS enrichment plot generation.',
            group: 'pipeline_parameter'
        }
        enable_annot_enrich: {
            description: 'Enables annotated regions enrichment analysis.',
            group: 'pipeline_parameter'
        }
        enable_jsd: {
            description: 'Enables Jensen-Shannon Distance (JSD) plot generation.',
            group: 'pipeline_parameter'
        }
        enable_compare_to_roadmap: {
            description: 'Enables comparison to Roadmap.',
            group: 'pipeline_parameter'
        }
        enable_gc_bias: {
            description: 'Enables GC bias calculation.',
            group: 'pipeline_parameter'
        }

        cutadapt_param: {
            description: 'Parameters for cutadapt.',
            group: 'adapter_trimming',
            help: 'It is -e 0.1 -m 5 by default (err_rate=0.1, min_trim_len=5). You can define any parameters that cutadapt supports.'
        }
        auto_detect_adapter: {
            description: 'Auto-detect/trim adapter sequences.',
            group: 'adapter_trimming',
            help: 'Can detect/trim three types of adapter sequences. Illumina: AGATCGGAAGAGC, Nextera: CTGTCTCTTATA, smallRNA: TGGAATTCTCGG.',
            example: true
        }
        adapter: {
            description: 'Adapter for all FASTQs.',
            group: 'adapter_trimming',
            help: 'Define if all FASTQs have the same adapter sequence. Otherwise define adapter sequence for individual FASTQ in atac.adapters_repX_R1 and atac.adapters_repX_R2 instead. Use atac.auto_detect_adapter if you want to detect adapters automatically. If all of your FASTQs are already trimmed then leave all adapter-related parameters undefined/empty.'
        }
        adapters_rep1_R1: {
            description: 'Adapter sequences for read1 FASTQs to be merged for a biological replicate 1.',
            group: 'adapter_trimming',
            help: 'Make sure that they are consistent with read2 FASTQs (atac.adapters_rep1_R2). You can combine this with atac.auto_detect_adapter. Pipeline will auto-detect/trim adapter sequences for null entry in this list. e.g. ["AAGGCCTT", null, "AAGGCCTT"].'
        }
        adapters_rep1_R2: {
            description: 'Adapter sequences for read2 FASTQs to be merged for a biological replicate 1.',
            group: 'adapter_trimming',
            help: 'Make sure that they are consistent with read1 FASTQs (atac.adapters_rep1_R1).'
        }
        adapters_rep2_R1: {
            description: 'Adapter sequences for read1 FASTQs to be merged for a biological replicate 2.',
            group: 'adapter_trimming',
            help: 'Make sure that they are consistent with read2 FASTQs (atac.adapters_rep2_R2).'
        }
        adapters_rep2_R2: {
            description: 'Adapter sequences for read2 FASTQs to be merged for a biological replicate 2.',
            group: 'adapter_trimming',
            help: 'Make sure that they are consistent with read1 FASTQs (atac.adapters_rep2_R1).'
        }
        adapters_rep3_R1: {
            description: 'Adapter sequences for read1 FASTQs to be merged for a biological replicate 3.',
            group: 'adapter_trimming',
            help: 'Make sure that they are consistent with read2 FASTQs (atac.adapters_rep3_R2).'
        }
        adapters_rep3_R2: {
            description: 'Adapter sequences for read2 FASTQs to be merged for a biological replicate 3.',
            group: 'adapter_trimming',
            help: 'Make sure that they are consistent with read1 FASTQs (atac.adapters_rep3_R1).'
        }
        adapters_rep4_R1: {
            description: 'Adapter sequences for read1 FASTQs to be merged for a biological replicate 4.',
            group: 'adapter_trimming',
            help: 'Make sure that they are consistent with read2 FASTQs (atac.adapters_rep4_R2).'
        }
        adapters_rep4_R2: {
            description: 'Adapter sequences for read2 FASTQs to be merged for a biological replicate 4.',
            group: 'adapter_trimming',
            help: 'Make sure that they are consistent with read1 FASTQs (atac.adapters_rep4_R1).'
        }
        adapters_rep5_R1: {
            description: 'Adapter sequences for read1 FASTQs to be merged for a biological replicate 5.',
            group: 'adapter_trimming',
            help: 'Make sure that they are consistent with read2 FASTQs (atac.adapters_rep5_R2).'
        }
        adapters_rep5_R2: {
            description: 'Adapter sequences for read2 FASTQs to be merged for a biological replicate 5.',
            group: 'adapter_trimming',
            help: 'Make sure that they are consistent with read1 FASTQs (atac.adapters_rep5_R1).'
        }
        adapters_rep6_R1: {
            description: 'Adapter sequences for read1 FASTQs to be merged for a biological replicate 6.',
            group: 'adapter_trimming',
            help: 'Make sure that they are consistent with read2 FASTQs (atac.adapters_rep6_R2).'
        }
        adapters_rep6_R2: {
            description: 'Adapter sequences for read2 FASTQs to be merged for a biological replicate 6.',
            group: 'adapter_trimming',
            help: 'Make sure that they are consistent with read1 FASTQs (atac.adapters_rep6_R1).'
        }
        adapters_rep7_R1: {
            description: 'Adapter sequences for read1 FASTQs to be merged for a biological replicate 7.',
            group: 'adapter_trimming',
            help: 'Make sure that they are consistent with read2 FASTQs (atac.adapters_rep7_R2).'
        }
        adapters_rep7_R2: {
            description: 'Adapter sequences for read2 FASTQs to be merged for a biological replicate 7.',
            group: 'adapter_trimming',
            help: 'Make sure that they are consistent with read1 FASTQs (atac.adapters_rep7_R1).'
        }
        adapters_rep8_R1: {
            description: 'Adapter sequences for read1 FASTQs to be merged for a biological replicate 8.',
            group: 'adapter_trimming',
            help: 'Make sure that they are consistent with read2 FASTQs (atac.adapters_rep8_R2).'
        }
        adapters_rep8_R2: {
            description: 'Adapter sequences for read2 FASTQs to be merged for a biological replicate 8.',
            group: 'adapter_trimming',
            help: 'Make sure that they are consistent with read1 FASTQs (atac.adapters_rep8_R1).'
        }
        adapters_rep9_R1: {
            description: 'Adapter sequences for read1 FASTQs to be merged for a biological replicate 9.',
            group: 'adapter_trimming',
            help: 'Make sure that they are consistent with read2 FASTQs (atac.adapters_rep9_R2).'
        }
        adapters_rep9_R2: {
            description: 'Adapter sequences for read2 FASTQs to be merged for a biological replicate 9.',
            group: 'adapter_trimming',
            help: 'Make sure that they are consistent with read1 FASTQs (atac.adapters_rep9_R1).'
        }
        adapters_rep10_R1: {
            description: 'Adapter sequences for read1 FASTQs to be merged for a biological replicate 10.',
            group: 'adapter_trimming',
            help: 'Make sure that they are consistent with read2 FASTQs (atac.adapters_rep10_R2).'
        }
        adapters_rep10_R2: {
            description: 'Adapter sequences for read2 FASTQs to be merged for a biological replicate 10.',
            group: 'adapter_trimming',
            help: 'Make sure that they are consistent with read1 FASTQs (atac.adapters_rep10_R1).'
        }

        multimapping: {
            description: 'Number of multimappers.',
            group: 'alignment',
            help: 'It is 4 by default. Set it to 0 if your sample does not have multimappers.'
        }
        dup_marker: {
            description: 'Marker for duplicate reads. picard or sambamba.',
            group: 'alignment',
            help: 'picard for Picard MarkDuplicates or sambamba for sambamba markdup.'
        }
        no_dup_removal: {
            description: 'Disable removal of duplicate reads during filtering BAM.',
            group: 'alignment',
            help: 'Duplicate reads are filtererd out during filtering BAMs to gerenate NODUP_BAM. This flag will keep all duplicate reads in NODUP_BAM. This flag does not affect naming of NODUP_BAM. NODUP_BAM will still have .nodup. suffix in its filename.'
        }
        mapq_thresh: {
            description: 'Threshold for low MAPQ reads removal.',
            group: 'alignment',
            help: 'Low MAPQ reads are filtered out while filtering BAM.'
        }
        filter_chrs: {
            description: 'List of chromosomes to be filtered out while filtering BAM.',
            group: 'alignment',
            help: 'It is ["chrM", "MT"] by default. Therefore, mitochondrial reads will be filtered out while filtering. Make it empty if you want to keep all reads.'
        }
        subsample_reads: {
            description: 'Subsample reads. Shuffle and subsample reads.',
            group: 'alignment',
            help: 'This affects all downstream analyses after filtering BAM. (e.g. all TAG-ALIGN files, peak-calling). Reads will be shuffled only if actual number of reads in BAM exceeds this number.  0 means disabled.'
        }
        xcor_subsample_reads: {
            description: 'Subsample reads for cross-corrlelation analysis only.',
            group: 'alignment',
            help: 'This does not affect downstream analyses after filtering BAM. It is for cross-correlation analysis only. 0 means disabled.'
        }
        read_len: {
            description: 'Read length per biological replicate.',
            group: 'alignment',
            help: 'Pipeline can estimate read length from FASTQs. If you start pipeline from other types (BAM, NODUP_BAM, TA, ...) than FASTQ. Then provide this for some analyses that require read length (e.g. TSS enrichment plot).'
        }
        pseudoreplication_random_seed: {
            description: 'Random seed (positive integer) used for pseudo-replication (shuffling reads in TAG-ALIGN and then split it into two).',
            group: 'alignment',
            help: 'Pseudo-replication (task spr) is done by using GNU "shuf --random-source=sha256(random_seed)". If this parameter == 0, then pipeline uses input TAG-ALIGN file\'s size (in bytes) for the random_seed.'
        }
        cap_num_peak: {
            description: 'Upper limit on the number of peaks.',
            group: 'peak_calling',
            help: 'Called peaks will be sorted in descending order of score and the number of peaks will be capped at this number by taking first N peaks.'
        }
        pval_thresh: {
            description: 'p-value Threshold for MACS2 peak caller.',
            group: 'peak_calling',
            help: 'macs2 callpeak -p'
        }
        smooth_win: {
            description: 'Size of smoothing windows for MACS2 peak caller.',
            group: 'peak_calling',
            help: 'This will be used for both generating MACS2 peaks/signal tracks.'
        }
        idr_thresh: {
            description: 'IDR threshold.',
            group: 'peak_calling'
        }
        align_cpu: {
            description: 'Number of cores for task align.',
            group: 'resource_parameter',
            help: 'Task align merges/crops/maps FASTQs.'
        }
        align_mem_factor: {
            description: 'Multiplication factor to determine memory required for task align.',
            group: 'resource_parameter',
            help: 'This factor will be multiplied to the size of FASTQs to determine required memory of instance (GCP/AWS) or job (HPCs).'
        }
        align_time_hr: {
            description: 'Walltime (h) required for task align.',
            group: 'resource_parameter',
            help: 'This is for HPCs only. e.g. SLURM, SGE, ...'
        }
        align_disk_factor: {
            description: 'Multiplication factor to determine persistent disk size for task align.',
            group: 'resource_parameter',
            help: 'This factor will be multiplied to the size of FASTQs to determine required disk size of instance on GCP/AWS.'
        }
        filter_cpu: {
            description: 'Number of cores for task filter.',
            group: 'resource_parameter',
            help: 'Task filter filters raw/unfilterd BAM to get filtered/deduped BAM.'
        }
        filter_mem_factor: {
            description: 'Multiplication factor to determine memory required for task filter.',
            group: 'resource_parameter',
            help: 'This factor will be multiplied to the size of BAMs to determine required memory of instance (GCP/AWS) or job (HPCs).'
        }
        filter_time_hr: {
            description: 'Walltime (h) required for task filter.',
            group: 'resource_parameter',
            help: 'This is for HPCs only. e.g. SLURM, SGE, ...'
        }
        filter_disk_factor: {
            description: 'Multiplication factor to determine persistent disk size for task filter.',
            group: 'resource_parameter',
            help: 'This factor will be multiplied to the size of BAMs to determine required disk size of instance on GCP/AWS.'
        }
        bam2ta_cpu: {
            description: 'Number of cores for task bam2ta.',
            group: 'resource_parameter',
            help: 'Task bam2ta converts filtered/deduped BAM in to TAG-ALIGN (6-col BED) format.'
        }
        bam2ta_mem_factor: {
            description: 'Multiplication factor to determine memory required for task bam2ta.',
            group: 'resource_parameter',
            help: 'This factor will be multiplied to the size of filtered BAMs to determine required memory of instance (GCP/AWS) or job (HPCs).'
        }
        bam2ta_time_hr: {
            description: 'Walltime (h) required for task bam2ta.',
            group: 'resource_parameter',
            help: 'This is for HPCs only. e.g. SLURM, SGE, ...'
        }
        bam2ta_disk_factor: {
            description: 'Multiplication factor to determine persistent disk size for task bam2ta.',
            group: 'resource_parameter',
            help: 'This factor will be multiplied to the size of filtered BAMs to determine required disk size of instance on GCP/AWS.'
        }
        spr_mem_factor: {
            description: 'Multiplication factor to determine memory required for task spr.',
            group: 'resource_parameter',
            help: 'This factor will be multiplied to the size of filtered BAMs to determine required memory of instance (GCP/AWS) or job (HPCs).'
        }
        spr_disk_factor: {
            description: 'Multiplication factor to determine persistent disk size for task spr.',
            group: 'resource_parameter',
            help: 'This factor will be multiplied to the size of filtered BAMs to determine required disk size of instance on GCP/AWS.'
        }
        jsd_cpu: {
            description: 'Number of cores for task jsd.',
            group: 'resource_parameter',
            help: 'Task jsd plots Jensen-Shannon distance and metrics related to it.'
        }
        jsd_mem_factor: {
            description: 'Multiplication factor to determine memory required for task jsd.',
            group: 'resource_parameter',
            help: 'This factor will be multiplied to the size of filtered BAMs to determine required memory of instance (GCP/AWS) or job (HPCs).'
        }
        jsd_time_hr: {
            description: 'Walltime (h) required for task jsd.',
            group: 'resource_parameter',
            help: 'This is for HPCs only. e.g. SLURM, SGE, ...'
        }
        jsd_disk_factor: {
            description: 'Multiplication factor to determine persistent disk size for task jsd.',
            group: 'resource_parameter',
            help: 'This factor will be multiplied to the size of filtered BAMs to determine required disk size of instance on GCP/AWS.'
        }
        xcor_cpu: {
            description: 'Number of cores for task xcor.',
            group: 'resource_parameter',
            help: 'Task xcor does cross-correlation analysis (including a plot) on subsampled TAG-ALIGNs.'
        }
        xcor_mem_factor: {
            description: 'Multiplication factor to determine memory required for task xcor.',
            group: 'resource_parameter',
            help: 'This factor will be multiplied to the size of TAG-ALIGNs (BEDs) to determine required memory of instance (GCP/AWS) or job (HPCs).'
        }
        xcor_time_hr: {
            description: 'Walltime (h) required for task xcor.',
            group: 'resource_parameter',
            help: 'This is for HPCs only. e.g. SLURM, SGE, ...'
        }
        xcor_disk_factor: {
            description: 'Multiplication factor to determine persistent disk size for task xcor.',
            group: 'resource_parameter',
            help: 'This factor will be multiplied to the size of TAG-ALIGNs (BEDs) to determine required disk size of instance on GCP/AWS.'
        }
        call_peak_cpu: {
            description: 'Number of cores for task call_peak. MACS2 is single-thread. No more than 2 is required.',
            group: 'resource_parameter',
            help: 'Task call_peak call peaks on TAG-ALIGNs by using MACS2 peak caller.'
        }
        call_peak_mem_factor: {
            description: 'Multiplication factor to determine memory required for task call_peak.',
            group: 'resource_parameter',
            help: 'This factor will be multiplied to the size of TAG-ALIGNs (BEDs) to determine required memory of instance (GCP/AWS) or job (HPCs).'
        }
        call_peak_time_hr: {
            description: 'Walltime (h) required for task call_peak.',
            group: 'resource_parameter',
            help: 'This is for HPCs only. e.g. SLURM, SGE, ...'
        }
        call_peak_disk_factor: {
            description: 'Multiplication factor to determine persistent disk size for task call_peak.',
            group: 'resource_parameter',
            help: 'This factor will be multiplied to the size of TAG-ALIGNs (BEDs) to determine required disk size of instance on GCP/AWS.'
        }
        macs2_signal_track_mem_factor: {
            description: 'Multiplication factor to determine memory required for task macs2_signal_track.',
            group: 'resource_parameter',
            help: 'This factor will be multiplied to the size of TAG-ALIGNs (BEDs) to determine required memory of instance (GCP/AWS) or job (HPCs).'
        }
        macs2_signal_track_time_hr: {
            description: 'Walltime (h) required for task macs2_signal_track.',
            group: 'resource_parameter',
            help: 'This is for HPCs only. e.g. SLURM, SGE, ...'
        }
        macs2_signal_track_disk_factor: {
            description: 'Multiplication factor to determine persistent disk size for task macs2_signal_track.',
            group: 'resource_parameter',
            help: 'This factor will be multiplied to the size of TAG-ALIGNs (BEDs) to determine required disk size of instance on GCP/AWS.'
        }
        preseq_mem_factor: {
            description: 'Multiplication factor to determine memory required for task preseq.',
            group: 'resource_parameter',
            help: 'This factor will be multiplied to the size of BAMs to determine required memory of instance (GCP/AWS) or job (HPCs).'
        }
        preseq_disk_factor: {
            description: 'Multiplication factor to determine persistent disk size for task preseq.',
            group: 'resource_parameter',
            help: 'This factor will be multiplied to the size of BAMs to determine required disk size of instance on GCP/AWS.'
        }
        filter_picard_java_heap: {
            description: 'Maximum Java heap (java -Xmx) in task filter.',
            group: 'resource_parameter',
            help: 'Maximum memory for Picard tools MarkDuplicates. If not defined, 90% of filter task\'s memory will be used.'
        }
        preseq_picard_java_heap: {
            description: 'Maximum Java heap (java -Xmx) in task preseq.',
            group: 'resource_parameter',
            help: 'Maximum memory for Picard tools EstimateLibraryComplexity. If not defined, 90% of preseq task\'s memory will be used.'
        }
        fraglen_stat_picard_java_heap: {
            description: 'Maximum Java heap (java -Xmx) in task fraglen_stat_pe (for paired end replicate only).',
            group: 'resource_parameter',
            help: 'Maximum memory for Picard tools CollectInsertSizeMetrics. If not defined, 90% of fraglen_stat tasks\'s memory will be used.'
        }
        gc_bias_picard_java_heap: {
            description: 'Maximum Java heap (java -Xmx) in task gc_bias.',
            group: 'resource_parameter',
            help: 'Maximum memory for Picard tools CollectGcBiasMetrics. If not defined, 90% of gc_bias tasks\'s memory will be used.'
        }
    }
    RuntimeEnvironment runtime_environment = {
        'docker': docker, 'singularity': singularity, 'conda': conda
    }
    RuntimeEnvironment runtime_environment_spp = {
        'docker': docker, 'singularity': singularity, 'conda': conda_spp
    }
    RuntimeEnvironment runtime_environment_macs2 = {
        'docker': docker, 'singularity': singularity, 'conda': conda_macs2
    }
    RuntimeEnvironment runtime_environment_python2 = {
        'docker': docker, 'singularity': singularity, 'conda': conda_python2
    }

    String aligner = 'bowtie2'
    String peak_caller = 'macs2'
    String peak_type = 'narrowPeak'
    
    # read genome data and paths
    if ( defined(genome_tsv) ) {
        call read_genome_tsv { input:
            genome_tsv = genome_tsv,
            runtime_environment = runtime_environment
        }
    }
    File ref_fa_ = select_first([ref_fa, read_genome_tsv.ref_fa])
    File bowtie2_idx_tar_ = select_first([bowtie2_idx_tar, read_genome_tsv.bowtie2_idx_tar])
    File bowtie2_mito_idx_tar_ = select_first([bowtie2_mito_idx_tar, read_genome_tsv.bowtie2_mito_idx_tar])
    File chrsz_ = select_first([chrsz, read_genome_tsv.chrsz])
    String gensz_ = select_first([gensz, read_genome_tsv.gensz])
    File? blacklist1_ = if defined(blacklist) then blacklist
        else read_genome_tsv.blacklist
    File? blacklist2_ = if defined(blacklist2) then blacklist2
        else read_genome_tsv.blacklist2        
    # merge multiple blacklists
    # two blacklists can have different number of columns (3 vs 6)
    # so we limit merged blacklist's columns to 3
    Array[File] blacklists = select_all([blacklist1_, blacklist2_])
    if ( length(blacklists) > 1 ) {
        call pool_ta as pool_blacklist { input:
            tas = blacklists,
            col = 3,
            runtime_environment = runtime_environment
        }
    }
    File? blacklist_ = if length(blacklists) > 1 then pool_blacklist.ta_pooled
        else if length(blacklists) > 0 then blacklists[0]
        else blacklist2_
    String mito_chr_name_ = select_first([mito_chr_name, read_genome_tsv.mito_chr_name])
    String regex_bfilt_peak_chr_name_ = select_first([regex_bfilt_peak_chr_name, read_genome_tsv.regex_bfilt_peak_chr_name])
    String genome_name_ = select_first([genome_name, read_genome_tsv.genome_name, basename(chrsz_)])

    # read additional annotation data
    File? tss_ = if defined(tss) then tss
        else read_genome_tsv.tss
    File? dnase_ = if defined(dnase) then dnase
        else read_genome_tsv.dnase
    File? prom_ = if defined(prom) then prom
        else read_genome_tsv.prom
    File? enh_ = if defined(enh) then enh
        else read_genome_tsv.enh
    File? reg2map_ = if defined(reg2map) then reg2map
        else read_genome_tsv.reg2map
    File? reg2map_bed_ = if defined(reg2map_bed) then reg2map_bed
        else read_genome_tsv.reg2map_bed
    File? roadmap_meta_ = if defined(roadmap_meta) then roadmap_meta
        else read_genome_tsv.roadmap_meta

    ### temp vars (do not define these)
    String aligner_ = aligner
    String peak_caller_ = peak_caller
    String peak_type_ = peak_type
    String idr_rank_ = if peak_caller_=='spp' then 'signal.value'
                        else if peak_caller_=='macs2' then 'p.value'
                        else 'p.value'
    Int cap_num_peak_ = cap_num_peak
    Int mapq_thresh_ = mapq_thresh

    # temporary 2-dim fastqs array [rep_id][merge_id]
    Array[Array[File]] fastqs_R1 = 
        if length(fastqs_rep10_R1)>0 then
            [fastqs_rep1_R1, fastqs_rep2_R1, fastqs_rep3_R1, fastqs_rep4_R1, fastqs_rep5_R1,
            fastqs_rep6_R1, fastqs_rep7_R1, fastqs_rep8_R1, fastqs_rep9_R1, fastqs_rep10_R1]
        else if length(fastqs_rep9_R1)>0 then
            [fastqs_rep1_R1, fastqs_rep2_R1, fastqs_rep3_R1, fastqs_rep4_R1, fastqs_rep5_R1,
            fastqs_rep6_R1, fastqs_rep7_R1, fastqs_rep8_R1, fastqs_rep9_R1]
        else if length(fastqs_rep8_R1)>0 then
            [fastqs_rep1_R1, fastqs_rep2_R1, fastqs_rep3_R1, fastqs_rep4_R1, fastqs_rep5_R1,
            fastqs_rep6_R1, fastqs_rep7_R1, fastqs_rep8_R1]
        else if length(fastqs_rep7_R1)>0 then
            [fastqs_rep1_R1, fastqs_rep2_R1, fastqs_rep3_R1, fastqs_rep4_R1, fastqs_rep5_R1,
            fastqs_rep6_R1, fastqs_rep7_R1]
        else if length(fastqs_rep6_R1)>0 then
            [fastqs_rep1_R1, fastqs_rep2_R1, fastqs_rep3_R1, fastqs_rep4_R1, fastqs_rep5_R1,
            fastqs_rep6_R1]
        else if length(fastqs_rep5_R1)>0 then
            [fastqs_rep1_R1, fastqs_rep2_R1, fastqs_rep3_R1, fastqs_rep4_R1, fastqs_rep5_R1]
        else if length(fastqs_rep4_R1)>0 then
            [fastqs_rep1_R1, fastqs_rep2_R1, fastqs_rep3_R1, fastqs_rep4_R1]
        else if length(fastqs_rep3_R1)>0 then
            [fastqs_rep1_R1, fastqs_rep2_R1, fastqs_rep3_R1]
        else if length(fastqs_rep2_R1)>0 then
            [fastqs_rep1_R1, fastqs_rep2_R1]
        else if length(fastqs_rep1_R1)>0 then
            [fastqs_rep1_R1]
        else []
    # no need to do that for R2 (R1 array will be used to determine presense of fastq for each rep)
    Array[Array[File]] fastqs_R2 = 
        [fastqs_rep1_R2, fastqs_rep2_R2, fastqs_rep3_R2, fastqs_rep4_R2, fastqs_rep5_R2,
        fastqs_rep6_R2, fastqs_rep7_R2, fastqs_rep8_R2, fastqs_rep9_R2, fastqs_rep10_R2]

    # temporary 2-dim adapters array [rep_id][merge_id]
    Array[Array[String]] adapters_R1 = 
        [adapters_rep1_R1, adapters_rep2_R1, adapters_rep3_R1, adapters_rep4_R1, adapters_rep5_R1,
        adapters_rep6_R1, adapters_rep7_R1, adapters_rep8_R1, adapters_rep9_R1, adapters_rep10_R1]
    Array[Array[String]] adapters_R2 = 
        [adapters_rep1_R2, adapters_rep2_R2, adapters_rep3_R2, adapters_rep4_R2, adapters_rep5_R2,
        adapters_rep6_R2, adapters_rep7_R2, adapters_rep8_R2, adapters_rep9_R2, adapters_rep10_R2]

    # temporary variables to get number of replicates
    #       WDLic implementation of max(A,B,C,...)
    Int num_rep_fastq = length(fastqs_R1)
    Int num_rep_bam = if length(bams)<num_rep_fastq then num_rep_fastq
        else length(bams)
    Int num_rep_nodup_bam = if length(nodup_bams)<num_rep_bam then num_rep_bam
        else length(nodup_bams)
    Int num_rep_ta = if length(tas)<num_rep_nodup_bam then num_rep_nodup_bam
        else length(tas)
    Int num_rep_peak = if length(peaks)<num_rep_ta then num_rep_ta
        else length(peaks)
    Int num_rep = num_rep_peak

    # sanity check for inputs
    if ( num_rep == 0 ) {
        call raise_exception as error_input_data  { input:
            msg = 'No FASTQ/BAM/TAG-ALIGN/PEAK defined in your input JSON. Check if your FASTQs are defined as "atac.fastqs_repX_RY". DO NOT MISS suffix _R1 even for single ended FASTQ.',
            runtime_environment = runtime_environment
        }
    }

    # align each replicate
    scatter(i in range(num_rep)) {
        # to override endedness definition for individual replicate
        #     paired_end will override paired_ends[i]
        Boolean paired_end_ = if !defined(paired_end) && i<length(paired_ends) then paired_ends[i]
            else select_first([paired_end])

        Boolean has_input_of_align = i<length(fastqs_R1) && length(fastqs_R1[i])>0
        Boolean has_output_of_align = i<length(bams)
        if ( has_input_of_align && !has_output_of_align ) {
            call align { input :
                fastqs_R1 = fastqs_R1[i],
                fastqs_R2 = fastqs_R2[i],
                adapter = adapter,
                adapters_R1 = adapters_R1[i],
                adapters_R2 = adapters_R2[i],
                paired_end = paired_end_,
                auto_detect_adapter = auto_detect_adapter,
                cutadapt_param = cutadapt_param,

                aligner = aligner_,
                mito_chr_name = mito_chr_name_,
                chrsz = chrsz_,
                multimapping = multimapping,
                idx_tar = bowtie2_idx_tar_,
                # resource
                cpu = align_cpu,
                mem_factor = align_mem_factor,
                time_hr = align_time_hr,
                disk_factor = align_disk_factor,
                runtime_environment = runtime_environment
            }
        }
        File? bam_ = if has_output_of_align then bams[i] else align.bam

        # mito only mapping to get frac mito qc
        Boolean has_input_of_align_mito = has_input_of_align &&
            defined(bowtie2_mito_idx_tar_)
        if ( has_input_of_align_mito ) {
            call align as align_mito { input :
                fastqs_R1 = fastqs_R1[i],
                fastqs_R2 = fastqs_R2[i],
                adapter = adapter,
                adapters_R1 = adapters_R1[i],
                adapters_R2 = adapters_R2[i],
                paired_end = paired_end_,
                auto_detect_adapter = auto_detect_adapter,
                cutadapt_param = cutadapt_param,

                aligner = aligner_,
                mito_chr_name = mito_chr_name_,
                chrsz = chrsz_,
                multimapping = multimapping,
                idx_tar = bowtie2_mito_idx_tar_,
                # resource
                cpu = align_cpu,
                mem_factor = align_mem_factor,
                time_hr = align_time_hr,
                disk_factor = align_disk_factor,
                runtime_environment = runtime_environment
            }
        }

        if ( defined(align.non_mito_samstat_qc) && defined(align_mito.samstat_qc) ) {
            call frac_mito { input:
                non_mito_samstat = align.non_mito_samstat_qc,
                mito_samstat = align_mito.samstat_qc,
                runtime_environment = runtime_environment
            }
        }

        Boolean has_input_of_filter = has_output_of_align || defined(align.bam)
        Boolean has_output_of_filter = i<length(nodup_bams)
        # skip if we already have output of this step
        if ( has_input_of_filter && !has_output_of_filter ) {
            call filter { input :
                bam = bam_,
                paired_end = paired_end_,
                dup_marker = dup_marker,
                mapq_thresh = mapq_thresh_,
                filter_chrs = filter_chrs,
                chrsz = chrsz_,
                no_dup_removal = no_dup_removal,
                multimapping = multimapping,
                mito_chr_name = mito_chr_name_,

                cpu = filter_cpu,
                mem_factor = filter_mem_factor,
                picard_java_heap = filter_picard_java_heap,
                time_hr = filter_time_hr,
                disk_factor = filter_disk_factor,
                runtime_environment = runtime_environment
            }
        }
        File? nodup_bam_ = if has_output_of_filter then nodup_bams[i] else filter.nodup_bam

        Boolean has_input_of_bam2ta = has_output_of_filter || defined(filter.nodup_bam)
        Boolean has_output_of_bam2ta = i<length(tas)
        if ( has_input_of_bam2ta && !has_output_of_bam2ta ) {
            call bam2ta { input :
                bam = nodup_bam_,
                disable_tn5_shift = if pipeline_type=='atac' then false else true,
                subsample = subsample_reads,
                paired_end = paired_end_,
                mito_chr_name = mito_chr_name_,

                cpu = bam2ta_cpu,
                mem_factor = bam2ta_mem_factor,
                time_hr = bam2ta_time_hr,
                disk_factor = bam2ta_disk_factor,
                runtime_environment = runtime_environment
            }
        }
        File? ta_ = if has_output_of_bam2ta then tas[i] else bam2ta.ta

        Boolean has_input_of_xcor = has_output_of_align || defined(align.bam)
        if ( has_input_of_xcor && enable_xcor ) {
            call filter as filter_no_dedup { input :
                bam = bam_,
                paired_end = paired_end_,
                dup_marker = dup_marker,
                mapq_thresh = mapq_thresh_,
                filter_chrs = filter_chrs,
                chrsz = chrsz_,
                no_dup_removal = true,
                multimapping = multimapping,
                mito_chr_name = mito_chr_name_,

                cpu = filter_cpu,
                mem_factor = filter_mem_factor,
                picard_java_heap = filter_picard_java_heap,
                time_hr = filter_time_hr,
                disk_factor = filter_disk_factor,
                runtime_environment = runtime_environment
            }
            call bam2ta as bam2ta_no_dedup { input :
                bam = filter_no_dedup.nodup_bam,  # output name is nodup but it's not deduped
                disable_tn5_shift = if pipeline_type=='atac' then false else true,
                subsample = 0,
                paired_end = paired_end_,
                mito_chr_name = mito_chr_name_,

                cpu = bam2ta_cpu,
                mem_factor = bam2ta_mem_factor,
                time_hr = bam2ta_time_hr,
                disk_factor = bam2ta_disk_factor,
                runtime_environment = runtime_environment
            }
            # subsample tagalign (non-mito) and cross-correlation analysis
            call xcor { input :
                ta = bam2ta_no_dedup.ta,
                subsample = xcor_subsample_reads,
                paired_end = paired_end_,
                mito_chr_name = mito_chr_name_,

                cpu = xcor_cpu,
                mem_factor = xcor_mem_factor,
                time_hr = xcor_time_hr,
                disk_factor = xcor_disk_factor,
                runtime_environment = runtime_environment_spp
            }
        }

        Boolean has_input_of_macs2_signal_track = has_output_of_bam2ta || defined(bam2ta.ta)
        if ( has_input_of_macs2_signal_track ) {
            # generate count signal track
            call macs2_signal_track { input :
                ta = ta_,
                gensz = gensz_,
                chrsz = chrsz_,
                pval_thresh = pval_thresh,
                smooth_win = smooth_win,

                mem_factor = macs2_signal_track_mem_factor,
                disk_factor = macs2_signal_track_disk_factor,
                time_hr = macs2_signal_track_time_hr,

                runtime_environment = runtime_environment_macs2
            }
        }

        Boolean has_input_of_call_peak = has_output_of_bam2ta || defined(bam2ta.ta)
        Boolean has_output_of_call_peak = i<length(peaks)
        if ( has_input_of_call_peak && !has_output_of_call_peak && !align_only ) {
            # call peaks on tagalign
            call call_peak { input :
                peak_caller = peak_caller_,
                peak_type = peak_type_,
                ta = ta_,
                gensz = gensz_,
                chrsz = chrsz_,
                cap_num_peak = cap_num_peak_,
                pval_thresh = pval_thresh,
                smooth_win = smooth_win,
                blacklist = blacklist_,
                regex_bfilt_peak_chr_name = regex_bfilt_peak_chr_name_,

                cpu = call_peak_cpu,
                mem_factor = call_peak_mem_factor,
                disk_factor = call_peak_disk_factor,
                time_hr = call_peak_time_hr,

                runtime_environment = if peak_caller == 'spp' then runtime_environment_spp
                    else if peak_caller == 'macs2' then runtime_environment_macs2
                    else runtime_environment
            }
        }
        File? peak_ = if has_output_of_call_peak then peaks[i] else call_peak.peak

        Boolean has_input_of_spr = has_output_of_bam2ta || defined(bam2ta.ta)
        if ( has_input_of_spr && !align_only && !true_rep_only ) {
            call spr { input :
                ta = ta_,
                paired_end = paired_end_,
                pseudoreplication_random_seed = pseudoreplication_random_seed,
                mem_factor = spr_mem_factor,
                disk_factor = spr_disk_factor,
                runtime_environment = runtime_environment
            }
        }

        Boolean has_input_of_call_peak_pr1 = defined(spr.ta_pr1)
        Boolean has_output_of_call_peak_pr1 = i<length(peaks_pr1)
        if ( has_input_of_call_peak_pr1 && !has_output_of_call_peak_pr1 &&
            !align_only && !true_rep_only ) {
            # call peaks on 1st pseudo replicated tagalign 
            call call_peak as call_peak_pr1 { input :
                peak_caller = peak_caller_,
                peak_type = peak_type_,
                ta = spr.ta_pr1,
                gensz = gensz_,
                chrsz = chrsz_,
                cap_num_peak = cap_num_peak_,
                pval_thresh = pval_thresh,
                smooth_win = smooth_win,
                blacklist = blacklist_,
                regex_bfilt_peak_chr_name = regex_bfilt_peak_chr_name_,

                cpu = call_peak_cpu,
                mem_factor = call_peak_mem_factor,
                disk_factor = call_peak_disk_factor,
                time_hr = call_peak_time_hr,

                runtime_environment = if peak_caller == 'spp' then runtime_environment_spp
                    else if peak_caller == 'macs2' then runtime_environment_macs2
                    else runtime_environment
            }
        }
        File? peak_pr1_ = if has_output_of_call_peak_pr1 then peaks_pr1[i]
            else call_peak_pr1.peak

        Boolean has_input_of_call_peak_pr2 = defined(spr.ta_pr2)
        Boolean has_output_of_call_peak_pr2 = i<length(peaks_pr2)
        if ( has_input_of_call_peak_pr2 && !has_output_of_call_peak_pr2 &&
            !align_only && !true_rep_only ) {
            # call peaks on 2nd pseudo replicated tagalign 
            call call_peak as call_peak_pr2 { input :
                peak_caller = peak_caller_,
                peak_type = peak_type_,
                ta = spr.ta_pr2,
                gensz = gensz_,
                chrsz = chrsz_,
                cap_num_peak = cap_num_peak_,
                pval_thresh = pval_thresh,
                smooth_win = smooth_win,
                blacklist = blacklist_,
                regex_bfilt_peak_chr_name = regex_bfilt_peak_chr_name_,

                cpu = call_peak_cpu,
                mem_factor = call_peak_mem_factor,
                disk_factor = call_peak_disk_factor,
                time_hr = call_peak_time_hr,

                runtime_environment = if peak_caller == 'spp' then runtime_environment_spp
                    else if peak_caller == 'macs2' then runtime_environment_macs2
                    else runtime_environment
            }
        }
        File? peak_pr2_ = if has_output_of_call_peak_pr2 then peaks_pr2[i]
            else call_peak_pr2.peak

        Boolean has_input_of_count_signal_track = has_output_of_bam2ta || defined(bam2ta.ta)
        if ( has_input_of_count_signal_track && enable_count_signal_track ) {
            # generate count signal track
            call count_signal_track { input :
                ta = ta_,
                chrsz = chrsz_,
                runtime_environment = runtime_environment
            }
        }
        # tasks factored out from ATAqC
        Boolean has_input_of_tss_enrich = defined(nodup_bam_) && defined(tss_) && (
            defined(align.read_len) || i<length(read_len) )
        if ( enable_tss_enrich && has_input_of_tss_enrich ) {
            call tss_enrich { input :
                read_len = if i<length(read_len) then read_len[i]
                    else align.read_len,
                nodup_bam = nodup_bam_,
                tss = tss_,
                chrsz = chrsz_,
                runtime_environment = runtime_environment_python2
            }
        }
        if ( enable_fraglen_stat && paired_end_ && defined(nodup_bam_) ) {
            call fraglen_stat_pe { input :
                nodup_bam = nodup_bam_,
                picard_java_heap = fraglen_stat_picard_java_heap,
                runtime_environment = runtime_environment
            }
        }
        if ( enable_preseq && defined(bam_) ) {
            call preseq { input :
                bam = bam_,
                paired_end = paired_end_,
                mem_factor = preseq_mem_factor,
                disk_factor = preseq_disk_factor,
                picard_java_heap = preseq_picard_java_heap,
                runtime_environment = runtime_environment
            }
        }
        if ( enable_gc_bias && defined(nodup_bam_) && defined(ref_fa_) ) {
            call gc_bias { input :
                nodup_bam = nodup_bam_,
                ref_fa = ref_fa_,
                picard_java_heap = gc_bias_picard_java_heap,
                runtime_environment = runtime_environment
            }
        }
        if ( enable_annot_enrich && defined(ta_) && defined(blacklist_) && defined(dnase_) && defined(prom_) && defined(enh_) ) {
            call annot_enrich { input :
                ta = ta_,
                blacklist = blacklist_,
                dnase = dnase_,
                prom = prom_,
                enh = enh_,
                runtime_environment = runtime_environment
            }
        }
        if ( enable_compare_to_roadmap && defined(macs2_signal_track.pval_bw) &&
             defined(reg2map_) && defined(roadmap_meta_) &&
             ( defined(reg2map_bed_) || defined(dnase_) ) ) {
            call compare_signal_to_roadmap { input :
                pval_bw = macs2_signal_track.pval_bw,
                dnase = dnase_,
                reg2map_bed = reg2map_bed_,
                reg2map = reg2map_,
                roadmap_meta = roadmap_meta_,
                runtime_environment = runtime_environment
            }
        }
    }

    # if there are TAs for ALL replicates then pool them
    Boolean has_all_inputs_of_pool_ta = length(select_all(ta_))==num_rep
    if ( has_all_inputs_of_pool_ta && num_rep>1 ) {
        # pool tagaligns from true replicates
        call pool_ta { input :
            tas = ta_,
            prefix = 'rep',

            runtime_environment = runtime_environment
        }
    }

    # if there are pr1 TAs for ALL replicates then pool them
    Boolean has_all_inputs_of_pool_ta_pr1 = length(select_all(spr.ta_pr1))==num_rep
    if ( has_all_inputs_of_pool_ta_pr1 && num_rep>1 && !align_only && !true_rep_only ) {
        # pool tagaligns from pseudo replicate 1
        call pool_ta as pool_ta_pr1 { input :
            tas = spr.ta_pr1,
            prefix = 'rep-pr1',
            runtime_environment = runtime_environment
        }
    }

    # if there are pr2 TAs for ALL replicates then pool them
    Boolean has_all_inputs_of_pool_ta_pr2 = length(select_all(spr.ta_pr2))==num_rep
    if ( has_all_inputs_of_pool_ta_pr1 && num_rep>1 && !align_only && !true_rep_only ) {
        # pool tagaligns from pseudo replicate 2
        call pool_ta as pool_ta_pr2 { input :
            tas = spr.ta_pr2,
            prefix = 'rep-pr2',
            runtime_environment = runtime_environment
        }
    }

    Boolean has_input_of_call_peak_pooled = defined(pool_ta.ta_pooled)
    Boolean has_output_of_call_peak_pooled = defined(peak_pooled)
    if ( has_input_of_call_peak_pooled && !has_output_of_call_peak_pooled &&
        !align_only && num_rep>1 ) {
        # call peaks on pooled replicate
        call call_peak as call_peak_pooled { input :
            peak_caller = peak_caller_,
            peak_type = peak_type_,
            ta = pool_ta.ta_pooled,
            gensz = gensz_,
            chrsz = chrsz_,
            cap_num_peak = cap_num_peak_,
            pval_thresh = pval_thresh,
            smooth_win = smooth_win,
            blacklist = blacklist_,
            regex_bfilt_peak_chr_name = regex_bfilt_peak_chr_name_,

            cpu = call_peak_cpu,
            mem_factor = call_peak_mem_factor,
            disk_factor = call_peak_disk_factor,
            time_hr = call_peak_time_hr,

            runtime_environment = if peak_caller == 'spp' then runtime_environment_spp
                else if peak_caller == 'macs2' then runtime_environment_macs2
                else runtime_environment
        }
    }
    File? peak_pooled_ = if has_output_of_call_peak_pooled then peak_pooled
        else call_peak_pooled.peak

    Boolean has_input_of_count_signal_track_pooled = defined(pool_ta.ta_pooled)
    if ( has_input_of_count_signal_track_pooled && enable_count_signal_track && num_rep>1 ) {
        call count_signal_track as count_signal_track_pooled { input :
            ta = pool_ta.ta_pooled,
            chrsz = chrsz_,
            runtime_environment = runtime_environment
        }
    }

    Boolean has_input_of_macs2_signal_track_pooled = defined(pool_ta.ta_pooled)
    if ( has_input_of_macs2_signal_track_pooled && num_rep>1 ) {
        call macs2_signal_track as macs2_signal_track_pooled { input :
            ta = pool_ta.ta_pooled,
            gensz = gensz_,
            chrsz = chrsz_,
            pval_thresh = pval_thresh,
            smooth_win = smooth_win,

            mem_factor = macs2_signal_track_mem_factor,
            disk_factor = macs2_signal_track_disk_factor,
            time_hr = macs2_signal_track_time_hr,
            runtime_environment = runtime_environment_macs2
        }
    }

    Boolean has_input_of_jsd = defined(blacklist_) &&
        length(select_all(nodup_bam_))==num_rep
    if ( has_input_of_jsd && num_rep > 0 && enable_jsd ) {
        # fingerprint and JS-distance plot
        call jsd { input :
            nodup_bams = nodup_bam_,
            blacklist = blacklist_,
            mapq_thresh = mapq_thresh_,

            cpu = jsd_cpu,
            mem_factor = jsd_mem_factor,
            time_hr = jsd_time_hr,
            disk_factor = jsd_disk_factor,
            runtime_environment = runtime_environment
        }
    }

    Boolean has_input_of_call_peak_ppr1 = defined(pool_ta_pr1.ta_pooled)
    Boolean has_output_of_call_peak_ppr1 = defined(peak_ppr1)
    if ( has_input_of_call_peak_ppr1 && !has_output_of_call_peak_ppr1 &&
        !align_only && !true_rep_only && num_rep>1 ) {
        # call peaks on 1st pooled pseudo replicates
        call call_peak as call_peak_ppr1 { input :
            peak_caller = peak_caller_,
            peak_type = peak_type_,
            ta = pool_ta_pr1.ta_pooled,
            gensz = gensz_,
            chrsz = chrsz_,
            cap_num_peak = cap_num_peak_,
            pval_thresh = pval_thresh,
            smooth_win = smooth_win,
            blacklist = blacklist_,
            regex_bfilt_peak_chr_name = regex_bfilt_peak_chr_name_,

            cpu = call_peak_cpu,
            mem_factor = call_peak_mem_factor,
            disk_factor = call_peak_disk_factor,
            time_hr = call_peak_time_hr,

            runtime_environment = if peak_caller == 'spp' then runtime_environment_spp
                else if peak_caller == 'macs2' then runtime_environment_macs2
                else runtime_environment
        }
    }
    File? peak_ppr1_ = if has_output_of_call_peak_ppr1 then peak_ppr1
        else call_peak_ppr1.peak

    Boolean has_input_of_call_peak_ppr2 = defined(pool_ta_pr2.ta_pooled)
    Boolean has_output_of_call_peak_ppr2 = defined(peak_ppr2)
    if ( has_input_of_call_peak_ppr2 && !has_output_of_call_peak_ppr2 &&
        !align_only && !true_rep_only && num_rep>1 ) {
        # call peaks on 2nd pooled pseudo replicates
        call call_peak as call_peak_ppr2 { input :
            peak_caller = peak_caller_,
            peak_type = peak_type_,
            ta = pool_ta_pr2.ta_pooled,
            gensz = gensz_,
            chrsz = chrsz_,
            cap_num_peak = cap_num_peak_,
            pval_thresh = pval_thresh,
            smooth_win = smooth_win,
            blacklist = blacklist_,
            regex_bfilt_peak_chr_name = regex_bfilt_peak_chr_name_,

            cpu = call_peak_cpu,
            mem_factor = call_peak_mem_factor,
            disk_factor = call_peak_disk_factor,
            time_hr = call_peak_time_hr,
            runtime_environment = if peak_caller == 'spp' then runtime_environment_spp
                else if peak_caller == 'macs2' then runtime_environment_macs2
                else runtime_environment
        }
    }
    File? peak_ppr2_ = if has_output_of_call_peak_ppr2 then peak_ppr2
        else call_peak_ppr2.peak

    # do IDR/overlap on all pairs of two replicates (i,j)
    #    where i and j are zero-based indices and 0 <= i < j < num_rep
    scatter( pair in cross(range(num_rep),range(num_rep)) ) {
        File? peak1_ = peak_[pair.left]
        File? peak2_ = peak_[pair.right]
        if ( !align_only && pair.left<pair.right ) {
            # pair.left = 0-based index of 1st replicate
            # pair.right = 0-based index of 2nd replicate
            # Naive overlap on every pair of true replicates
            call overlap { input :
                prefix = 'rep'+(pair.left+1)+'_vs_rep'+(pair.right+1),
                peak1 = peak1_,
                peak2 = peak2_,
                peak_pooled = peak_pooled_,
                peak_type = peak_type_,
                blacklist = blacklist_,
                chrsz = chrsz_,
                regex_bfilt_peak_chr_name = regex_bfilt_peak_chr_name_,
                ta = pool_ta.ta_pooled,
                runtime_environment = runtime_environment
            }
        }
        if ( enable_idr && !align_only && pair.left<pair.right ) {
            # pair.left = 0-based index of 1st replicate
            # pair.right = 0-based index of 2nd replicate
            # IDR on every pair of true replicates
            call idr { input :
                prefix = 'rep'+(pair.left+1)+'_vs_rep'+(pair.right+1),
                peak1 = peak1_,
                peak2 = peak2_,
                peak_pooled = peak_pooled_,
                idr_thresh = idr_thresh,
                peak_type = peak_type_,
                rank = idr_rank_,
                blacklist = blacklist_,
                chrsz = chrsz_,
                regex_bfilt_peak_chr_name = regex_bfilt_peak_chr_name_,
                ta = pool_ta.ta_pooled,
                runtime_environment = runtime_environment
            }
        }
    }

    # overlap on pseudo-replicates (pr1, pr2) for each true replicate
    if ( !align_only && !true_rep_only ) {
        scatter( i in range(num_rep) ) {
            call overlap as overlap_pr { input :
                prefix = 'rep'+(i+1)+'-pr1_vs_rep'+(i+1)+'-pr2',
                peak1 = peak_pr1_[i],
                peak2 = peak_pr2_[i],
                peak_pooled = peak_[i],
                peak_type = peak_type_,
                blacklist = blacklist_,
                chrsz = chrsz_,
                regex_bfilt_peak_chr_name = regex_bfilt_peak_chr_name_,
                ta = ta_[i],
                runtime_environment = runtime_environment
            }
        }
    }

    if ( !align_only && !true_rep_only && enable_idr ) {
        scatter( i in range(num_rep) ) {
            # IDR on pseduo replicates
            call idr as idr_pr { input :
                prefix = 'rep'+(i+1)+'-pr1_vs_rep'+(i+1)+'-pr2',
                peak1 = peak_pr1_[i],
                peak2 = peak_pr2_[i],
                peak_pooled = peak_[i],
                idr_thresh = idr_thresh,
                peak_type = peak_type_,
                rank = idr_rank_,
                blacklist = blacklist_,
                chrsz = chrsz_,
                regex_bfilt_peak_chr_name = regex_bfilt_peak_chr_name_,
                ta = ta_[i],
                runtime_environment = runtime_environment
            }
        }
    }

    if ( !align_only && !true_rep_only && num_rep>1 ) {
        # Naive overlap on pooled pseudo replicates
        call overlap as overlap_ppr { input :
            prefix = 'pooled-pr1_vs_pooled-pr2',
            peak1 = peak_ppr1_,
            peak2 = peak_ppr2_,
            peak_pooled = peak_pooled_,
            peak_type = peak_type_,
            blacklist = blacklist_,
            chrsz = chrsz_,
            regex_bfilt_peak_chr_name = regex_bfilt_peak_chr_name_,
            ta = pool_ta.ta_pooled,
            runtime_environment = runtime_environment
        }
    }

    if ( !align_only && !true_rep_only && num_rep>1 ) {
        # IDR on pooled pseduo replicates
        call idr as idr_ppr { input :
            prefix = 'pooled-pr1_vs_pooled-pr2',
            peak1 = peak_ppr1_,
            peak2 = peak_ppr2_,
            peak_pooled = peak_pooled_,
            idr_thresh = idr_thresh,
            peak_type = peak_type_,
            rank = idr_rank_,
            blacklist = blacklist_,
            chrsz = chrsz_,
            regex_bfilt_peak_chr_name = regex_bfilt_peak_chr_name_,
            ta = pool_ta.ta_pooled,
            runtime_environment = runtime_environment
        }
    }

    # reproducibility QC for overlap/IDR peaks
    if ( !align_only && !true_rep_only && num_rep > 0 ) {
        # reproducibility QC for overlapping peaks
        call reproducibility as reproducibility_overlap { input :
            prefix = 'overlap',
            peaks = select_all(overlap.bfilt_overlap_peak),
            peaks_pr = if defined(overlap_pr.bfilt_overlap_peak) then select_first([overlap_pr.bfilt_overlap_peak]) else [],
            peak_ppr = overlap_ppr.bfilt_overlap_peak,
            peak_type = peak_type_,
            chrsz = chrsz_,
            runtime_environment = runtime_environment
        }
    }

    if ( !align_only && !true_rep_only && num_rep > 0 && enable_idr ) {
        # reproducibility QC for IDR peaks
        call reproducibility as reproducibility_idr { input :
            prefix = 'idr',
            peaks = select_all(idr.bfilt_idr_peak),
            peaks_pr = if defined(idr_pr.bfilt_idr_peak) then select_first([idr_pr.bfilt_idr_peak]) else [],
            peak_ppr = idr_ppr.bfilt_idr_peak,
            peak_type = peak_type_,
            chrsz = chrsz_,
            runtime_environment = runtime_environment
        }
    }

    # Generate final QC report and JSON
    call qc_report { input :
        pipeline_ver = pipeline_ver,
        title = title,
        description = description,
        genome = genome_name_,
        multimapping = multimapping,
        paired_ends = paired_end_,
        pipeline_type = pipeline_type,
        aligner = aligner_,
        no_dup_removal = no_dup_removal,
        peak_caller = peak_caller_,
        cap_num_peak = cap_num_peak_,
        idr_thresh = idr_thresh,
        pval_thresh = pval_thresh,
        xcor_subsample_reads = xcor_subsample_reads,

        samstat_qcs = select_all(align.samstat_qc),
        nodup_samstat_qcs = select_all(filter.samstat_qc),

        frac_mito_qcs = select_all(frac_mito.frac_mito_qc),
        dup_qcs = select_all(filter.dup_qc),
        lib_complexity_qcs = select_all(filter.lib_complexity_qc),
        xcor_plots = select_all(xcor.plot_png),
        xcor_scores = select_all(xcor.score),

        jsd_plot = jsd.plot,
        jsd_qcs = if defined(jsd.jsd_qcs) then select_first([jsd.jsd_qcs]) else [],

        frip_qcs = select_all(call_peak.frip_qc),
        frip_qcs_pr1 = select_all(call_peak_pr1.frip_qc),
        frip_qcs_pr2 = select_all(call_peak_pr2.frip_qc),

        frip_qc_pooled = call_peak_pooled.frip_qc,
        frip_qc_ppr1 = call_peak_ppr1.frip_qc,
        frip_qc_ppr2 = call_peak_ppr2.frip_qc,

        idr_plots = select_all(idr.idr_plot),
        idr_plots_pr = if defined(idr_pr.idr_plot) then select_first([idr_pr.idr_plot]) else [],
        idr_plot_ppr = idr_ppr.idr_plot,
        frip_idr_qcs = select_all(idr.frip_qc),
        frip_idr_qcs_pr = if defined(idr_pr.frip_qc) then select_first([idr_pr.frip_qc]) else [],
        frip_idr_qc_ppr = idr_ppr.frip_qc,
        frip_overlap_qcs = select_all(overlap.frip_qc),
        frip_overlap_qcs_pr = if defined(overlap_pr.frip_qc) then select_first([overlap_pr.frip_qc]) else [],
        frip_overlap_qc_ppr = overlap_ppr.frip_qc,
        idr_reproducibility_qc = reproducibility_idr.reproducibility_qc,
        overlap_reproducibility_qc = reproducibility_overlap.reproducibility_qc,

        annot_enrich_qcs = select_all(annot_enrich.annot_enrich_qc),
        tss_enrich_qcs = select_all(tss_enrich.tss_enrich_qc),
        tss_large_plots = select_all(tss_enrich.tss_large_plot),
        roadmap_compare_plots = select_all(compare_signal_to_roadmap.roadmap_compare_plot),
        fraglen_dist_plots = select_all(fraglen_stat_pe.fraglen_dist_plot),
        fraglen_nucleosomal_qcs = select_all(fraglen_stat_pe.nucleosomal_qc),
        gc_plots = select_all(gc_bias.gc_plot),
        preseq_plots = select_all(preseq.preseq_plot),
        picard_est_lib_size_qcs = select_all(preseq.picard_est_lib_size_qc),

        peak_region_size_qcs = select_all(call_peak.peak_region_size_qc),
        peak_region_size_plots = select_all(call_peak.peak_region_size_plot),
        num_peak_qcs = select_all(call_peak.num_peak_qc),

        idr_opt_peak_region_size_qc = reproducibility_idr.peak_region_size_qc,
        idr_opt_peak_region_size_plot = reproducibility_overlap.peak_region_size_plot,
        idr_opt_num_peak_qc = reproducibility_idr.num_peak_qc,

        overlap_opt_peak_region_size_qc = reproducibility_overlap.peak_region_size_qc,
        overlap_opt_peak_region_size_plot = reproducibility_overlap.peak_region_size_plot,
        overlap_opt_num_peak_qc = reproducibility_overlap.num_peak_qc,

        runtime_environment = runtime_environment
    }

    output {
        File report = qc_report.report
        File qc_json = qc_report.qc_json
        Boolean qc_json_ref_match = qc_report.qc_json_ref_match
    }
}

task align {
    input {        
        # for task trim_adapter
        Array[File] fastqs_R1         # [merge_id]
        Array[File] fastqs_R2

        String? adapter     # adapter for all fastqs,
                            #    this will override individual adapters in adapters_R1/R2
        Array[String] adapters_R1
        Array[String] adapters_R2
        Boolean paired_end
        Boolean auto_detect_adapter
        String cutadapt_param

        # for task align
        String aligner
        String mito_chr_name
        File chrsz            # 2-col chromosome sizes file
        File idx_tar        # reference index tar or tar.gz
        Int multimapping

        # resource
        Int cpu
        Float mem_factor
        Int time_hr
        Float disk_factor

        RuntimeEnvironment runtime_environment
    }
    Float input_file_size_gb = size(fastqs_R1, "G") + size(fastqs_R2, "G")
    Float mem_gb = 5.0 + size(idx_tar, "G") + mem_factor * input_file_size_gb
    Float samtools_mem_gb = 0.8 * mem_gb
    Int disk_gb = round(40.0 + disk_factor * input_file_size_gb)

    # tmp vars for task trim_adapter
    Array[Array[File]] tmp_fastqs = if paired_end then transpose([fastqs_R1, fastqs_R2])
                else transpose([fastqs_R1])
    Array[Array[String]] tmp_adapters = if paired_end then transpose([adapters_R1, adapters_R2])
                else transpose([adapters_R1])
    command {
        set -e

        # check if pipeline dependencies can be found
        if [[ -z "$(which encode_task_trim_adapter.py 2> /dev/null || true)" ]]
        then
          echo -e "\n* Error: pipeline environment (docker, singularity or conda) not found." 1>&2
          exit 3
        fi

        # trim adapter
        python3 $(which encode_task_trim_adapter.py) \
            ${write_tsv(tmp_fastqs)} \
            ${'--adapter ' + adapter} \
            --adapters ${write_tsv(tmp_adapters)} \
            ${if paired_end then '--paired-end' else ''} \
            ${if auto_detect_adapter then '--auto-detect-adapter' else ''} \
            --cutadapt-param ' ${cutadapt_param}' \
            ${'--nth ' + cpu}

        # align on trimmed/merged fastqs
        if [ '${aligner}' == 'bowtie2' ]; then
            python3 $(which encode_task_bowtie2.py) \
                ${idx_tar} \
                R1/*.fastq.gz \
                ${if paired_end then 'R2/*.fastq.gz' else ''} \
                ${if paired_end then '--paired-end' else ''} \
                ${'--multimapping ' + multimapping} \
                ${'--mem-gb ' + samtools_mem_gb} \
                ${'--nth ' + cpu}
        fi

        python3 $(which encode_task_post_align.py) \
            R1/*.fastq.gz $(ls *.bam) \
            ${'--mito-chr-name ' + mito_chr_name} \
            ${'--chrsz ' + chrsz} \
            ${'--mem-gb ' + samtools_mem_gb} \
            ${'--nth ' + cpu}
        rm -rf R1 R2
    }
    output {
        File bam = glob('*.bam')[0]
        File bai = glob('*.bai')[0]
        File samstat_qc = glob('*.samstats.qc')[0]
        File non_mito_samstat_qc = glob('non_mito/*.samstats.qc')[0]
        File read_len_log = glob('*.read_length.txt')[0]
        Int read_len = read_int(read_len_log)
    }
    runtime {
        cpu : cpu
        memory : '${mem_gb} GB'
        time : time_hr
        disks : 'local-disk ${disk_gb} SSD'
        preemptible: 0

        docker : runtime_environment.docker
        singularity : runtime_environment.singularity
        conda : runtime_environment.conda
    }
}

task frac_mito {
    input {
        File? non_mito_samstat
        File? mito_samstat

        RuntimeEnvironment runtime_environment
    }

    command {
        set -e
        python3 $(which encode_task_frac_mito.py) \
            ${non_mito_samstat} ${mito_samstat}
    }
    output {
        File frac_mito_qc = glob('*.frac_mito.qc')[0]
    }
    runtime {
        cpu : 1
        memory : '4 GB'
        time : 4
        disks : 'local-disk 10 SSD'

        docker : runtime_environment.docker
        singularity : runtime_environment.singularity
        conda : runtime_environment.conda
    }
}

task filter {
    input {
        File? bam
        Boolean paired_end
        Int multimapping
        String dup_marker             # picard.jar MarkDuplicates (picard) or 
                                    # sambamba markdup (sambamba)
        Int mapq_thresh                # threshold for low MAPQ reads removal
        Array[String] filter_chrs     # chrs to be removed from final (nodup/filt) BAM
        File chrsz                    # 2-col chromosome sizes file
        Boolean no_dup_removal         # no dupe reads removal when filtering BAM
        String mito_chr_name

        Int cpu
        Float mem_factor
        String? picard_java_heap
        Int time_hr
        Float disk_factor
        # runtime environment
        RuntimeEnvironment runtime_environment
    }
    Float input_file_size_gb = size(bam, "G")
    Float picard_java_heap_factor = 0.9
    Float mem_gb = 6.0 + mem_factor * input_file_size_gb
    Float samtools_mem_gb = 0.8 * mem_gb
    Int disk_gb = round(20.0 + disk_factor * input_file_size_gb)

    command {
        set -e
        python3 $(which encode_task_filter.py) \
            ${bam} \
            ${if paired_end then '--paired-end' else ''} \
            ${'--multimapping ' + multimapping} \
            ${'--dup-marker ' + dup_marker} \
            ${'--mapq-thresh ' + mapq_thresh} \
            --filter-chrs ${sep=' ' filter_chrs} \
            ${'--chrsz ' + chrsz} \
            ${if no_dup_removal then '--no-dup-removal' else ''} \
            ${'--mito-chr-name ' + mito_chr_name} \
            ${'--mem-gb ' + samtools_mem_gb} \
            ${'--nth ' + cpu} \
            ${'--picard-java-heap ' + if defined(picard_java_heap) then picard_java_heap else (round(mem_gb * picard_java_heap_factor) + 'G')}
    }
    output {
        File nodup_bam = glob('*.bam')[0]
        File nodup_bai = glob('*.bai')[0]
        File samstat_qc = glob('*.samstats.qc')[0]
        File dup_qc = glob('*.dup.qc')[0]
        File lib_complexity_qc = glob('*.lib_complexity.qc')[0]
    }
    runtime {
        cpu : cpu
        memory : '${mem_gb} GB'
        time : time_hr
        disks : 'local-disk ${disk_gb} SSD'

        docker : runtime_environment.docker
        singularity : runtime_environment.singularity
        conda : runtime_environment.conda
    }
}

task bam2ta {
    input {
        File? bam
        Boolean paired_end
        Boolean disable_tn5_shift     # no tn5 shifting (it's for dnase-seq)
        String mito_chr_name         # mito chromosome name
        Int subsample                 # number of reads to subsample TAGALIGN
                                    # this affects all downstream analysis
        Int cpu
        Float mem_factor
        Int time_hr
        Float disk_factor
        # runtime environment
        RuntimeEnvironment runtime_environment
    }
    Float input_file_size_gb = size(bam, "G")
    Float mem_gb = 4.0 + mem_factor * input_file_size_gb
    Float samtools_mem_gb = 0.8 * mem_gb
    Int disk_gb = round(20.0 + disk_factor * input_file_size_gb)

    command {
        set -e
        python3 $(which encode_task_bam2ta.py) \
            ${bam} \
            ${if paired_end then '--paired-end' else ''} \
            ${if disable_tn5_shift then '--disable-tn5-shift' else ''} \
            ${'--mito-chr-name ' + mito_chr_name} \
            ${'--subsample ' + subsample} \
            ${'--mem-gb ' + samtools_mem_gb} \
            ${'--nth ' + cpu}
    }
    output {
        File ta = glob('*.tagAlign.gz')[0]
    }
    runtime {
        cpu : cpu
        memory : '${mem_gb} GB'
        time : time_hr
        disks : 'local-disk ${disk_gb} SSD'

        docker : runtime_environment.docker
        singularity : runtime_environment.singularity
        conda : runtime_environment.conda
    }
}

task spr {
    input {
        File? ta
        Boolean paired_end
        Int pseudoreplication_random_seed

        Float mem_factor
        Float disk_factor

        # runtime environment
        RuntimeEnvironment runtime_environment
    }
    Float input_file_size_gb = size(ta, "G")
    Float mem_gb = 4.0 + mem_factor * input_file_size_gb
    Int disk_gb = round(20.0 + disk_factor * input_file_size_gb)

    command {
        set -e
        python3 $(which encode_task_spr.py) \
            ${ta} \
            ${'--pseudoreplication-random-seed ' + pseudoreplication_random_seed} \
            ${if paired_end then '--paired-end' else ''}
    }
    output {
        File ta_pr1 = glob('*.pr1.tagAlign.gz')[0]
        File ta_pr2 = glob('*.pr2.tagAlign.gz')[0]
    }
    runtime {
        cpu : 1
        memory : '${mem_gb} GB'
        time : 4
        disks : 'local-disk ${disk_gb} SSD'

        docker : runtime_environment.docker
        singularity : runtime_environment.singularity
        conda : runtime_environment.conda
    }
}

task pool_ta {
    input {
        Array[File?] tas     # TAG-ALIGNs to be merged
        Int? col             # number of columns in pooled TA
        String? prefix         # basename prefix

        # runtime environment
        RuntimeEnvironment runtime_environment
    }
    command {
        set -e
        python3 $(which encode_task_pool_ta.py) \
            ${sep=' ' select_all(tas)} \
            ${'--prefix ' + prefix} \
            ${'--col ' + col}
    }
    output {
        File ta_pooled = glob('*.tagAlign.gz')[0]
    }
    runtime {
        cpu : 1
        memory : '8 GB'
        time : 4
        disks : 'local-disk 100 SSD'

        docker : runtime_environment.docker
        singularity : runtime_environment.singularity
        conda : runtime_environment.conda
    }
}

task xcor {
    input {
        File? ta
        Boolean paired_end
        String mito_chr_name
        Int subsample  # number of reads to subsample TAGALIGN
                    # this will be used for xcor only
                    # will not affect any downstream analysis
        Int cpu
        Float mem_factor
        Int time_hr
        Float disk_factor

        # runtime environment
        RuntimeEnvironment runtime_environment
    }
    Float input_file_size_gb = size(ta, "G")
    Float mem_gb = 8.0 + mem_factor * input_file_size_gb
    Int disk_gb = round(20.0 + disk_factor * input_file_size_gb)

    command {
        set -e
        python3 $(which encode_task_xcor.py) \
            ${ta} \
            ${if paired_end then '--paired-end' else ''} \
            ${'--mito-chr-name ' + mito_chr_name} \
            ${'--subsample ' + subsample} \
            --speak=0 \
            ${'--nth ' + cpu}
    }
    output {
        File plot_pdf = glob('*.cc.plot.pdf')[0]
        File plot_png = glob('*.cc.plot.png')[0]
        File score = glob('*.cc.qc')[0]
    }
    runtime {
        cpu : cpu
        memory : '${mem_gb} GB'
        time : time_hr
        disks : 'local-disk ${disk_gb} SSD'

        docker : runtime_environment.docker
        singularity : runtime_environment.singularity
        conda : runtime_environment.conda
    }
}

task jsd {
    input {
        Array[File?] nodup_bams
        File? blacklist
        Int mapq_thresh

        Int cpu
        Float mem_factor
        Int time_hr
        Float disk_factor

        # runtime environment
        RuntimeEnvironment runtime_environment
    }
    Float input_file_size_gb = size(nodup_bams, "G")
    Float mem_gb = 5.0 + mem_factor * input_file_size_gb
    Int disk_gb = round(20.0 + disk_factor * input_file_size_gb)

    command {
        set -e
        python3 $(which encode_task_jsd.py) \
            ${sep=' ' select_all(nodup_bams)} \
            ${'--mapq-thresh '+ mapq_thresh} \
            ${'--blacklist '+ blacklist} \
            ${'--nth ' + cpu}
    }
    output {
        File plot = glob('*.png')[0]
        Array[File] jsd_qcs = glob('*.jsd.qc')
    }
    runtime {
        cpu : cpu
        memory : '${mem_gb} GB'
        time : time_hr
        disks : 'local-disk ${disk_gb} SSD'

        docker : runtime_environment.docker
        singularity : runtime_environment.singularity
        conda : runtime_environment.conda
    }
}

task count_signal_track {
    input {
        File? ta             # tag-align
        File chrsz            # 2-col chromosome sizes file

        # runtime environment
        RuntimeEnvironment runtime_environment
    }
    Float mem_gb = 8.0
    command {
        set -e
        python3 $(which encode_task_count_signal_track.py) \
            ${ta} \
            ${'--chrsz ' + chrsz} \
            ${'--mem-gb ' + mem_gb}
    }
    output {
        File pos_bw = glob('*.positive.bigwig')[0]
        File neg_bw = glob('*.negative.bigwig')[0]
    }
    runtime {
        cpu : 1
        memory : '${mem_gb} GB'
        time : 4
        disks : 'local-disk 50 SSD'

        docker : runtime_environment.docker
        singularity : runtime_environment.singularity
        conda : runtime_environment.conda
    }
}

task call_peak {
    input {
        String peak_caller
        String peak_type

        File? ta
        String gensz        # Genome size (sum of entries in 2nd column of 
                            # chr. sizes file, or hs for human, ms for mouse)
        File chrsz            # 2-col chromosome sizes file
        Int cap_num_peak    # cap number of raw peaks called from MACS2
        Float pval_thresh      # p.value threshold
        Int smooth_win         # size of smoothing window
        File? blacklist     # blacklist BED to filter raw peaks
        String? regex_bfilt_peak_chr_name

        Int cpu
        Float mem_factor
        Int time_hr
        Float disk_factor

        # runtime environment
        RuntimeEnvironment runtime_environment
    }
    Float input_file_size_gb = size(ta, "G")
    Float mem_gb = 4.0 + mem_factor * input_file_size_gb
    Int disk_gb = round(20.0 + disk_factor * input_file_size_gb)

    command {
        set -e

        if [ '${peak_caller}' == 'macs2' ]; then
            python3 $(which encode_task_macs2_atac.py) \
                ${ta} \
                ${'--gensz ' + gensz} \
                ${'--chrsz ' + chrsz} \
                ${'--cap-num-peak ' + cap_num_peak} \
                ${'--pval-thresh '+ pval_thresh} \
                ${'--smooth-win '+ smooth_win} \
                ${'--mem-gb ' + mem_gb}
        fi

        python3 $(which encode_task_post_call_peak_atac.py) \
            $(ls *Peak.gz) \
            ${'--ta ' + ta} \
            ${'--regex-bfilt-peak-chr-name \'' + regex_bfilt_peak_chr_name + '\''} \
            ${'--chrsz ' + chrsz} \
            ${'--peak-type ' + peak_type} \
            ${'--blacklist ' + blacklist}
    }
    output {
        File peak = glob('*[!.][!b][!f][!i][!l][!t].'+peak_type+'.gz')[0]
        # generated by post_call_peak py
        File bfilt_peak = glob('*.bfilt.'+peak_type+'.gz')[0]
        File bfilt_peak_bb = glob('*.bfilt.'+peak_type+'.bb')[0]
        File bfilt_peak_starch = glob('*.bfilt.'+peak_type+'.starch')[0]
        File bfilt_peak_hammock = glob('*.bfilt.'+peak_type+'.hammock.gz*')[0]
        File bfilt_peak_hammock_tbi = glob('*.bfilt.'+peak_type+'.hammock.gz*')[1]
        File frip_qc = glob('*.frip.qc')[0]
        File peak_region_size_qc = glob('*.peak_region_size.qc')[0]
        File peak_region_size_plot = glob('*.peak_region_size.png')[0]
        File num_peak_qc = glob('*.num_peak.qc')[0]
    }
    runtime {
        cpu : if peak_caller == 'macs2' then 1 else cpu
        memory : '${mem_gb} GB'
        time : time_hr
        disks : 'local-disk ${disk_gb} SSD'
        preemptible: 0

        docker : runtime_environment.docker
        singularity : runtime_environment.singularity
        conda : runtime_environment.conda
    }
}

task macs2_signal_track {
    input {
        File? ta
        String gensz        # Genome size (sum of entries in 2nd column of 
                            # chr. sizes file, or hs for human, ms for mouse)
        File chrsz            # 2-col chromosome sizes file
        Float pval_thresh      # p.value threshold
        Int smooth_win         # size of smoothing window

        Float mem_factor
        Int time_hr
        Float disk_factor

        # runtime environment
        RuntimeEnvironment runtime_environment
    }
    Float input_file_size_gb = size(ta, "G")
    Float mem_gb = 4.0 + mem_factor * input_file_size_gb
    Int disk_gb = round(20.0 + disk_factor * input_file_size_gb)

    command {
        set -e
        python3 $(which encode_task_macs2_signal_track_atac.py) \
            ${ta} \
            ${'--gensz '+ gensz} \
            ${'--chrsz ' + chrsz} \
            ${'--pval-thresh '+ pval_thresh} \
            ${'--smooth-win '+ smooth_win} \
            ${'--mem-gb ' + mem_gb}
    }
    output {
        File pval_bw = glob('*.pval.signal.bigwig')[0]
        File fc_bw = glob('*.fc.signal.bigwig')[0]
    }
    runtime {
        cpu : 1
        memory : '${mem_gb} GB'
        time : time_hr
        disks : 'local-disk ${disk_gb} SSD'
        preemptible: 0

        docker : runtime_environment.docker
        singularity : runtime_environment.singularity
        conda : runtime_environment.conda
    }
}

task idr {
    input {
        String prefix         # prefix for IDR output file
        File? peak1
        File? peak2
        File? peak_pooled
        Float idr_thresh
        File? blacklist     # blacklist BED to filter raw peaks
        String regex_bfilt_peak_chr_name
        # parameters to compute FRiP
        File? ta            # to calculate FRiP
        File chrsz            # 2-col chromosome sizes file
        String peak_type
        String rank

        # runtime environment
        RuntimeEnvironment runtime_environment
    }
    command {
        set -e
        touch null
        python3 $(which encode_task_idr.py) \
            ${peak1} ${peak2} ${peak_pooled} \
            ${'--prefix ' + prefix} \
            ${'--idr-thresh ' + idr_thresh} \
            ${'--peak-type ' + peak_type} \
            --idr-rank ${rank} \
            ${'--chrsz ' + chrsz} \
            ${'--blacklist '+ blacklist} \
            ${'--regex-bfilt-peak-chr-name \'' + regex_bfilt_peak_chr_name + '\''} \
            ${'--ta ' + ta}
    }
    output {
        File idr_peak = glob('*[!.][!b][!f][!i][!l][!t].'+peak_type+'.gz')[0]
        File bfilt_idr_peak = glob('*.bfilt.'+peak_type+'.gz')[0]
        File bfilt_idr_peak_bb = glob('*.bfilt.'+peak_type+'.bb')[0]
        File bfilt_idr_peak_starch = glob('*.bfilt.'+peak_type+'.starch')[0]
        File bfilt_idr_peak_hammock = glob('*.bfilt.'+peak_type+'.hammock.gz*')[0]
        File bfilt_idr_peak_hammock_tbi = glob('*.bfilt.'+peak_type+'.hammock.gz*')[1]
        File idr_plot = glob('*.txt.png')[0]
        File idr_unthresholded_peak = glob('*.txt.gz')[0]
        File idr_log = glob('*.idr*.log')[0]
        File frip_qc = if defined(ta) then glob('*.frip.qc')[0] else glob('null')[0]
    }
    runtime {
        cpu : 1
        memory : '4 GB'
        time : 4
        disks : 'local-disk 50 SSD'

        docker : runtime_environment.docker
        singularity : runtime_environment.singularity
        conda : runtime_environment.conda
    }
}

task overlap {
    input {
        String prefix         # prefix for IDR output file
        File? peak1
        File? peak2
        File? peak_pooled
        File? blacklist     # blacklist BED to filter raw peaks
        String regex_bfilt_peak_chr_name
        File? ta        # to calculate FRiP
        File chrsz            # 2-col chromosome sizes file
        String peak_type

        # runtime environment
        RuntimeEnvironment runtime_environment
    }
    command {
        set -e
        touch null 
        python3 $(which encode_task_overlap.py) \
            ${peak1} ${peak2} ${peak_pooled} \
            ${'--prefix ' + prefix} \
            ${'--peak-type ' + peak_type} \
            ${'--chrsz ' + chrsz} \
            ${'--blacklist '+ blacklist} \
            --nonamecheck \
            ${'--regex-bfilt-peak-chr-name \'' + regex_bfilt_peak_chr_name + '\''} \
            ${'--ta ' + ta}
    }
    output {
        File overlap_peak = glob('*[!.][!b][!f][!i][!l][!t].'+peak_type+'.gz')[0]
        File bfilt_overlap_peak = glob('*.bfilt.'+peak_type+'.gz')[0]
        File bfilt_overlap_peak_bb = glob('*.bfilt.'+peak_type+'.bb')[0]
        File bfilt_overlap_peak_starch = glob('*.bfilt.'+peak_type+'.starch')[0]
        File bfilt_overlap_peak_hammock = glob('*.bfilt.'+peak_type+'.hammock.gz*')[0]
        File bfilt_overlap_peak_hammock_tbi = glob('*.bfilt.'+peak_type+'.hammock.gz*')[1]
        File frip_qc = if defined(ta) then glob('*.frip.qc')[0] else glob('null')[0]
    }
    runtime {
        cpu : 1
        memory : '4 GB'
        time : 4
        disks : 'local-disk 50 SSD'

        docker : runtime_environment.docker
        singularity : runtime_environment.singularity
        conda : runtime_environment.conda
    }
}

task reproducibility {
    input {
        String prefix
        Array[File] peaks # peak files from pair of true replicates
                            # in a sorted order. for example of 4 replicates,
                            # 1,2 1,3 1,4 2,3 2,4 3,4.
                            # x,y means peak file from rep-x vs rep-y
        Array[File] peaks_pr    # peak files from pseudo replicates
        File? peak_ppr            # Peak file from pooled pseudo replicate.
        String peak_type
        File chrsz            # 2-col chromosome sizes file

        # runtime environment
        RuntimeEnvironment runtime_environment
    }
    command {
        set -e
        python3 $(which encode_task_reproducibility.py) \
            ${sep=' ' peaks} \
            --peaks-pr ${sep=' ' peaks_pr} \
            ${'--peak-ppr '+ peak_ppr} \
            --prefix ${prefix} \
            ${'--peak-type ' + peak_type} \
            ${'--chrsz ' + chrsz}
    }
    output {
        File optimal_peak = glob('*optimal_peak.*.gz')[0]
        File optimal_peak_bb = glob('*optimal_peak.*.bb')[0]
        File optimal_peak_starch = glob('*optimal_peak.*.starch')[0]
        File optimal_peak_hammock = glob('*optimal_peak.*.hammock.gz*')[0]
        File optimal_peak_hammock_tbi = glob('*optimal_peak.*.hammock.gz*')[1]
        File conservative_peak = glob('*conservative_peak.*.gz')[0]
        File conservative_peak_bb = glob('*conservative_peak.*.bb')[0]
        File conservative_peak_starch = glob('*conservative_peak.*.starch')[0]
        File conservative_peak_hammock = glob('*conservative_peak.*.hammock.gz*')[0]
        File conservative_peak_hammock_tbi = glob('*conservative_peak.*.hammock.gz*')[1]
        File reproducibility_qc = glob('*reproducibility.qc')[0]
        # QC metrics for optimal peak
        File peak_region_size_qc = glob('*.peak_region_size.qc')[0]
        File peak_region_size_plot = glob('*.peak_region_size.png')[0]
        File num_peak_qc = glob('*.num_peak.qc')[0]
    }
    runtime {
        cpu : 1
        memory : '4 GB'
        time : 4
        disks : 'local-disk 50 SSD'

        docker : runtime_environment.docker
        singularity : runtime_environment.singularity
        conda : runtime_environment.conda
    }
}

task preseq {
    input {
        File? bam
        Boolean paired_end

        Float mem_factor
        Float disk_factor
        String? picard_java_heap
        File? null

        # runtime environment
        RuntimeEnvironment runtime_environment
    }
    Float input_file_size_gb = size(bam, "G")
    Float mem_gb = 4.0 + mem_factor * input_file_size_gb
    Float samtools_mem_gb = 0.8 * mem_gb
    Int disk_gb = round(20.0 + disk_factor * input_file_size_gb)
    Float picard_java_heap_factor = 0.9

    command {
        set -e
        python3 $(which encode_task_preseq.py) \
            ${if paired_end then '--paired-end' else ''} \
            ${'--bam ' + bam} \
            ${'--mem-gb ' + samtools_mem_gb} \
            ${'--picard-java-heap ' + if defined(picard_java_heap) then picard_java_heap else (round(mem_gb * picard_java_heap_factor) + 'G')}
        ${if !paired_end then 'touch null.picard_est_lib_size' else ''}
    }
    output {
        File? picard_est_lib_size_qc = if paired_end then 
            glob('*.picard_est_lib_size.qc')[0] else null
        File preseq_plot = glob('*.preseq.png')[0]
        File preseq_log = glob('*.preseq.log')[0]
    }
    runtime {
        cpu : 1
        memory : '${mem_gb} GB'
        time : 4
        disks : 'local-disk ${disk_gb} SSD'

        docker : runtime_environment.docker
        singularity : runtime_environment.singularity
        conda : runtime_environment.conda
    }
}

task annot_enrich {
    input {
        # Fraction of Reads In Annotated Regions
        File? ta
        File? blacklist
        File? dnase
        File? prom
        File? enh

        # runtime environment
        RuntimeEnvironment runtime_environment
    }
    command {
        set -e
        python3 $(which encode_task_annot_enrich.py) \
            ${'--ta ' + ta} \
            ${'--blacklist ' + blacklist} \
            ${'--dnase ' + dnase} \
            ${'--prom ' + prom} \
            ${'--enh ' + enh}
    }
    output {
        File annot_enrich_qc = glob('*.annot_enrich.qc')[0]
    }
    runtime {
        cpu : 1
        memory : '8 GB'
        time : 4
        disks : 'local-disk 50 SSD'

        docker : runtime_environment.docker
        singularity : runtime_environment.singularity
        conda : runtime_environment.conda
    }
}

task tss_enrich {
    # based on metaseq, which is still in python2
    # python2 environment is required for this task
    input {
        Int? read_len
        File? nodup_bam
        File? tss
        File chrsz

        # runtime environment
        RuntimeEnvironment runtime_environment
    }
    command {
        set -e
        python2 $(which encode_task_tss_enrich.py) \
            ${'--read-len ' + read_len} \
            ${'--nodup-bam ' + nodup_bam} \
            ${'--chrsz ' + chrsz} \
            ${'--tss ' + tss}
    }
    output {
        File tss_plot = glob('*.tss_enrich.png')[0]
        File tss_large_plot = glob('*.large_tss_enrich.png')[0]
        File tss_enrich_qc = glob('*.tss_enrich.qc')[0]
        Float tss_enrich = read_float(tss_enrich_qc)
    }
    runtime {
        cpu : 1
        memory : '8 GB'
        time : 4
        disks : 'local-disk 150 SSD'

        docker : runtime_environment.docker
        singularity : runtime_environment.singularity
        conda : runtime_environment.conda
    }
}

task fraglen_stat_pe {
    # for PE only
    input {
        File? nodup_bam
        String? picard_java_heap

        # runtime environment
        RuntimeEnvironment runtime_environment
    }
    Float input_file_size_gb = size(nodup_bam, "G")
    Float mem_gb = 8.0
    Float picard_java_heap_factor = 0.9

    command {
        set -e
        python3 $(which encode_task_fraglen_stat_pe.py) \
            ${'--nodup-bam ' + nodup_bam} \
            ${'--picard-java-heap ' + if defined(picard_java_heap) then picard_java_heap else (round(mem_gb * picard_java_heap_factor) + 'G')}
    }
    output {
        File nucleosomal_qc = glob('*nucleosomal.qc')[0]
        File fraglen_dist_plot = glob('*fraglen_dist.png')[0]
    }
    runtime {
        cpu : 1
        memory : '${mem_gb} GB'
        time : 6
        disks : 'local-disk 150 SSD'

        docker : runtime_environment.docker
        singularity : runtime_environment.singularity
        conda : runtime_environment.conda
    }
}

task gc_bias {
    input {
        File? nodup_bam
        File ref_fa

        String? picard_java_heap

        # runtime environment
        RuntimeEnvironment runtime_environment
    }
    Float mem_factor = 0.3
    Float input_file_size_gb = size(nodup_bam, "G")
    Float mem_gb = 4.0 + mem_factor * input_file_size_gb
    Float picard_java_heap_factor = 0.9

    command {
        set -e
        python3 $(which encode_task_gc_bias.py) \
            ${'--nodup-bam ' + nodup_bam} \
            ${'--ref-fa ' + ref_fa} \
            ${'--picard-java-heap ' + if defined(picard_java_heap) then picard_java_heap else (round(mem_gb * picard_java_heap_factor) + 'G')}
    }
    output {
        File gc_plot = glob('*.gc_plot.png')[0]
        File gc_log = glob('*.gc.txt')[0]
    }
    runtime {
        cpu : 1
        memory : '${mem_gb} GB'
        time : 6
        disks : 'local-disk 250 SSD'

        docker : runtime_environment.docker
        singularity : runtime_environment.singularity
        conda : runtime_environment.conda
    }
}

task compare_signal_to_roadmap {
    input {
        File? pval_bw
        File? dnase
        File? reg2map_bed
        File? reg2map
        File? roadmap_meta

        # runtime environment
        RuntimeEnvironment runtime_environment
    }
    command {
        set -e
        python3 $(which encode_task_compare_signal_to_roadmap.py) \
            ${'--bigwig ' + pval_bw} \
            ${'--dnase ' + dnase} \
            ${'--reg2map-bed ' + reg2map_bed} \
            ${'--reg2map ' + reg2map} \
            ${'--roadmap-meta ' + roadmap_meta}
    }
    output {
        File roadmap_compare_plot = glob('*roadmap_compare_plot.png')[0]
        File roadmap_compare_log = glob('*roadmap_compare.log')[0]
    }
    runtime {
        cpu : 1
        memory : '8 GB'
        time : 4
        disks : 'local-disk 100 SSD'

        docker : runtime_environment.docker
        singularity : runtime_environment.singularity
        conda : runtime_environment.conda
    }
}

task qc_report {
    input {
        String pipeline_ver
        String title
        String description
        String? genome
        # workflow params
        Int multimapping
        Array[Boolean] paired_ends
        String pipeline_type
        String aligner
        Boolean no_dup_removal
        String peak_caller
        Int cap_num_peak
        Float idr_thresh
        Float pval_thresh
        Int xcor_subsample_reads
        # QCs
        Array[File] frac_mito_qcs
        Array[File] samstat_qcs
        Array[File] nodup_samstat_qcs
        Array[File] dup_qcs
        Array[File] lib_complexity_qcs
        Array[File] xcor_plots
        Array[File] xcor_scores
        File? jsd_plot
        Array[File] jsd_qcs
        Array[File] idr_plots
        Array[File] idr_plots_pr
        File? idr_plot_ppr
        Array[File] frip_qcs
        Array[File] frip_qcs_pr1
        Array[File] frip_qcs_pr2
        File? frip_qc_pooled
        File? frip_qc_ppr1
        File? frip_qc_ppr2
        Array[File] frip_idr_qcs
        Array[File] frip_idr_qcs_pr
        File? frip_idr_qc_ppr
        Array[File] frip_overlap_qcs
        Array[File] frip_overlap_qcs_pr
        File? frip_overlap_qc_ppr
        File? idr_reproducibility_qc
        File? overlap_reproducibility_qc

        Array[File] annot_enrich_qcs
        Array[File] tss_enrich_qcs
        Array[File] tss_large_plots
        Array[File] roadmap_compare_plots
        Array[File] fraglen_dist_plots
        Array[File] fraglen_nucleosomal_qcs
        Array[File] gc_plots
        Array[File] preseq_plots
        Array[File] picard_est_lib_size_qcs

        Array[File] peak_region_size_qcs
        Array[File] peak_region_size_plots
        Array[File] num_peak_qcs

        File? idr_opt_peak_region_size_qc
        File? idr_opt_peak_region_size_plot
        File? idr_opt_num_peak_qc

        File? overlap_opt_peak_region_size_qc
        File? overlap_opt_peak_region_size_plot
        File? overlap_opt_num_peak_qc

        File? qc_json_ref

        # runtime environment
        RuntimeEnvironment runtime_environment
    }
    command {
        set -e
        python3 $(which encode_task_qc_report.py) \
            --pipeline-prefix atac \
            ${'--pipeline-ver ' + pipeline_ver} \
            ${"--title '" + sub(title,"'","_") + "'"} \
            ${"--desc '" + sub(description,"'","_") + "'"} \
            ${'--genome ' + genome} \
            ${'--multimapping ' + multimapping} \
            --paired-ends ${sep=' ' paired_ends} \
            --pipeline-type ${pipeline_type} \
            --aligner ${aligner} \
            ${if (no_dup_removal) then '--no-dup-removal ' else ''} \
            --peak-caller ${peak_caller} \
            ${'--cap-num-peak ' + cap_num_peak} \
            --idr-thresh ${idr_thresh} \
            --pval-thresh ${pval_thresh} \
            --xcor-subsample-reads ${xcor_subsample_reads} \
            --frac-mito-qcs ${sep='_:_' frac_mito_qcs} \
            --samstat-qcs ${sep='_:_' samstat_qcs} \
            --nodup-samstat-qcs ${sep='_:_' nodup_samstat_qcs} \
            --dup-qcs ${sep='_:_' dup_qcs} \
            --lib-complexity-qcs ${sep='_:_' lib_complexity_qcs} \
            --xcor-plots ${sep='_:_' xcor_plots} \
            --xcor-scores ${sep='_:_' xcor_scores} \
            --idr-plots ${sep='_:_' idr_plots} \
            --idr-plots-pr ${sep='_:_' idr_plots_pr} \
            ${'--jsd-plot ' + jsd_plot} \
            --jsd-qcs ${sep='_:_' jsd_qcs} \
            ${'--idr-plot-ppr ' + idr_plot_ppr} \
            --frip-qcs ${sep='_:_' frip_qcs} \
            --frip-qcs-pr1 ${sep='_:_' frip_qcs_pr1} \
            --frip-qcs-pr2 ${sep='_:_' frip_qcs_pr2} \
            ${'--frip-qc-pooled ' + frip_qc_pooled} \
            ${'--frip-qc-ppr1 ' + frip_qc_ppr1} \
            ${'--frip-qc-ppr2 ' + frip_qc_ppr2} \
            --frip-idr-qcs ${sep='_:_' frip_idr_qcs} \
            --frip-idr-qcs-pr ${sep='_:_' frip_idr_qcs_pr} \
            ${'--frip-idr-qc-ppr ' + frip_idr_qc_ppr} \
            --frip-overlap-qcs ${sep='_:_' frip_overlap_qcs} \
            --frip-overlap-qcs-pr ${sep='_:_' frip_overlap_qcs_pr} \
            ${'--frip-overlap-qc-ppr ' + frip_overlap_qc_ppr} \
            ${'--idr-reproducibility-qc ' + idr_reproducibility_qc} \
            ${'--overlap-reproducibility-qc ' + overlap_reproducibility_qc} \
            --annot-enrich-qcs ${sep='_:_' annot_enrich_qcs} \
            --tss-enrich-qcs ${sep='_:_' tss_enrich_qcs} \
            --tss-large-plots ${sep='_:_' tss_large_plots} \
            --roadmap-compare-plots ${sep='_:_' roadmap_compare_plots} \
            --fraglen-dist-plots ${sep='_:_' fraglen_dist_plots} \
            --fraglen-nucleosomal-qcs ${sep='_:_' fraglen_nucleosomal_qcs} \
            --gc-plots ${sep='_:_' gc_plots} \
            --preseq-plots ${sep='_:_' preseq_plots} \
            --picard-est-lib-size-qcs ${sep='_:_' picard_est_lib_size_qcs} \
            --peak-region-size-qcs ${sep='_:_' peak_region_size_qcs} \
            --peak-region-size-plots ${sep='_:_' peak_region_size_plots} \
            --num-peak-qcs ${sep='_:_' num_peak_qcs} \
            ${'--idr-opt-peak-region-size-qc ' + idr_opt_peak_region_size_qc} \
            ${'--idr-opt-peak-region-size-plot ' + idr_opt_peak_region_size_plot} \
            ${'--idr-opt-num-peak-qc ' + idr_opt_num_peak_qc} \
            ${'--overlap-opt-peak-region-size-qc ' + overlap_opt_peak_region_size_qc} \
            ${'--overlap-opt-peak-region-size-plot ' + overlap_opt_peak_region_size_plot} \
            ${'--overlap-opt-num-peak-qc ' + overlap_opt_num_peak_qc} \
            --out-qc-html qc.html \
            --out-qc-json qc.json \
            ${'--qc-json-ref ' + qc_json_ref}
    }
    output {
        File report = glob('*qc.html')[0]
        File qc_json = glob('*qc.json')[0]
        Boolean qc_json_ref_match = read_string('qc_json_ref_match.txt')=='True'
    }
    runtime {
        cpu : 1
        memory : '4 GB'
        time : 4
        disks : 'local-disk 50 SSD'

        docker : runtime_environment.docker
        singularity : runtime_environment.singularity
        conda : runtime_environment.conda
    }
}

task read_genome_tsv {
    input {
        File? genome_tsv
        String? null_s

        RuntimeEnvironment runtime_environment
    }
    command <<<
        echo "$(basename ~{genome_tsv})" > genome_name
        # create empty files for all entries
        touch ref_fa bowtie2_idx_tar chrsz gensz blacklist blacklist2
        touch ref_mito_fa
        touch bowtie2_mito_idx_tar
        touch tss tss_enrich # for backward compatibility
        touch dnase prom enh reg2map reg2map_bed roadmap_meta
        touch mito_chr_name
        touch regex_bfilt_peak_chr_name

        python <<CODE
        import os
        with open('~{genome_tsv}','r') as fp:
            for line in fp:
                arr = line.strip('\n').split('\t')
                if arr:
                    key, val = arr
                    with open(key,'w') as fp2:
                        fp2.write(val)
        CODE
    >>>
    output {
        String? genome_name = read_string('genome_name')
        String? ref_fa = if size('ref_fa')==0 then null_s else read_string('ref_fa')
        String? ref_mito_fa = if size('ref_mito_fa')==0 then null_s else read_string('ref_mito_fa')
        String? bowtie2_idx_tar = if size('bowtie2_idx_tar')==0 then null_s else read_string('bowtie2_idx_tar')
        String? bowtie2_mito_idx_tar = if size('bowtie2_mito_idx_tar')==0 then null_s else read_string('bowtie2_mito_idx_tar')
        String? chrsz = if size('chrsz')==0 then null_s else read_string('chrsz')
        String? gensz = if size('gensz')==0 then null_s else read_string('gensz')
        String? blacklist = if size('blacklist')==0 then null_s else read_string('blacklist')
        String? blacklist2 = if size('blacklist2')==0 then null_s else read_string('blacklist2')
        String? mito_chr_name = if size('mito_chr_name')==0 then null_s else read_string('mito_chr_name')
        String? regex_bfilt_peak_chr_name = if size('regex_bfilt_peak_chr_name')==0 then 'chr[\\dXY]+'
            else read_string('regex_bfilt_peak_chr_name')
        String? tss = if size('tss')!=0 then read_string('tss')
            else if size('tss_enrich')!=0 then read_string('tss_enrich') else null_s
        String? dnase = if size('dnase')==0 then null_s else read_string('dnase')
        String? prom = if size('prom')==0 then null_s else read_string('prom')
        String? enh = if size('enh')==0 then null_s else read_string('enh')
        String? reg2map = if size('reg2map')==0 then null_s else read_string('reg2map')
        String? reg2map_bed = if size('reg2map_bed')==0 then null_s else read_string('reg2map_bed')
        String? roadmap_meta = if size('roadmap_meta')==0 then null_s else read_string('roadmap_meta')
    }
    runtime {
        maxRetries : 0
        cpu : 1
        memory : '2 GB'
        time : 4
        disks : 'local-disk 10 SSD'

        docker : runtime_environment.docker
        singularity : runtime_environment.singularity
        conda : runtime_environment.conda
    }
}

task raise_exception {
    input {
        String msg

        # runtime environment
        RuntimeEnvironment runtime_environment
    }
    command {
        echo -e "\n* Error: ${msg}\n" >&2
        exit 2
    }
    output {
        String error_msg = '${msg}'
    }
    runtime {
        maxRetries : 0
        cpu : 1
        memory : '2 GB'
        time : 4
        disks : 'local-disk 10 SSD'

        docker : runtime_environment.docker
        singularity : runtime_environment.singularity
        conda : runtime_environment.conda
    }
}
