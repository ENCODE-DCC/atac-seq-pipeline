version 1.0

workflow atac {
    String pipeline_ver = 'dev-v1.7.2'

    meta {
        author: 'Jin wook Lee (leepc12@gmail.com) at ENCODE-DCC'
        description: 'ATAC-Seq/DNase-Seq pipeline'
        specification_document: 'https://docs.google.com/document/d/1f0Cm4vRyDQDu0bMehHD7P7KOMxTOP-HiNoIvL1VcBt8/edit?usp=sharing'

        caper_docker: 'quay.io/encode-dcc/atac-seq-pipeline:dev-v1.7.2'
        caper_singularity: 'docker://quay.io/encode-dcc/atac-seq-pipeline:dev-v1.7.2'
        croo_out_def: 'https://storage.googleapis.com/encode-pipeline-output-definition/atac.croo.v4.json'
    }
    input {
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
        Array[File?] bams = []
        Array[File?] nodup_bams = []
        Array[File?] tas = []
        Array[File?] peaks = []
        Array[File?] peaks_pr1 = []
        Array[File?] peaks_pr2 = []
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

        # group: peak_calling
        Int cap_num_peak = 300000
        Float pval_thresh = 0.01
        Int smooth_win = 150
        Float idr_thresh = 0.05

        # group: resource_parameter
        Int align_cpu = 4
        Int align_mem_mb = 20000
        Int align_time_hr = 48
        String align_disks = 'local-disk 400 HDD'

        Int filter_cpu = 2
        Int filter_mem_mb = 20000
        Int filter_time_hr = 24
        String filter_disks = 'local-disk 400 HDD'

        Int bam2ta_cpu = 2
        Int bam2ta_mem_mb = 10000
        Int bam2ta_time_hr = 6
        String bam2ta_disks = 'local-disk 100 HDD'

        Int spr_mem_mb = 16000

        Int jsd_cpu = 2
        Int jsd_mem_mb = 12000
        Int jsd_time_hr = 6
        String jsd_disks = 'local-disk 200 HDD'

        Int xcor_cpu = 2
        Int xcor_mem_mb = 16000
        Int xcor_time_hr = 6
        String xcor_disks = 'local-disk 100 HDD'

        Int call_peak_cpu = 1
        Int call_peak_mem_mb = 16000
        Int call_peak_time_hr = 24
        String call_peak_disks = 'local-disk 200 HDD'

        Int macs2_signal_track_mem_mb = 16000
        Int macs2_signal_track_time_hr = 24
        String macs2_signal_track_disks = 'local-disk 400 HDD'

        Int preseq_mem_mb = 16000

        String? filter_picard_java_heap
        String? preseq_picard_java_heap
        String? fraglen_stat_picard_java_heap
        String? gc_bias_picard_java_heap
    }

    parameter_meta {
        title: {
            description: 'Experiment title.',
            group: 'pipeline_metadata',
        }
        description: {
            description: 'Experiment description.',
            group: 'pipeline_metadata',
        }
        genome_tsv: {
            description: 'Reference genome database TSV.',
            group: 'reference_genome',
            help: 'This TSV files includes all genome specific parameters (e.g. reference FASTA, bowtie2 index). You can still invidiaully define any parameters in it. Parameters defined in input JSON will override those defined in genome TSV.'
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
            help: 'Setting this on means that all replicates are paired ended. For mixed samples, use atac.paired_ends array instead.'
        }
        paired_ends: {
            description: 'Sequencing endedness array (for mixed SE/PE datasets).',
            group: 'input_genomic_data',
            help: 'Whether each biological replicate is paired ended or not.'
        }
        fastqs_rep1_R1: {
            description: 'Read1 FASTQs to be merged for a biological replicate 1.',
            group: 'input_genomic_data',
            help: 'Define if you want to start pipeline from FASTQs files. Pipeline can start from any type of inputs (e.g. FASTQs, BAMs, ...). Choose one type and fill paramters for that type and leave other undefined. Especially for FASTQs, we have individual variable for each biological replicate to allow FASTQs of technical replicates can be merged. Make sure that they are consistent with read2 FASTQs (atac.fastqs_rep1_R2). These FASTQs are usually technical replicates to be merged.'
        }
        fastqs_rep1_R1: {
            description: 'Read2 FASTQs to be merged for a biological replicate 1.',
            group: 'input_genomic_data',
            help: 'Make sure that they are consistent with read1 FASTQs (atac.fastqs_rep1_R1). These FASTQs are usually technical replicates to be merged.'
        }
        fastqs_rep2_R1: {
            description: 'Read1 FASTQs to be merged for a biological replicate 2.',
            group: 'input_genomic_data',
            help: 'Make sure that they are consistent with read2 FASTQs (atac.fastqs_rep2_R2). These FASTQs are usually technical replicates to be merged.'
        }
        fastqs_rep2_R1: {
            description: 'Read2 FASTQs to be merged for a biological replicate 2.',
            group: 'input_genomic_data',
            help: 'Make sure that they are consistent with read1 FASTQs (atac.fastqs_rep2_R1). These FASTQs are usually technical replicates to be merged.'
        }
        fastqs_rep3_R1: {
            description: 'Read1 FASTQs to be merged for a biological replicate 3.',
            group: 'input_genomic_data',
            help: 'Make sure that they are consistent with read2 FASTQs (atac.fastqs_rep3_R2). These FASTQs are usually technical replicates to be merged.'
        }
        fastqs_rep3_R1: {
            description: 'Read2 FASTQs to be merged for a biological replicate 3.',
            group: 'input_genomic_data',
            help: 'Make sure that they are consistent with read1 FASTQs (atac.fastqs_rep3_R1). These FASTQs are usually technical replicates to be merged.'
        }
        fastqs_rep4_R1: {
            description: 'Read1 FASTQs to be merged for a biological replicate 4.',
            group: 'input_genomic_data',
            help: 'Make sure that they are consistent with read2 FASTQs (atac.fastqs_rep4_R2). These FASTQs are usually technical replicates to be merged.'
        }
        fastqs_rep4_R1: {
            description: 'Read2 FASTQs to be merged for a biological replicate 4.',
            group: 'input_genomic_data',
            help: 'Make sure that they are consistent with read1 FASTQs (atac.fastqs_rep4_R1). These FASTQs are usually technical replicates to be merged.'
        }
        fastqs_rep5_R1: {
            description: 'Read1 FASTQs to be merged for a biological replicate 5.',
            group: 'input_genomic_data',
            help: 'Make sure that they are consistent with read2 FASTQs (atac.fastqs_rep5_R2). These FASTQs are usually technical replicates to be merged.'
        }
        fastqs_rep5_R1: {
            description: 'Read2 FASTQs to be merged for a biological replicate 5.',
            group: 'input_genomic_data',
            help: 'Make sure that they are consistent with read1 FASTQs (atac.fastqs_rep5_R1). These FASTQs are usually technical replicates to be merged.'
        }
        fastqs_rep6_R1: {
            description: 'Read1 FASTQs to be merged for a biological replicate 6.',
            group: 'input_genomic_data',
            help: 'Make sure that they are consistent with read2 FASTQs (atac.fastqs_rep6_R2). These FASTQs are usually technical replicates to be merged.'
        }
        fastqs_rep6_R1: {
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
            help: 'The only difference of two types is that TN5 shifting of TAG-ALIGN is done for atac. TAG-ALIGN is in 6-col BED format. It is a simplified version of BAM.'
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
        enable_tss_enrich: {
            description: 'Enables annotation enrichment analysis.',
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
            help: 'Can detect three types of adapter sequences. Illumina: AGATCGGAAGAGC, Nextera: CTGTCTCTTATA, smallRNA: TGGAATTCTCGG.'
        }
        adapter: {
            description: 'Adapter for all FASTQs.',
            group: 'adapter_trimming',
            help: 'Define if all FASTQs have the same adapter sequence. Otherwise define adapter sequence for individual FASTQ in atac.adapters_repX_R1 and atac.adapters_repX_R2 instead. Use atac.auto_detect_adapter if you want to detect adapters automatically.'
        }
        adapters_rep1_R1: {
            description: 'Adapter sequences for read1 FASTQs to be merged for a biological replicate 1.',
            group: 'adapter_trimming',
            help: 'Make sure that they are consistent with read2 FASTQs (atac.adapters_rep1_R2). You can combine this with atac.auto_detect_adapter. Pipeline will auto-detect/trim adapter sequences for null entry in this list. e.g. ["AAGGCCTT", null, "AAGGCCTT"].'
        }
        adapters_rep1_R1: {
            description: 'Adapter sequences for read2 FASTQs to be merged for a biological replicate 1.',
            group: 'adapter_trimming',
            help: 'Make sure that they are consistent with read1 FASTQs (atac.adapters_rep1_R1).'
        }
        adapters_rep2_R1: {
            description: 'Adapter sequences for read1 FASTQs to be merged for a biological replicate 2.',
            group: 'adapter_trimming',
            help: 'Make sure that they are consistent with read2 FASTQs (atac.adapters_rep2_R2).'
        }
        adapters_rep2_R1: {
            description: 'Adapter sequences for read2 FASTQs to be merged for a biological replicate 2.',
            group: 'adapter_trimming',
            help: 'Make sure that they are consistent with read1 FASTQs (atac.adapters_rep2_R1).'
        }
        adapters_rep3_R1: {
            description: 'Adapter sequences for read1 FASTQs to be merged for a biological replicate 3.',
            group: 'adapter_trimming',
            help: 'Make sure that they are consistent with read2 FASTQs (atac.adapters_rep3_R2).'
        }
        adapters_rep3_R1: {
            description: 'Adapter sequences for read2 FASTQs to be merged for a biological replicate 3.',
            group: 'adapter_trimming',
            help: 'Make sure that they are consistent with read1 FASTQs (atac.adapters_rep3_R1).'
        }
        adapters_rep4_R1: {
            description: 'Adapter sequences for read1 FASTQs to be merged for a biological replicate 4.',
            group: 'adapter_trimming',
            help: 'Make sure that they are consistent with read2 FASTQs (atac.adapters_rep4_R2).'
        }
        adapters_rep4_R1: {
            description: 'Adapter sequences for read2 FASTQs to be merged for a biological replicate 4.',
            group: 'adapter_trimming',
            help: 'Make sure that they are consistent with read1 FASTQs (atac.adapters_rep4_R1).'
        }
        adapters_rep5_R1: {
            description: 'Adapter sequences for read1 FASTQs to be merged for a biological replicate 5.',
            group: 'adapter_trimming',
            help: 'Make sure that they are consistent with read2 FASTQs (atac.adapters_rep5_R2).'
        }
        adapters_rep5_R1: {
            description: 'Adapter sequences for read2 FASTQs to be merged for a biological replicate 5.',
            group: 'adapter_trimming',
            help: 'Make sure that they are consistent with read1 FASTQs (atac.adapters_rep5_R1).'
        }
        adapters_rep6_R1: {
            description: 'Adapter sequences for read1 FASTQs to be merged for a biological replicate 6.',
            group: 'adapter_trimming',
            help: 'Make sure that they are consistent with read2 FASTQs (atac.adapters_rep6_R2).'
        }
        adapters_rep6_R1: {
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
            description: 'Threshold for low MAPQ reads removal',
            group: 'alignment',
            help: 'Low MAPQ reads are filtered out while filtering BAM.',
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
        align_mem_mb: {
            description: 'Memory (MB) required for task align.',
            group: 'resource_parameter',
            help: 'This will be used for determining instance type for GCP/AWS or job\'s maximum required memory for HPC.'
        }
        align_time_hr: {
            description: 'Walltime (h) required for task align.',
            group: 'resource_parameter',
            help: 'This is for HPCs only. e.g. SLURM, SGE, ...'
        }
        align_disks: {
            description: 'Persistent disk size to store intermediate files and outputs of task align.',
            group: 'resource_parameter',
            help: 'This is for GCP/AWS only.'
        }
        filter_cpu: {
            description: 'Number of cores for task filter.',
            group: 'resource_parameter',
            help: 'Task filter filters raw/unfilterd BAM to get filtered/deduped BAM.'
        }
        filter_mem_mb: {
            description: 'Memory (MB) required for task filter.',
            group: 'resource_parameter',
            help: 'This will be used for determining instance type for GCP/AWS or job\'s maximum required memory for HPC.'
        }
        filter_time_hr: {
            description: 'Walltime (h) required for task filter.',
            group: 'resource_parameter',
            help: 'This is for HPCs only. e.g. SLURM, SGE, ...'
        }
        filter_disks: {
            description: 'Persistent disk size to store intermediate files and outputs of task filter.',
            group: 'resource_parameter',
            help: 'This is for GCP/AWS only.'
        }
        bam2ta_cpu: {
            description: 'Number of cores for task bam2ta.',
            group: 'resource_parameter',
            help: 'Task bam2ta converts filtered/deduped BAM in to TAG-ALIGN (6-col BED) format.'
        }
        bam2ta_mem_mb: {
            description: 'Memory (MB) required for task bam2ta.',
            group: 'resource_parameter',
            help: 'This will be used for determining instance type for GCP/AWS or job\'s maximum required memory for HPC.'
        }
        bam2ta_time_hr: {
            description: 'Walltime (h) required for task bam2ta.',
            group: 'resource_parameter',
            help: 'This is for HPCs only. e.g. SLURM, SGE, ...'
        }
        bam2ta_disks: {
            description: 'Persistent disk size to store intermediate files and outputs of task bam2ta.',
            group: 'resource_parameter',
            help: 'This is for GCP/AWS only.'
        }
        spr_mem_mb: {
            description: 'Memory (MB) required for task spr.',
            group: 'resource_parameter',
            help: 'Task spr generates two pseudo-replicates from each biological replicate. This task uses Unix shuf/sort, which can take significant amount of memory.'
        }
        jsd_cpu: {
            description: 'Number of cores for task jsd.',
            group: 'resource_parameter',
            help: 'Task jsd plots Jensen-Shannon distance and metrics related to it.'
        }
        jsd_mem_mb: {
            description: 'Memory (MB) required for task jsd.',
            group: 'resource_parameter',
            help: 'This will be used for determining instance type for GCP/AWS or job\'s maximum required memory for HPC.'
        }
        jsd_time_hr: {
            description: 'Walltime (h) required for task jsd.',
            group: 'resource_parameter',
            help: 'This is for HPCs only. e.g. SLURM, SGE, ...'
        }
        jsd_disks: {
            description: 'Persistent disk size to store intermediate files and outputs of task jsd.',
            group: 'resource_parameter',
            help: 'This is for GCP/AWS only.'
        }
        xcor_cpu: {
            description: 'Number of cores for task xcor.',
            group: 'resource_parameter',
            help: 'Task xcor does cross-correlation analysis (including a plot) on subsampled TAG-ALIGNs.'
        }
        xcor_mem_mb: {
            description: 'Memory (MB) required for task xcor.',
            group: 'resource_parameter',
            help: 'This will be used for determining instance type for GCP/AWS or job\'s maximum required memory for HPC.'
        }
        xcor_time_hr: {
            description: 'Walltime (h) required for task xcor.',
            group: 'resource_parameter',
            help: 'This is for HPCs only. e.g. SLURM, SGE, ...'
        }
        xcor_disks: {
            description: 'Persistent disk size to store intermediate files and outputs of task xcor.',
            group: 'resource_parameter',
            help: 'This is for GCP/AWS only.'
        }
        call_peak_cpu: {
            description: 'Number of cores for task call_peak.',
            group: 'resource_parameter',
            help: 'Task call_peak call peaks on TAG-ALIGNs by using MACS2 peak caller.'
        }
        call_peak_mem_mb: {
            description: 'Memory (MB) required for task call_peak.',
            group: 'resource_parameter',
            help: 'This will be used for determining instance type for GCP/AWS or job\'s maximum required memory for HPC.'
        }
        call_peak_time_hr: {
            description: 'Walltime (h) required for task call_peak.',
            group: 'resource_parameter',
            help: 'This is for HPCs only. e.g. SLURM, SGE, ...'
        }
        call_peak_disks: {
            description: 'Persistent disk size to store intermediate files and outputs of task call_peak.',
            group: 'resource_parameter',
            help: 'This is for GCP/AWS only.'
        }
        macs2_signal_track_mem_mb: {
            description: 'Memory (MB) required for task macs2_signal_track.',
            group: 'resource_parameter',
            help: 'Task macs2_signal_track_mem_mb uses MACS2 to get signal tracks for p-val (suffixed with .pval.) and fold enrichment (suffixed with .fc.). This will be used for determining instance type for GCP/AWS or job\'s maximum required memory for HPC.'
        }
        macs2_signal_track_time_hr: {
            description: 'Walltime (h) required for task macs2_signal_track.',
            group: 'resource_parameter',
            help: 'This is for HPCs only. e.g. SLURM, SGE, ...'
        }
        macs2_signal_track_disks: {
            description: 'Persistent disk size to store intermediate files and outputs of task macs2_signal_track.',
            group: 'resource_parameter',
            help: 'This is for GCP/AWS only.'
        }
        preseq_mem_mb: {
            description: 'Memory (MB) required for task preseq.',
            group: 'resource_parameter',
            help: 'Task preseq estimates the future yield of a genomic library. This will be used for determining instance type for GCP/AWS or job\'s maximum required memory for HPC.'
        }
        filter_picard_java_heap: {
            description: 'Maximum Java heap (java -Xmx) in task filter.',
            group: 'resource_parameter',
            help: 'Maximum memory for Picard tools MarkDuplicates. If not defined, 90% of filter_mem_mb will be used.'
        }
        preseq_picard_java_heap: {
            description: 'Maximum Java heap (java -Xmx) in task preseq.',
            group: 'resource_parameter',
            help: 'Maximum memory for Picard tools EstimateLibraryComplexity. If not defined, 90% of preseq_mem_mb will be used.'
        }
        fraglen_stat_picard_java_heap: {
            description: 'Maximum Java heap (java -Xmx) in task fraglen_stat_pe (for paired end replicate only).',
            group: 'resource_parameter',
            help: 'Maximum memory for Picard tools CollectInsertSizeMetrics. If not defined, 90% of 8000 MB will be used.'
        }
        gc_bias_picard_java_heap: {
            description: 'Maximum Java heap (java -Xmx) in task gc_bias.',
            group: 'resource_parameter',
            help: 'Maximum memory for Picard tools CollectGcBiasMetrics. If not defined, 90% of 10000 MB will be used.'
        }
    }

    String aligner = 'bowtie2'
    String peak_caller = 'macs2'
    String peak_type = 'narrowPeak'
    
    # read genome data and paths
    if ( defined(genome_tsv) ) {
        call read_genome_tsv { input: genome_tsv = genome_tsv }
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
            msg = 'No FASTQ/BAM/TAG-ALIGN/PEAK defined in your input JSON. Check if your FASTQs are defined as "atac.fastqs_repX_RY". DO NOT MISS suffix _R1 even for single ended FASTQ.'
        }
    }

    # align each replicate
    scatter(i in range(num_rep)) {
        # to override endedness definition for individual replicate
        #     paired_end will override paired_ends[i]
        Boolean paired_end_ = if !defined(paired_end) && i<length(paired_ends) then paired_ends[i]
            else select_first([paired_end])

        Boolean has_input_of_align = i<length(fastqs_R1) && length(fastqs_R1[i])>0
        Boolean has_output_of_align = i<length(bams) && defined(bams[i])
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
                mem_mb = align_mem_mb,
                time_hr = align_time_hr,
                disks = align_disks,
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
                mem_mb = align_mem_mb,
                time_hr = align_time_hr,
                disks = align_disks,
            }
        }

        if ( defined(align.non_mito_samstat_qc) && defined(align_mito.samstat_qc) ) {
            call frac_mito { input:
                non_mito_samstat = align.non_mito_samstat_qc,
                mito_samstat = align_mito.samstat_qc,
            }
        }

        Boolean has_input_of_filter = has_output_of_align || defined(align.bam)
        Boolean has_output_of_filter = i<length(nodup_bams) && defined(nodup_bams[i])
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
                mem_mb = filter_mem_mb,
                picard_java_heap = filter_picard_java_heap,
                time_hr = filter_time_hr,
                disks = filter_disks,
            }
        }
        File? nodup_bam_ = if has_output_of_filter then nodup_bams[i] else filter.nodup_bam

        Boolean has_input_of_bam2ta = has_output_of_filter || defined(filter.nodup_bam)
        Boolean has_output_of_bam2ta = i<length(tas) && defined(tas[i])
        if ( has_input_of_bam2ta && !has_output_of_bam2ta ) {
            call bam2ta { input :
                bam = nodup_bam_,
                disable_tn5_shift = if pipeline_type=='atac' then false else true,
                subsample = subsample_reads,
                paired_end = paired_end_,
                mito_chr_name = mito_chr_name_,

                cpu = bam2ta_cpu,
                mem_mb = bam2ta_mem_mb,
                time_hr = bam2ta_time_hr,
                disks = bam2ta_disks,
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
                mem_mb = filter_mem_mb,
                picard_java_heap = filter_picard_java_heap,
                time_hr = filter_time_hr,
                disks = filter_disks,
            }
            call bam2ta as bam2ta_no_dedup { input :
                bam = filter_no_dedup.nodup_bam,  # output name is nodup but it's not deduped
                disable_tn5_shift = if pipeline_type=='atac' then false else true,
                subsample = 0,
                paired_end = paired_end_,
                mito_chr_name = mito_chr_name_,

                cpu = bam2ta_cpu,
                mem_mb = bam2ta_mem_mb,
                time_hr = bam2ta_time_hr,
                disks = bam2ta_disks,
            }
            # subsample tagalign (non-mito) and cross-correlation analysis
            call xcor { input :
                ta = bam2ta_no_dedup.ta,
                subsample = xcor_subsample_reads,
                paired_end = paired_end_,
                mito_chr_name = mito_chr_name_,

                cpu = xcor_cpu,
                mem_mb = xcor_mem_mb,
                time_hr = xcor_time_hr,
                disks = xcor_disks,
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

                mem_mb = macs2_signal_track_mem_mb,
                disks = macs2_signal_track_disks,
                time_hr = macs2_signal_track_time_hr,
            }
        }

        Boolean has_input_of_call_peak = has_output_of_bam2ta || defined(bam2ta.ta)
        Boolean has_output_of_call_peak = i<length(peaks) && defined(peaks[i])
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
                mem_mb = call_peak_mem_mb,
                disks = call_peak_disks,
                time_hr = call_peak_time_hr,
            }
        }
        File? peak_ = if has_output_of_call_peak then peaks[i] else call_peak.peak

        Boolean has_input_of_spr = has_output_of_bam2ta || defined(bam2ta.ta)
        if ( has_input_of_spr && !align_only && !true_rep_only ) {
            call spr { input :
                ta = ta_,
                paired_end = paired_end_,
                mem_mb = spr_mem_mb,
            }
        }

        Boolean has_input_of_call_peak_pr1 = defined(spr.ta_pr1)
        Boolean has_output_of_call_peak_pr1 = i<length(peaks_pr1) && defined(peaks_pr1[i])
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
                mem_mb = call_peak_mem_mb,
                disks = call_peak_disks,
                time_hr = call_peak_time_hr,
            }
        }
        File? peak_pr1_ = if has_output_of_call_peak_pr1 then peaks_pr1[i]
            else call_peak_pr1.peak

        Boolean has_input_of_call_peak_pr2 = defined(spr.ta_pr2)
        Boolean has_output_of_call_peak_pr2 = i<length(peaks_pr2) && defined(peaks_pr2[i])
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
                mem_mb = call_peak_mem_mb,
                disks = call_peak_disks,
                time_hr = call_peak_time_hr,
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
            }
        }
        # tasks factored out from ATAqC
        Boolean has_input_of_tss_enrich = defined(nodup_bam_) && defined(tss_) && (
            defined(align.read_len) || i<length(read_len) && defined(read_len[i]) )
        if ( enable_tss_enrich && has_input_of_tss_enrich ) {
            call tss_enrich { input :
                read_len = if i<length(read_len) && defined(read_len[i]) then read_len[i]
                    else align.read_len,
                nodup_bam = nodup_bam_,
                tss = tss_,
                chrsz = chrsz_,
            }
        }
        if ( enable_fraglen_stat && paired_end_ && defined(nodup_bam_) ) {
            call fraglen_stat_pe { input :
                nodup_bam = nodup_bam_,
                picard_java_heap = fraglen_stat_picard_java_heap,
            }
        }
        if ( enable_preseq && defined(bam_) ) {
            call preseq { input :
                bam = bam_,
                paired_end = paired_end_,
                mem_mb = preseq_mem_mb,
                picard_java_heap = preseq_picard_java_heap,
            }
        }
        if ( enable_gc_bias && defined(nodup_bam_) && defined(ref_fa_) ) {
            call gc_bias { input :
                nodup_bam = nodup_bam_,
                ref_fa = ref_fa_,
                picard_java_heap = gc_bias_picard_java_heap,
            }
        }
        if ( enable_annot_enrich && defined(ta_) && defined(blacklist_) && defined(dnase_) && defined(prom_) && defined(enh_) ) {
            call annot_enrich { input :
                ta = ta_,
                blacklist = blacklist_,
                dnase = dnase_,
                prom = prom_,
                enh = enh_,
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
        }
    }

    # if there are pr1 TAs for ALL replicates then pool them
    Boolean has_all_inputs_of_pool_ta_pr1 = length(select_all(spr.ta_pr1))==num_rep
    if ( has_all_inputs_of_pool_ta_pr1 && num_rep>1 && !align_only && !true_rep_only ) {
        # pool tagaligns from pseudo replicate 1
        call pool_ta as pool_ta_pr1 { input :
            tas = spr.ta_pr1,
            prefix = 'rep-pr1',
        }
    }

    # if there are pr2 TAs for ALL replicates then pool them
    Boolean has_all_inputs_of_pool_ta_pr2 = length(select_all(spr.ta_pr2))==num_rep
    if ( has_all_inputs_of_pool_ta_pr1 && num_rep>1 && !align_only && !true_rep_only ) {
        # pool tagaligns from pseudo replicate 2
        call pool_ta as pool_ta_pr2 { input :
            tas = spr.ta_pr2,
            prefix = 'rep-pr2',
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
            mem_mb = call_peak_mem_mb,
            disks = call_peak_disks,
            time_hr = call_peak_time_hr,
        }
    }
    File? peak_pooled_ = if has_output_of_call_peak_pooled then peak_pooled
        else call_peak_pooled.peak

    Boolean has_input_of_count_signal_track_pooled = defined(pool_ta.ta_pooled)
    if ( has_input_of_count_signal_track_pooled && enable_count_signal_track && num_rep>1 ) {
        call count_signal_track as count_signal_track_pooled { input :
            ta = pool_ta.ta_pooled,
            chrsz = chrsz_,
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

            mem_mb = macs2_signal_track_mem_mb,
            disks = macs2_signal_track_disks,
            time_hr = macs2_signal_track_time_hr,
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
            mem_mb = jsd_mem_mb,
            time_hr = jsd_time_hr,
            disks = jsd_disks,
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
            mem_mb = call_peak_mem_mb,
            disks = call_peak_disks,
            time_hr = call_peak_time_hr,
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
            mem_mb = call_peak_mem_mb,
            disks = call_peak_disks,
            time_hr = call_peak_time_hr,
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
        }
    }

    # reproducibility QC for overlap/IDR peaks
    if ( !align_only && !true_rep_only && num_rep > 0 ) {
        # reproducibility QC for overlapping peaks
        call reproducibility as reproducibility_overlap { input :
            prefix = 'overlap',
            peaks = select_all(overlap.bfilt_overlap_peak),
            peaks_pr = overlap_pr.bfilt_overlap_peak,
            peak_ppr = overlap_ppr.bfilt_overlap_peak,
            peak_type = peak_type_,
            chrsz = chrsz_,
        }
    }

    if ( !align_only && !true_rep_only && num_rep > 0 && enable_idr ) {
        # reproducibility QC for IDR peaks
        call reproducibility as reproducibility_idr { input :
            prefix = 'idr',
            peaks = select_all(idr.bfilt_idr_peak),
            peaks_pr = idr_pr.bfilt_idr_peak,
            peak_ppr = idr_ppr.bfilt_idr_peak,
            peak_type = peak_type_,
            chrsz = chrsz_,
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
        jsd_qcs = jsd.jsd_qcs,

        frip_qcs = select_all(call_peak.frip_qc),
        frip_qcs_pr1 = select_all(call_peak_pr1.frip_qc),
        frip_qcs_pr2 = select_all(call_peak_pr2.frip_qc),

        frip_qc_pooled = call_peak_pooled.frip_qc,
        frip_qc_ppr1 = call_peak_ppr1.frip_qc,
        frip_qc_ppr2 = call_peak_ppr2.frip_qc,

        idr_plots = select_all(idr.idr_plot),
        idr_plots_pr = idr_pr.idr_plot,
        idr_plot_ppr = idr_ppr.idr_plot,
        frip_idr_qcs = select_all(idr.frip_qc),
        frip_idr_qcs_pr = idr_pr.frip_qc,
        frip_idr_qc_ppr = idr_ppr.frip_qc,
        frip_overlap_qcs = select_all(overlap.frip_qc),
        frip_overlap_qcs_pr = overlap_pr.frip_qc,
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
        Int mem_mb
        Int time_hr
        String disks
    }

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
          echo -e "\n* Error: pipeline dependencies not found." 1>&2
          echo 'Conda users: Did you activate Conda environment (conda activate encode-atac-seq-pipeline)?' 1>&2
          echo '    Or did you install Conda and environment correctly (bash scripts/install_conda_env.sh)?' 1>&2
          echo 'GCP/AWS/Docker users: Did you add --docker flag to Caper command line arg?' 1>&2
          echo 'Singularity users: Did you add --singularity flag to Caper command line arg?' 1>&2
          echo -e "\n" 1>&2
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
                ${'--nth ' + cpu}
        fi

        python3 $(which encode_task_post_align.py) \
            R1/*.fastq.gz $(ls *.bam) \
            ${'--mito-chr-name ' + mito_chr_name} \
            ${'--chrsz ' + chrsz} \
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
        memory : '${mem_mb} MB'
        time : time_hr
        disks : disks
        preemptible: 0
    }
}

task frac_mito {
    input {
        File? non_mito_samstat
        File? mito_samstat
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
        memory : '8000 MB'
        time : 1
        disks : 'local-disk 100 HDD'
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
        Int mem_mb
        String? picard_java_heap
        Int time_hr
        String disks
    }
    Float picard_java_heap_factor = 0.9

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
            ${'--nth ' + cpu} \
            ${'--picard-java-heap ' + if defined(picard_java_heap) then picard_java_heap else (round(mem_mb * picard_java_heap_factor) + 'M')}
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
        memory : '${mem_mb} MB'
        time : time_hr
        disks : disks
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
        Int mem_mb
        Int time_hr
        String disks
    }

    command {
        set -e
        python3 $(which encode_task_bam2ta.py) \
            ${bam} \
            ${if paired_end then '--paired-end' else ''} \
            ${if disable_tn5_shift then '--disable-tn5-shift' else ''} \
            ${'--mito-chr-name ' + mito_chr_name} \
            ${'--subsample ' + subsample} \
            ${'--nth ' + cpu}
    }
    output {
        File ta = glob('*.tagAlign.gz')[0]
    }
    runtime {
        cpu : cpu
        memory : '${mem_mb} MB'
        time : time_hr
        disks : disks
    }
}

task spr {
    input {
        File? ta
        Boolean paired_end

        Int mem_mb
    }

    command {
        set -e
        python3 $(which encode_task_spr.py) \
            ${ta} \
            ${if paired_end then '--paired-end' else ''}
    }
    output {
        File ta_pr1 = glob('*.pr1.tagAlign.gz')[0]
        File ta_pr2 = glob('*.pr2.tagAlign.gz')[0]
    }
    runtime {
        cpu : 1
        memory : '${mem_mb} MB'
        time : 1
        disks : 'local-disk 50 HDD'
    }
}

task pool_ta {
    input {
        Array[File?] tas     # TAG-ALIGNs to be merged
        Int? col             # number of columns in pooled TA
        String? prefix         # basename prefix
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
        memory : '4000 MB'
        time : 1
        disks : 'local-disk 50 HDD'
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
        Int mem_mb
        Int time_hr
        String disks
    }

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
        memory : '${mem_mb} MB'
        time : time_hr
        disks : disks
    }
}

task jsd {
    input {
        Array[File?] nodup_bams
        File? blacklist
        Int mapq_thresh

        Int cpu
        Int mem_mb
        Int time_hr
        String disks
    }

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
        memory : '${mem_mb} MB'
        time : time_hr
        disks : disks
    }
}

task count_signal_track {
    input {
        File? ta             # tag-align
        File chrsz            # 2-col chromosome sizes file
    }
    command {
        set -e
        python3 $(which encode_task_count_signal_track.py) \
            ${ta} \
            ${'--chrsz ' + chrsz}
    }
    output {
        File pos_bw = glob('*.positive.bigwig')[0]
        File neg_bw = glob('*.negative.bigwig')[0]
    }
    runtime {
        cpu : 1
        memory : '8000 MB'
        time : 4
        disks : 'local-disk 50 HDD'
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
        Int mem_mb
        Int time_hr
        String disks
    }
    command {
        set -e

        if [ '${peak_caller}' == 'macs2' ]; then
            python3 $(which encode_task_macs2_atac.py) \
                ${ta} \
                ${'--gensz ' + gensz} \
                ${'--chrsz ' + chrsz} \
                ${'--cap-num-peak ' + cap_num_peak} \
                ${'--pval-thresh '+ pval_thresh} \
                ${'--smooth-win '+ smooth_win}
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
        File bfilt_peak_hammock = glob('*.bfilt.'+peak_type+'.hammock.gz*')[0]
        File bfilt_peak_hammock_tbi = glob('*.bfilt.'+peak_type+'.hammock.gz*')[1]
        File frip_qc = glob('*.frip.qc')[0]
        File peak_region_size_qc = glob('*.peak_region_size.qc')[0]
        File peak_region_size_plot = glob('*.peak_region_size.png')[0]
        File num_peak_qc = glob('*.num_peak.qc')[0]
    }
    runtime {
        cpu : if peak_caller == 'macs2' then 1 else cpu
        memory : '${mem_mb} MB'
        time : time_hr
        disks : disks
        preemptible: 0
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

        Int mem_mb
        Int time_hr
        String disks
    }
    command {
        set -e
        python3 $(which encode_task_macs2_signal_track_atac.py) \
            ${ta} \
            ${'--gensz '+ gensz} \
            ${'--chrsz ' + chrsz} \
            ${'--pval-thresh '+ pval_thresh} \
            ${'--smooth-win '+ smooth_win}
    }
    output {
        File pval_bw = glob('*.pval.signal.bigwig')[0]
        File fc_bw = glob('*.fc.signal.bigwig')[0]
    }
    runtime {
        cpu : 1
        memory : '${mem_mb} MB'
        time : time_hr
        disks : disks
        preemptible: 0
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
        File bfilt_idr_peak_hammock = glob('*.bfilt.'+peak_type+'.hammock.gz*')[0]
        File bfilt_idr_peak_hammock_tbi = glob('*.bfilt.'+peak_type+'.hammock.gz*')[1]
        File idr_plot = glob('*.txt.png')[0]
        File idr_unthresholded_peak = glob('*.txt.gz')[0]
        File idr_log = glob('*.idr*.log')[0]
        File frip_qc = if defined(ta) then glob('*.frip.qc')[0] else glob('null')[0]
    }
    runtime {
        cpu : 1
        memory : '8000 MB'
        time : 1
        disks : 'local-disk 50 HDD'
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
        File bfilt_overlap_peak_hammock = glob('*.bfilt.'+peak_type+'.hammock.gz*')[0]
        File bfilt_overlap_peak_hammock_tbi = glob('*.bfilt.'+peak_type+'.hammock.gz*')[1]
        File frip_qc = if defined(ta) then glob('*.frip.qc')[0] else glob('null')[0]
    }
    runtime {
        cpu : 1
        memory : '4000 MB'
        time : 1
        disks : 'local-disk 50 HDD'
    }
}

task reproducibility {
    input {
        String prefix
        Array[File] peaks # peak files from pair of true replicates
                            # in a sorted order. for example of 4 replicates,
                            # 1,2 1,3 1,4 2,3 2,4 3,4.
                            # x,y means peak file from rep-x vs rep-y
        Array[File]? peaks_pr    # peak files from pseudo replicates
        File? peak_ppr            # Peak file from pooled pseudo replicate.
        String peak_type
        File chrsz            # 2-col chromosome sizes file
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
        File optimal_peak_hammock = glob('*optimal_peak.*.hammock.gz*')[0]
        File optimal_peak_hammock_tbi = glob('*optimal_peak.*.hammock.gz*')[1]
        File conservative_peak = glob('*conservative_peak.*.gz')[0]
        File conservative_peak_bb = glob('*conservative_peak.*.bb')[0]
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
        memory : '4000 MB'
        time : 1
        disks : 'local-disk 50 HDD'
    }
}

task preseq {
    input {
        File? bam
        Boolean paired_end

        Int mem_mb
        String? picard_java_heap
        File? null
    }
    Float picard_java_heap_factor = 0.9

    command {
        set -e
        python3 $(which encode_task_preseq.py) \
            ${if paired_end then '--paired-end' else ''} \
            ${'--bam ' + bam} \
            ${'--picard-java-heap ' + if defined(picard_java_heap) then picard_java_heap else (round(mem_mb * picard_java_heap_factor) + 'M')}
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
        memory : '${mem_mb} MB'
        time : 1
        disks : 'local-disk 100 HDD'
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
        memory : '8000 MB'
        time : 1
        disks : 'local-disk 100 HDD'
    }
}

task tss_enrich {
    input {
        Int? read_len
        File? nodup_bam
        File? tss
        File chrsz
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
        memory : '8000 MB'
        time : 1
        disks : 'local-disk 100 HDD'
    }
}

task fraglen_stat_pe {
    # for PE only
    input {
        File? nodup_bam
        String? picard_java_heap
    }
    Int mem_mb = 8000
    Float picard_java_heap_factor = 0.9

    command {
        set -e
        python3 $(which encode_task_fraglen_stat_pe.py) \
            ${'--nodup-bam ' + nodup_bam} \
            ${'--picard-java-heap ' + if defined(picard_java_heap) then picard_java_heap else (round(mem_mb * picard_java_heap_factor) + 'M')}
    }
    output {
        File nucleosomal_qc = glob('*nucleosomal.qc')[0]
        File fraglen_dist_plot = glob('*fraglen_dist.png')[0]
    }
    runtime {
        cpu : 1
        memory : '${mem_mb} MB'
        time : 6
        disks : 'local-disk 100 HDD'
    }
}

task gc_bias {
    input {
        File? nodup_bam
        File ref_fa

        String? picard_java_heap
    }
    Int mem_mb = 10000
    Float picard_java_heap_factor = 0.9

    command {
        set -e
        python3 $(which encode_task_gc_bias.py) \
            ${'--nodup-bam ' + nodup_bam} \
            ${'--ref-fa ' + ref_fa} \
            ${'--picard-java-heap ' + if defined(picard_java_heap) then picard_java_heap else (round(mem_mb * picard_java_heap_factor) + 'M')}
    }
    output {
        File gc_plot = glob('*.gc_plot.png')[0]
        File gc_log = glob('*.gc.txt')[0]
    }
    runtime {
        cpu : 1
        memory : '${mem_mb} MB'
        time : 6
        disks : 'local-disk 100 HDD'
    }
}

task compare_signal_to_roadmap {
    input {
        File? pval_bw
        File? dnase
        File? reg2map_bed
        File? reg2map
        File? roadmap_meta
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
        memory : '8000 MB'
        time : 1
        disks : 'local-disk 100 HDD'
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
        Array[File]? jsd_qcs
        Array[File] idr_plots
        Array[File]? idr_plots_pr
        File? idr_plot_ppr
        Array[File] frip_qcs
        Array[File] frip_qcs_pr1
        Array[File] frip_qcs_pr2
        File? frip_qc_pooled
        File? frip_qc_ppr1
        File? frip_qc_ppr2
        Array[File] frip_idr_qcs
        Array[File]? frip_idr_qcs_pr
        File? frip_idr_qc_ppr
        Array[File] frip_overlap_qcs
        Array[File]? frip_overlap_qcs_pr
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
    }
    command {
        set -e
        python3 $(which encode_task_qc_report.py) \
            ${'--pipeline-ver ' + pipeline_ver} \
            ${"--title '" + sub(title,"'","_") + "'"} \
            ${"--desc '" + sub(description,"'","_") + "'"} \
            ${'--genome ' + genome} \
            ${'--multimapping ' + multimapping} \
            --paired-ends ${sep=' ' paired_ends} \
            --pipeline-type ${pipeline_type} \
            --aligner ${aligner} \
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
        memory : '4000 MB'
        time : 1
        disks : 'local-disk 50 HDD'
    }
}

task read_genome_tsv {
    input {
        File? genome_tsv
        String? null_s
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
        memory : '4000 MB'
        time : 1
        disks : 'local-disk 50 HDD'
    }
}

task raise_exception {
    input {
        String msg
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
    }
}
