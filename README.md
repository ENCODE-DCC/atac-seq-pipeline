ENCODE ATAC-seq pipeline
===================================================

# Directories
* `backends/` : Backend configuration files (`.conf`)
* `workflow_opts/` : Workflow option files (`.json`)
* `examples/` : input JSON examples (SE and PE)
* `genome/` : genome data TSV files
* `src/` : Python script for each task in WDL
* `installers/` : dependency/genome data installers for systems (Local, SGE and SLURM) without docker support
* `docker_image/` : Dockerfile and MySQL DB initialization script

# Usage

See [Usage](https://github.com/encode-dcc/wdl-pipelines/blob/master/USAGE.md).

# Input JSON for `atac.wdl`

Optional parameters and flags are marked with `?`.

1) Reference genome

    Currently supported genomes:

    * hg38: ENCODE [GRCh38_no_alt_analysis_set_GCA_000001405](https://www.encodeproject.org/files/GRCh38_no_alt_analysis_set_GCA_000001405.15/@@download/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta.gz)
    * mm10: ENCODE [mm10_no_alt_analysis_set_ENCODE](https://www.encodeproject.org/files/mm10_no_alt_analysis_set_ENCODE/@@download/mm10_no_alt_analysis_set_ENCODE.fasta.gz)
    * hg19: ENCODE [GRCh37/hg19](http://hgdownload.cse.ucsc.edu/goldenpath/hg19/encodeDCC/referenceSequences/male.hg19.fa.gz)
    * mm9: [mm9, NCBI Build 37](http://hgdownload.cse.ucsc.edu/goldenPath/mm9/bigZips/mm9.2bit)

    This TSV file has all genome specific data parameters and file path/URIs. Choose one of TSVs in `genome` directory.

    * `"atac.genome_tsv"` : TSV file path/URI.

2) Input genome data files

    Choose any genome data type you want to start with and do not define others. For FASTQs and their adapters, we provide two ways to define them since DNANexus website UI only supports 1-dim array as inputs. Choose between `fastqs` or `fastqs_rep[REP_ID]_R[READ_END_ID]` for your preference. The pipeline supports up to 4 replicates.

    * `"atac.fastqs"`? : 3-dimensional array with FASTQ file path/URI.
        - 1st dimension: replicate ID
        - 2nd dimension: merge ID (this dimension will be reduced after merging FASTQs)
        - 3rd dimension: endedness ID (0 for SE and 0,1 for PE)
    * `"atac.fastqs_rep1_R1"`? : Array of FASTQ file to be merged for rep1-R1.
    * `"atac.fastqs_rep1_R2"`? : Array of FASTQ file to be merged for rep1-R2. Do not define if your FASTQ is single ended.
    * `"atac.fastqs_rep2_R1"`? : Array of FASTQ file to be merged for rep2-R1. Do not define if you don't have replicate 2.
    * `"atac.fastqs_rep2_R2"`? : Array of FASTQ file to be merged for rep2-R2. Do not define if you don't have replicate 2.
    * `"atac.fastqs_rep3_R1"`? : Array of FASTQ file to be merged for rep3-R1. Do not define if you don't have replicate 3.
    * `"atac.fastqs_rep3_R2"`? : Array of FASTQ file to be merged for rep3-R2. Do not define if you don't have replicate 3.
    * `"atac.fastqs_rep4_R1"`? : Array of FASTQ file to be merged for rep4-R1. Do not define if you don't have replicate 4.
    * `"atac.fastqs_rep4_R2"`? : Array of FASTQ file to be merged for rep4-R2. Do not define if you don't have replicate 4.
    * `"atac.bams"`? : Array of raw (unfiltered) BAM file path/URI.
        - 1st dimension: replicate ID
    * `"atac.nodup_bams"`? : Array of filtered (deduped) BAM file path/URI.
        - 1st dimension: replicate ID
    * `"atac.tas"`? : Array of TAG-ALIGN file path/URI.
        - 1st dimension: replicate ID
    * `"atac.peaks"`? : Array of NARROWPEAK file path/URI.
        - 1st dimension: replicate ID
    * `"atac.peaks_pr1"`? : Array of NARROWPEAK file path/URI for 1st self pseudo replicate of replicate ID.
        - 1st dimension: replicate ID
    * `"atac.peaks_pr2"`? : Array of NARROWPEAK file path/URI for 2nd self pseudo replicate of replicate ID.
        - 1st dimension: replicate ID
    * `"atac.peak_ppr1"`? : NARROWPEAK file path/URI for pooled 1st pseudo replicates.
    * `"atac.peak_ppr2"`? : NARROWPEAK file path/URI for pooled 2nd pseudo replicates.
    * `"atac.peak_pooled"`? : NARROWPEAK file path/URI for pooled replicate.

    If starting from peaks then always define `"atac.peaks"`. Define `"atac.peaks_pr1"`, `"atac.peaks_pr2"`, `"atac.peak_pooled"`, `"atac.peak_ppr1"` and `"atac.peak_ppr2"` according to the following rules:

    ```
    if num_rep>1:
        if true_rep_only: peak_pooled, 
        else: peaks_pr1[], peaks_pr2[], peak_pooled, peak_ppr1, peak_ppr2
    else:
        if true_rep_only: "not the case!"
        else: peaks_pr1[], peaks_pr2[]
    ```

3) Pipeline settings

    Pipeline type (ATAC-Seq or DNase-Seq) : The only difference between two types is TN5 shifting.

    * `"atac.pipeline_type` : `atac` for ATAC-Seq. `dnase` for DNase-Seq. 

    Input data endedness.

    * `"atac.paired_end"` : Set it as `true` if input data are paired end, otherwise `false`.

    Other important settings.

    * `"atac.align_only"`? : Disable all downstream analysis after mapping.
    * `"atac.multimapping"`? : Multimapping reads.
    * `"atac.true_rep_only"`? : Set it as `true` to disable all analyses (including IDR, naive-overlap and reproducibility QC) related to pseudo replicates. This flag suppresses `"atac.enable_idr"`.
    * `"atac.disable_xcor`? : Disable cross-correlation analysis.

4) Adapter trimmer settings

    Structure/dimension of `"atac.adapters` must match with that of `"atac.fastqs"`. If no adapters are given then do not define `"atac.adapters"` in `input.json`. If some adapters are known then define them in `"atac.adapters"` and leave other entries empty (`""`) while keeping the same structure/dimension as in `"atac.fastqs"`. All undefined/non-empty adapters will be trimmed without auto detection.
    
    * `"atac.trim_adapter.auto_detect_adapter"` : Set it as `true` to automatically detect/trim adapters for empty entries in `"atac.adapters"`. There will be no auto detection for non-empty entries it. If `"atac.adapters"`t is not defined then all adapters will be detected/trimmed for all fastqs.
    * `"atac.trim_adapter.min_trim_len"`? : Minimum trim length for `cutadapt -m`.
    * `"atac.trim_adapter.err_rate"`? : Maximum allowed adapter error rate for `cutadapt -e`.

5) Bowtie2 settings

    * `"atac.bowtie2.score_min"`? : Min. acceptable alignment score function w.r.t read length.

6) Filter/dedup (post-alignment) settings

    * `"atac.filter.dup_marker"`? : Dup marker. Choose between `picard` (default) and `sambamba`.
    * `"atac.filter.mapq_thresh"`? : Threshold for low MAPQ reads removal.
    * `"atac.filter.no_dup_removal"`? : No dup reads removal when filtering BAM.

7) BAM-2-TAGALIGN settings

    Pipeline filters out chrM reads by default.

    * `"atac.bam2ta.regex_grep_v_ta"`? : Perl-style regular expression pattern to remove matching reads from TAGALIGN (default: `chrM`).
    * `"atac.bam2ta.subsample"`? : Number of reads to subsample TAGALIGN. Subsampled TAGALIGN will be used for all downstream analysis (MACS2, IDR, naive-overlap).

8) Cross correlation analysis settings

    * `"atac.xcor.subsample"`? : Number of reads to subsample TAGALIGN. Only one end (R1) will be used for cross correlation analysis. This will not affect downstream analysis.

9) MACS2 settings

    **DO NOT DEFINE MACS2 PARAMETERS IN `"atac.macs2"` SCOPE**. All MACS2 parameters must be defined in `"atac"` scope.

    * `"atac.cap_num_peak"`? : Cap number of raw peaks called from MACS2.
    * `"atac.pval_thresh"`? : P-value threshold.
    * `"atac.smooth_win"`? : Size of smoothing window.

10) IDR settings

    **DO NOT DEFINE IDR PARAMETERS IN `"atac.idr"` SCOPE**. All IDR parameters must be defined in `"atac"` scope.

    * `"atac.enable_idr"`? : Set it as `true` to enable IDR on raw peaks.
    * `"atac.idr_thresh"`? : IDR threshold.

11) Resources

    **RESOURCES DEFINED IN `input.json` ARE PER TASK**. For example, if you have FASTQs for 2 replicates (2 tasks) and set `cpu` for `bowtie2` task as 4 then total number of cpu cores to map FASTQs is 2 x 4 = 8.

    CPU (`cpu`), memory (`mem_mb`) settings are used for submitting jobs to cluster engines (SGE and SLURM) and Cloud platforms (Google Cloud Platform, AWS, ...). VM instance type on cloud platforms will be automatically chosen according to each task's `cpu` and `mem_mb`. Number of cores for tasks without `cpu` parameter is fixed at 1.

    * `"atac.trim_adapter.cpu"`? : Number of cores for `trim_adapter` (default: 2).
    * `"atac.bowtie2.cpu"`? : Number of cores for `bowtie2` (default: 4).
    * `"atac.filter.cpu"`? : Number of cores for `filter` (default: 2).
    * `"atac.bam2ta.cpu"`? : Number of cores for `bam2ta` (default: 2).
    * `"atac.xcor.cpu"`? : Number of cores for `xcor` (default: 2).
    * `"atac.trim_adapter.mem_mb"`? : Max. memory limit in MB for `trim_adapter` (default: 10000).
    * `"atac.bowtie2.mem_mb"`? : Max. memory limit in MB for `bowtie2` (default: 20000).
    * `"atac.filter.mem_mb"`? : Max. memory limit in MB for `filter` (default: 20000).
    * `"atac.bam2ta.mem_mb"`? : Max. memory limit in MB for `bam2ta` (default: 10000).
    * `"atac.spr.mem_mb"`? : Max. memory limit in MB for `spr` (default: 12000).
    * `"atac.xcor.mem_mb"`? : Max. memory limit in MB for `xcor` (default: 10000).
    * `"atac.macs2_mem_mb"`? : Max. memory limit in MB for `macs2` (default: 16000).

    Disks (`disks`) is used for Cloud platforms (Google Cloud Platforms, AWS, ...).

    * `"atac.trim_adapter.disks"`? : Disks for `trim_adapter` (default: "local-disk 100 HDD").
    * `"atac.bowtie2.disks"`? : Disks for `bowtie2` (default: "local-disk 100 HDD").
    * `"atac.filter.disks"`? : Disks for `filter` (default: "local-disk 100 HDD").
    * `"atac.bam2ta.disks"`? : Disks for `bam2ta` (default: "local-disk 100 HDD").
    * `"atac.xcor.disks"`? : Disks for `xcor` (default: "local-disk 100 HDD").
    * `"atac.macs2_disks"`? : Disks for `macs2` (default: "local-disk 100 HDD").

    Walltime (`time`) settings (for SGE and SLURM only).

    * `"atac.trim_adapter.time_hr"`? : Walltime for `trim_adapter` (default: 24).
    * `"atac.bowtie2.time_hr"`? : Walltime for `bowtie2` (default: 48).
    * `"atac.filter.time_hr"`? : Walltime for `filter` (default: 24).
    * `"atac.bam2ta.time_hr"`? : Walltime for `bam2ta` (default: 6).
    * `"atac.xcor.time_hr"`? : Walltime for `xcor` (default: 6).
    * `"atac.macs2_time_hr"`? : Walltime for `macs2` (default: 24).

12) QC report HTML/JSON

    * `"atac.qc_report.name"`? : Name of sample.
    * `"atac.qc_report.desc"`? : Description for sample.

# Output directory structure

All output filenames keep prefixes from corresponding input filenames. For example. If you have started from `REP1.fastq.gz` and `REP2.fastq.gz` then corresponding alignment log for each replicate has a filename of `REP1.flagstat.qc` and `REP2.flagstat.qc`, respectively.

Some summarizing QC files do not have any prefix. Find `qc.json` and `qc.html` for final QC and HTML report.

1) `DNANexus`: If you choose to use `dxWDL` and run pipelines on DNANexus platform, then output will be stored on the specified output directory without any subdirectories.

2) `Cromwell`: Otherwise `Cromwell` will store outputs for each task under `cromwell-executions/[WORKFLOW_ID]/call-[TASK_NAME]/shard-[IDX]`. For all tasks except `idr` and `overlap`, `[IDX]` means a zero-based index for each replicate but for tasks `idr` and `overlap` it stands for a zero-based index for all possible pair of replicates. For example, you have 3 replicates and all possible combination of two replicates are `[(rep1,rep2), (rep1,rep3), (rep2,rep3)]`. Therefore, `call-idr/shard-2` should be an output directory for the pair of replicate 2 and 3.

REP_PAIR_ID

`cromwell-executions/[WORKFLOW_ID]/`

call-bowtie2/shard-[REP_ID]

For more details, refer to the file table section in an HTML report generated by the pipeline. Files marked as (E) are outputs to be uploaded during ENCODE accession.
```
out                               # root dir. of outputs
│
├ *report.html                    #  HTML report
├ *tracks.json                    #  Tracks datahub (JSON) for WashU browser
├ ENCODE_summary.json             #  Metadata of all datafiles and QC results
│
├ align                           #  mapped alignments
│ ├ rep1                          #   for true replicate 1 
│ │ ├ *.trim.fastq.gz             #    adapter-trimmed fastq
│ │ ├ *.bam                       #    raw bam
│ │ ├ *.nodup.bam (E)             #    filtered and deduped bam
│ │ ├ *.tagAlign.gz               #    tagAlign (bed6) generated from filtered bam
│ │ ├ *.tn5.tagAlign.gz           #    TN5 shifted tagAlign for ATAC pipeline (not for DNase pipeline)
│ │ └ *.*M.tagAlign.gz            #    subsampled tagAlign for cross-corr. analysis
│ ├ rep2                          #   for true repilicate 2
│ ...
│ ├ pooled_rep                    #   for pooled replicate
│ ├ pseudo_reps                   #   for self pseudo replicates
│ │ ├ rep1                        #    for replicate 1
│ │ │ ├ pr1                       #     for self pseudo replicate 1 of replicate 1
│ │ │ ├ pr2                       #     for self pseudo replicate 2 of replicate 1
│ │ ├ rep2                        #    for repilicate 2
│ │ ...                           
│ └ pooled_pseudo_reps            #   for pooled pseudo replicates
│   ├ ppr1                        #    for pooled pseudo replicate 1 (rep1-pr1 + rep2-pr1 + ...)
│   └ ppr2                        #    for pooled pseudo replicate 2 (rep1-pr2 + rep2-pr2 + ...)
│
├ peak                             #  peaks called
│ └ macs2                          #   peaks generated by MACS2
│   ├ rep1                         #    for replicate 1
│   │ ├ *.narrowPeak.gz            #     narrowPeak (p-val threshold = 0.01)
│   │ ├ *.filt.narrowPeak.gz (E)   #     blacklist filtered narrowPeak 
│   │ ├ *.narrowPeak.bb (E)        #     narrowPeak bigBed
│   │ ├ *.narrowPeak.hammock.gz    #     narrowPeak track for WashU browser
│   │ ├ *.pval0.1.narrowPeak.gz    #     narrowPeak (p-val threshold = 0.1)
│   │ └ *.pval0.1.*K.narrowPeak.gz #     narrowPeak (p-val threshold = 0.1) with top *K peaks
│   ├ rep2                         #    for replicate 2
│   ...
│   ├ pseudo_reps                          #   for self pseudo replicates
│   ├ pooled_pseudo_reps                   #   for pooled pseudo replicates
│   ├ overlap                              #   naive-overlapped peaks
│   │ ├ *.naive_overlap.narrowPeak.gz      #     naive-overlapped peak
│   │ └ *.naive_overlap.filt.narrowPeak.gz #     naive-overlapped peak after blacklist filtering
│   └ idr                           #   IDR thresholded peaks
│     ├ true_reps                   #    for replicate 1
│     │ ├ *.narrowPeak.gz           #     IDR thresholded narrowPeak
│     │ ├ *.filt.narrowPeak.gz (E)  #     IDR thresholded narrowPeak (blacklist filtered)
│     │ └ *.12-col.bed.gz           #     IDR thresholded narrowPeak track for WashU browser
│     ├ pseudo_reps                 #    for self pseudo replicates
│     │ ├ rep1                      #    for replicate 1
│     │ ...
│     ├ optimal_set                 #    optimal IDR thresholded peaks
│     │ └ *.filt.narrowPeak.gz (E)  #     IDR thresholded narrowPeak (blacklist filtered)
│     ├ conservative_set            #    optimal IDR thresholded peaks
│     │ └ *.filt.narrowPeak.gz (E)  #     IDR thresholded narrowPeak (blacklist filtered)
│     ├ pseudo_reps                 #    for self pseudo replicates
│     └ pooled_pseudo_reps          #    for pooled pseudo replicate
│
│   
│ 
├ qc                              #  QC logs
│ ├ *IDR_final.qc                 #   Final IDR QC
│ ├ rep1                          #   for true replicate 1
│ │ ├ *.align.log                 #    Bowtie2 mapping stat log
│ │ ├ *.dup.qc                    #    Picard (or sambamba) MarkDuplicate QC log
│ │ ├ *.pbc.qc                    #    PBC QC
│ │ ├ *.nodup.flagstat.qc         #    Flagstat QC for filtered bam
│ │ ├ *M.cc.qc                    #    Cross-correlation analysis score for tagAlign
│ │ ├ *M.cc.plot.pdf/png          #    Cross-correlation analysis plot for tagAlign
│ │ └ *_qc.html/txt               #    ATAQC report
│ ...
│
├ signal                          #  signal tracks
│ ├ macs2                         #   signal tracks generated by MACS2
│ │ ├ rep1                        #    for true replicate 1 
│ │ │ ├ *.pval.signal.bigwig (E)  #     signal track for p-val
│ │ │ └ *.fc.signal.bigwig   (E)  #     signal track for fold change
│ ...
│ └ pooled_rep                    #   for pooled replicate
│ 
├ report                          # files for HTML report
└ meta                            # text files containing md5sum of output files and other metadata
```