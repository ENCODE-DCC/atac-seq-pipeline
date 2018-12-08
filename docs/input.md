Input JSON
==========

We provide two template JSON files for both single ended and paired-end samples. We recommend to use one of these input JSON files instead of that used in the tutorial section. These template JSON files include all parameters of the pipeline with default values defined.

* [template](../examples/template_se.json) for single ended sample
* [template](../examples/template_pe.json) for paired-end sample

An input JSON file includes all input parameters and metadata for running pipelines. 3) and 4) are optional so that our pipeline will use default values if they are not defined. However, 1) and 2) are mandatory.

1)<sup>*</sup> Reference genome (hg38, mm10, hg19, ...).
2)<sup>*</sup> Input data file paths/URIs.
3) Pipeline parameters.
4) Resource settings for jobs.

Let us take a quick look at the following template JSON. A JSON file does not allow comments in it but we added some to help you understand each parameter.
```javascript
{
    ////////// 1) Reference genome //////////

    // Download or build reference genome database and pick a TSV from it.
    "atac.genome_tsv" : "/path_to_genome_data/hg38/hg38.tsv",

    ////////// 2) Input data files paths/URIs //////////

    // sample's endedness
    "atac.paired_end" : true,

    // If you start from FASTQs then define these, otherwise remove from this file.
    // You can define up to 6 replicates.
    // FASTQs in an array will be merged after trimming adapters.
    // For example, 
    // "rep1_R1_L1.fastq.gz", "rep1_R1_L2.fastq.gz" and "rep1_R1_L3.fastq.gz" will be merged together.
    "atac.fastqs_rep1_R1" : [ "rep1_R1_L1.fastq.gz", "rep1_R1_L2.fastq.gz", "rep1_R1_L3.fastq.gz" ],
    "atac.fastqs_rep1_R2" : [ "rep1_R2_L1.fastq.gz", "rep1_R2_L2.fastq.gz", "rep1_R2_L3.fastq.gz" ],
    "atac.fastqs_rep2_R1" : [ "rep2_R1_L1.fastq.gz", "rep2_R1_L2.fastq.gz" ],
    "atac.fastqs_rep2_R2" : [ "rep2_R2_L1.fastq.gz", "rep2_R2_L2.fastq.gz" ],

    // If you start from BAMs then define these, otherwise remove from this file.
    // You can define up to 6 replicates. The following example array has two replicates.
    "atac.bams" : [
        "raw_rep1.bam",
        "raw_rep2.bam"
    ],

    // If you start from filtered/deduped BAMs then define these, otherwise remove from this file.
    // You can define up to 6 replicates. The following example array has two replicates.
    "atac.nodup_bams" : [
        "nodup_rep1.bam",
        "nodup_rep2.bam"
    ],

    // If you start from TAG-ALIGNs then define these, otherwise remove from this file.
    // You can define up to 6 replicates. The following example array has two replicates.
    "atac.tas" : [
        "rep1.tagAlign.gz",
        "rep2.tagAlign.gz"
    ],

    // You can use auto-detection for adapters.
    // List of adapters can be detected:
    //   AGATCGGAAGAGC (Illumina), CTGTCTCTTATA (Nextera) and TGGAATTCTCGG (smallRNA)
    "atac.auto_detect_adapter" : false,

    // If you don't start from FASTQs, remove these adapter arrays from this file.
    // else if you chooose to use auto-detection for adapters, then remove adapter arrays from this file.
    // Otherwise define adapters for each FASTQ.
    // Adapters should have the same dimension as FASTQs.
    "atac.adapters_rep1_R1" : [ "AATTCCGG", "AATTCCGG", "AATTCCGG" ],
    "atac.adapters_rep2_R2" : [ "AATTCCGG", "AATTCCGG" ],
    "atac.adapters_rep1_R1" : [ "AATTCCGG", "AATTCCGG", "AATTCCGG" ],
    "atac.adapters_rep2_R2" : [ "AATTCCGG", "AATTCCGG" ],

    ////////// 3) Pipeline parameters //////////

    // Pipeline title and description
    "atac.title" : "Example (paired end)",
    "atac.description" : "This is a template input JSON for paired ended sample.",

    // Pipeline type (atac or dnase). DNase-Seq pipeline is not an official ENCODE pipeline yet.
    "atac.pipeline_type" : "atac",

    // Pipeline will not proceed to post alignment steps (peak-calling, ...).
    // You will get QC report for alignment only.
    "atac.align_only" : false,

    // Pipeline will not generate pseudo replicates and will skip all analyses related to them.
    "atac.true_rep_only" : false,

    // cutadapt (trim_adapter) parameters
    "atac.cutadapt_min_trim_len" : 5,
    "atac.cutadapt_err_rate" : 0.1,

    // multimapping reads
    "atac.multimapping" : 0,

    // bowtie2 parameters
    "atac.bowtie2_score_min" : "",

    // Choose a dup marker between picard and sambamba
    // picard is recommended, use sambamba only when picard fails.
    "atac.dup_marker" : "picard",

    // Threshold for mapped reads quality (samtools view -q)
    "atac.mapq_thresh" : 30,

    // Skip dup removal in a BAM filtering stage.
    "atac.no_dup_removal" : false,

    // Regular expression to filter out reads
    // Any read that matches with this reg-ex pattern will be removed from outputs
    "atac.regex_filter_reads" : "chrM",

    // Subsample reads
    // Subsampled reads will be used for all downsteam analyses including peak-calling
    "atac.subsample_reads" : 0,

    // Cross-correlation analysis
    "atac.enable_xcor" : false,

    // Subsample reads for cross-corr. analysis only
    // Subsampled reads will be used for cross-corr. analysis only
    "atac.xcor_subsample_reads" : 25000000,

    // Keep irregular chromosome names
    // Use this for custom genomes without canonical chromosome names (chr1, chrX, ...)
    "atac.keep_irregular_chr_in_bfilt_peak" : false,        

    // Cap number of peaks called from a peak-caller (MACS2)
    "atac.cap_num_peak" : 300000,
    // P-value threshold for MACS2 (macs2 callpeak -p)
    "atac.pval_thresh" : 0.01,
    // Smoothing window for MACS2 (macs2 callpeak --shift -smooth_win/2 --extsize smooth_win)
    "atac.smooth_win" : 150,

    // IDR (irreproducible discovery rate)
    "atac.enable_idr" : false,
    // Threshold for IDR
    "atac.idr_thresh" : 0.1,

    // ATAqC (annotation-based analysis which include TSS enrichment and etc.)
    "atac.disable_ataqc" : false,

    ////////// 4) Resource settings //////////

    // Set of resources defined here is PER REPLICATE.
    // Therefore, total number of cores used will "atac.bowtie2_cpu" x [NUMBER_OF_REPLICATES]
    // because bowtie2 is a bottlenecking job of the pipeline.
    // Use this total number of cores if you manually qsub or sbatch your job (using local mode of our pipeline).
    // "disks" is used for Google Cloud and DNANexus only.

    "atac.trim_adapter_cpu" : 2,
    "atac.trim_adapter_mem_mb" : 12000,
    "atac.trim_adapter_time_hr" : 24,
    "atac.trim_adapter_disks" : "local-disk 100 HDD",

    "atac.bowtie2_cpu" : 4,
    "atac.bowtie2_mem_mb" : 20000,
    "atac.bowtie2_time_hr" : 48,
    "atac.bowtie2_disks" : "local-disk 100 HDD",

    "atac.filter_cpu" : 2,
    "atac.filter_mem_mb" : 20000,
    "atac.filter_time_hr" : 24,
    "atac.filter_disks" : "local-disk 100 HDD",

    "atac.bam2ta_cpu" : 2,
    "atac.bam2ta_mem_mb" : 10000,
    "atac.bam2ta_time_hr" : 6,
    "atac.bam2ta_disks" : "local-disk 100 HDD",

    "atac.spr_mem_mb" : 16000,

    "atac.xcor_cpu" : 2,
    "atac.xcor_mem_mb" : 16000,
    "atac.xcor_time_hr" : 6,
    "atac.xcor_disks" : "local-disk 100 HDD",

    "atac.macs2_mem_mb" : 16000,
    "atac.macs2_time_hr" : 24,
    "atac.macs2_disks" : "local-disk 100 HDD",

    "atac.ataqc_mem_mb" : 16000,
    "atac.ataqc_mem_java_mb" : 16000,
    "atac.ataqc_time_hr" : 24,
    "atac.ataqc_disks" : "local-disk 100 HDD"
}
```

## Reference genome

We currently support 4 genomes. You can also [build a genome database for your own genome](build_genome_database.md).

|genome|source|built from|
|-|-|-|
|hg38|ENCODE|[GRCh38_no_alt_analysis_set_GCA_000001405](https://www.encodeproject.org/files/GRCh38_no_alt_analysis_set_GCA_000001405.15/@@download/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta.gz)|
|mm10|ENCODE|[mm10_no_alt_analysis_set_ENCODE](https://www.encodeproject.org/files/mm10_no_alt_analysis_set_ENCODE/@@download/mm10_no_alt_analysis_set_ENCODE.fasta.gz)|
|hg19|UCSC|[GRCh37/hg19](http://hgdownload.cse.ucsc.edu/goldenpath/hg19/encodeDCC/referenceSequences/male.hg19.fa.gz)|
|mm9|UCSC|[mm9, NCBI Build 37](<http://hgdownload.cse.ucsc.edu/goldenPath/mm9/bigZips/mm9.2bit>)|

Choose one TSV file for `"atac.genome_tsv"` in your input JSON. `[GENOME]` should be `hg38`, `mm10`, `hg19` or `mm9`.

|platform|path/URI|
|-|-|
|Google Cloud Platform|`gs://encode-pipeline-genome-data/[GENOME]_google.tsv`|
|DNANexus (CLI)|`dx://project-BKpvFg00VBPV975PgJ6Q03v6:data/pipeline-genome-data/[GENOME]_dx.tsv`|
|DNANexus (CLI, Azure)|`dx://project-XXXXXXXXXXXXXX:data/pipeline-genome-data/[GENOME]_dx.tsv`|
|DNANExus (Web)|Choose `[GENOME]_dx.tsv` from [here](https://platform.dnanexus.com/projects/BKpvFg00VBPV975PgJ6Q03v6/data/pipeline-genome-data)|
|DNANExus (Web, Azure)|Choose `[GENOME]_dx.tsv` from [here](https://platform.dnanexus.com/projects/XXXXXXXXXXXXXX/data/pipeline-genome-data)|
|Stanford Sherlock|`genome/scg/[GENOME]_scg.tsv`|
|Stanford SCG|`genome/sherlock/[GENOME]_sherlock.tsv`|
|Local/SLURM/SGE|You need to [build a genome database](build_genome_database.md). |

## Input data file

Choose any data type (FASTQ, BAM, nodup/filtered BAM, TAG-ALIGN and PEAK) you want and DO NOT define arrays for other types. For FASTQs and their corresponding adapter arrays, we provide two ways to define them since DNANexus web UI supports up to an 1-dim array. Choose between 3-dim `fastqs` or 1-dim `fastqs_rep[REP_ID]_R[READ_END_ID]` according to your preference. The pipeline supports up to 6 replicates.

* `"atac.fastqs"` : 3-dimensional array with FASTQ file path/URI.
    - 1st dimension: replicate ID
    - 2nd dimension: merge ID (this dimension will be reduced after merging FASTQs)
    - 3rd dimension: endedness ID (0 for SE and 0,1 for PE)
* `"atac.fastqs_rep1_R1"` : Array of FASTQ file to be merged for rep1-R1.
* `"atac.fastqs_rep1_R2"` : Array of FASTQ file to be merged for rep1-R2. Do not define if your FASTQ is single ended.
* `"atac.fastqs_rep2_R1"` : Array of FASTQ file to be merged for rep2-R1. Do not define if you don't have replicate 2.
* `"atac.fastqs_rep2_R2"` : Array of FASTQ file to be merged for rep2-R2. Do not define if you don't have replicate 2.
* `"atac.fastqs_rep3_R1"` : Array of FASTQ file to be merged for rep3-R1. Do not define if you don't have replicate 3.
* `"atac.fastqs_rep3_R2"` : Array of FASTQ file to be merged for rep3-R2. Do not define if you don't have replicate 3.
* `"atac.fastqs_rep4_R1"` : Array of FASTQ file to be merged for rep4-R1. Do not define if you don't have replicate 4.
* `"atac.fastqs_rep4_R2"` : Array of FASTQ file to be merged for rep4-R2. Do not define if you don't have replicate 4.
* `"atac.bams"` : Array of raw (unfiltered) BAM file path/URI.
    - 1st dimension: replicate ID
* `"atac.nodup_bams"` : Array of filtered (deduped) BAM file path/URI.
    - 1st dimension: replicate ID
* `"atac.tas"` : Array of TAG-ALIGN file path/URI.
    - 1st dimension: replicate ID
* `"atac.peaks"` : Array of NARROWPEAK file path/URI.
    - 1st dimension: replicate ID
* `"atac.peaks_pr1"` : Array of NARROWPEAK file path/URI for 1st self pseudo replicate of replicate ID.
    - 1st dimension: replicate ID
* `"atac.peaks_pr2"` : Array of NARROWPEAK file path/URI for 2nd self pseudo replicate of replicate ID.
    - 1st dimension: replicate ID
* `"atac.peak_ppr1"` : NARROWPEAK file path/URI for pooled 1st pseudo replicates.
* `"atac.peak_ppr2"` : NARROWPEAK file path/URI for pooled 2nd pseudo replicates.
* `"atac.peak_pooled"` : NARROWPEAK file path/URI for pooled replicate.

If starting from peaks then always define `"atac.peaks"`. Define `"atac.peaks_pr1"`, `"atac.peaks_pr2"`, `"atac.peak_pooled"`, `"atac.peak_ppr1"` and `"atac.peak_ppr2"` according to the following rules:

```
if num_rep>1:
    if true_rep_only: peak_pooled, 
    else: peaks_pr1[], peaks_pr2[], peak_pooled, peak_ppr1, peak_ppr2
else:
    if true_rep_only: "not the case!"
    else: peaks_pr1[], peaks_pr2[]
```

## Pipeline parameters

1. General

    Pipeline type (ATAC-Seq or DNase-Seq) : The only difference between two types is TN5 shifting for TAG-ALIGN outputs.

    * `"atac.pipeline_type` : `atac` for ATAC-Seq. `dnase` for DNase-Seq. 

    Input data endedness.

    * `"atac.paired_end"` : Set it as `true` if input data are paired end, otherwise `false`.

    Other optional settings.

    * `"atac.align_only"` : (optional) Disable all downstream analysis (peak calling, ...) after mapping.
    * `"atac.multimapping"` : (optional) Multimapping reads.
    * `"atac.true_rep_only"` : (optional) Set it as `true` to disable all analyses (including IDR, naive-overlap and reproducibility QC) related to pseudo replicates. This flag suppresses `"atac.enable_idr"`.
    * `"atac.enable_xcor` : (optional) Enable cross-correlation analysis.

    * `"atac.qc_report.name"` : (optional) Name of sample.
    * `"atac.qc_report.desc"` : (optional) Description for sample.

2. Adapter trimmer settings

    Structure/dimension of `"atac.adapters` must match with that of `"atac.fastqs"`. If no adapters are given then do not define `"atac.adapters"` in `input.json`. If some adapters are known then define them in `"atac.adapters"` and leave other entries empty (`""`) while keeping the same structure/dimension as in `"atac.fastqs"`. All undefined/non-empty adapters will be trimmed without auto detection.
    
    * `"atac.auto_detect_adapter"` : (optional) Set it as `true` to automatically detect/trim adapters for empty entries in `"atac.adapters"`. There will be no auto detection for non-empty entries it. If `"atac.adapters"` is not defined then all adapters will be detected/trimmed for all fastqs.
    * `"atac.cutadapt_min_trim_len"` : (optional) Minimum trim length for `cutadapt -m` (default: 5).
    * `"atac.cutadapt_err_rate"` : (optional) Maximum allowed adapter error rate for `cutadapt -e` (default: 0.1).

3. Bowtie2 settings (remove a prefix `atac.` for DNANexus CLI).

    * `"atac.bowtie2.score_min"` : (optional) Min. acceptable alignment score function w.r.t read length.

4. Filter/dedup (post-alignment) settings (remove a prefix `atac.` for DNANexus CLI).

    * `"atac.filter.dup_marker"` : (optional) Dup marker. Choose between `picard` (default) and `sambamba`.
    * `"atac.filter.mapq_thresh"` : (optional) Threshold for low MAPQ reads removal (default: 30).
    * `"atac.filter.no_dup_removal"` : (optional) No dup reads removal when filtering BAM.

5. BAM-2-TAGALIGN settings (remove a prefix `atac.` for DNANexus CLI).

    Pipeline filters out chrM reads by default.

    * `"atac.bam2ta.regex_grep_v_ta"` : (optional) Perl-style regular expression pattern to remove matching reads from TAGALIGN (default: `chrM`).
    * `"atac.bam2ta.subsample"` : (optional) Number of reads to subsample TAGALIGN. Subsampled TAGALIGN will be used for all downstream analysis (MACS2, IDR, naive-overlap).

6. Cross correlation analysis settings (remove a prefix `atac.` for DNANexus CLI).

    * `"atac.xcor.subsample"` : (optional) Number of reads to subsample TAGALIGN.

7. MACS2 settings

    **DO NOT DEFINE MACS2 PARAMETERS IN `"atac.macs2"` SCOPE**. All MACS2 parameters must be defined in `"atac"` scope.

    * `"atac.cap_num_peak"` : (optional) Cap number of raw peaks called from MACS2.
    * `"atac.pval_thresh"` : (optional) P-value threshold  (default: 0.01).
    * `"atac.smooth_win"` : (optional) Size of smoothing window (default: 150).

8. IDR settings

    **DO NOT DEFINE IDR PARAMETERS IN `"atac.idr"` SCOPE**. All IDR parameters must be defined in `"atac"` scope.

    * `"atac.enable_idr"` : (optional) Set it as `true` to enable IDR on raw peaks.
    * `"atac.idr_thresh"` : (optional) IDR threshold (default: 0.05).

9. ATAQC (annotation based analysis) settings

    * `"atac.disable_ataqc"` : (optional) Set it as `true` to disable ATAQC.

## Resource

**RESOURCES DEFINED IN AN INPUT JSON ARE PER TASK**. For example, if you have FASTQs for 2 replicates (2 tasks) and set `cpu` for `bowtie2` task as 4 then total number of cpu cores to map FASTQs is 2 x 4 = 8.

CPU (`cpu`), memory (`mem_mb`) settings are used for submitting jobs to cluster engines (SGE and SLURM) and Cloud platforms (Google Cloud Platform, AWS, ...). VM instance type on cloud platforms will be automatically chosen according to each task's `cpu` and `mem_mb`. Number of cores for tasks without `cpu` parameter is fixed at 1.

* `"atac.trim_adapter_cpu"` : (optional) Number of cores for `trim_adapter` (default: 2).
* `"atac.bowtie2_cpu"` : (optional) Number of cores for `bowtie2` (default: 4).
* `"atac.filter_cpu"` : (optional) Number of cores for `filter` (default: 2).
* `"atac.bam2ta_cpu"` : (optional) Number of cores for `bam2ta` (default: 2).
* `"atac.xcor_cpu"` : (optional) Number of cores for `xcor` (default: 2).
* `"atac.trim_adapter_mem_mb"` : (optional) Max. memory limit in MB for `trim_adapter` (default: 10000).
* `"atac.bowtie2_mem_mb"` : (optional) Max. memory limit in MB for `bowtie2` (default: 20000).
* `"atac.filter_mem_mb"` : (optional) Max. memory limit in MB for `filter` (default: 20000).
* `"atac.bam2ta_mem_mb"` : (optional) Max. memory limit in MB for `bam2ta` (default: 10000).
* `"atac.spr_mem_mb"` : (optional) Max. memory limit in MB for `spr` (default: 12000).
* `"atac.xcor_mem_mb"` : (optional) Max. memory limit in MB for `xcor` (default: 10000).
* `"atac.macs2_mem_mb"` : (optional) Max. memory limit in MB for `macs2` (default: 16000).
* `"atac.ataqc_mem_mb"` : (optional) Max. memory limit in MB for `ATAQC` (default: 16000).
* `"atac.ataqc_mem_java_mb"` : (optional) Max. JAVA heap limit in MB for `ATAQC` (default: 16000).

Disks (`disks`) is used for Cloud platforms (Google Cloud Platforms, AWS, ...).

* `"atac.trim_adapter_disks"` : (optional) Disks for `trim_adapter` (default: "local-disk 100 HDD").
* `"atac.bowtie2_disks"` : (optional) Disks for `bowtie2` (default: "local-disk 100 HDD").
* `"atac.filter_disks"` : (optional) Disks for `filter` (default: "local-disk 100 HDD").
* `"atac.bam2ta_disks"` : (optional) Disks for `bam2ta` (default: "local-disk 100 HDD").
* `"atac.xcor_disks"` : (optional) Disks for `xcor` (default: "local-disk 100 HDD").
* `"atac.macs2_disks"` : (optional) Disks for `macs2` (default: "local-disk 100 HDD").

Walltime (`time`) settings (for SGE and SLURM only).

* `"atac.trim_adapter_time_hr"` : (optional) Walltime for `trim_adapter` (default: 24).
* `"atac.bowtie2_time_hr"` : (optional) Walltime for `bowtie2` (default: 48).
* `"atac.filter_time_hr"` : (optional) Walltime for `filter` (default: 24).
* `"atac.bam2ta_time_hr"` : (optional) Walltime for `bam2ta` (default: 6).
* `"atac.xcor_time_hr"` : (optional) Walltime for `xcor` (default: 6).
* `"atac.macs2_time_hr"` : (optional) Walltime for `macs2` (default: 24).
* `"atac.ataqc_time_hr"` : (optional) Walltime for `ATAQC` (default: 24).

