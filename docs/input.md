Input JSON
==========

An input JSON file includes all input parameters and metadata for running pipelines:

1) Reference genome (hg38, mm10, hg19, ...) and genome specific parameters (indices, ...).
2) Input data file paths/URIs (FASTQs, BAMs, TAG-ALIGNs, ...).
3) Pipeline parameters.
4) Resource for instances/jobs.

## For DNANexus CLI users

dxWDL (DNANexus CLI for WDL) does not support definition of task level variables with a prefix `atac.` in an input JSON file. Therefore, `atac.[TASK_NAME].[VAR_NAME]` should be replaced with `[TASK_NAME].[VAR_NAME]`. Simply remove a prefix `atac.` for task level variables. BUT DO NOT REMOVE it for workflow level variables. For example, `atac.qc_report.name` is a task (task `qc_report` in a workflow `atac`) level variable so it should be replaced with `qc_report.name`. But `atac.genome_tsv` is a workflow (`atac`) level variable, so you need to keep it the same. This is the only difference between DNANexus CLI and other platforms.

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
|DNANExus (Web)|Choose `[GENOME]_dx.tsv` from [here](https://platform.dnanexus.com/projects/BKpvFg00VBPV975PgJ6Q03v6/data/pipeline-genome-data)|
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
    * `"atac.disable_xcor` : (optional) Disable cross-correlation analysis.

    * `"atac.qc_report.name"` : (optional) Name of sample.
    * `"atac.qc_report.desc"` : (optional) Description for sample.

2. Adapter trimmer settings

    Structure/dimension of `"atac.adapters` must match with that of `"atac.fastqs"`. If no adapters are given then do not define `"atac.adapters"` in `input.json`. If some adapters are known then define them in `"atac.adapters"` and leave other entries empty (`""`) while keeping the same structure/dimension as in `"atac.fastqs"`. All undefined/non-empty adapters will be trimmed without auto detection.
    
    * `"atac.trim_adapter.auto_detect_adapter"` : (optional) Set it as `true` to automatically detect/trim adapters for empty entries in `"atac.adapters"`. There will be no auto detection for non-empty entries it. If `"atac.adapters"` is not defined then all adapters will be detected/trimmed for all fastqs.
    * `"atac.trim_adapter.min_trim_len"` : (optional) Minimum trim length for `cutadapt -m` (default: 5).
    * `"atac.trim_adapter.err_rate"` : (optional) Maximum allowed adapter error rate for `cutadapt -e` (default: 0.1).

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

* `"atac.trim_adapter.cpu"` : (optional) Number of cores for `trim_adapter` (default: 2).
* `"atac.bowtie2.cpu"` : (optional) Number of cores for `bowtie2` (default: 4).
* `"atac.filter.cpu"` : (optional) Number of cores for `filter` (default: 2).
* `"atac.bam2ta.cpu"` : (optional) Number of cores for `bam2ta` (default: 2).
* `"atac.xcor.cpu"` : (optional) Number of cores for `xcor` (default: 2).
* `"atac.trim_adapter.mem_mb"` : (optional) Max. memory limit in MB for `trim_adapter` (default: 10000).
* `"atac.bowtie2.mem_mb"` : (optional) Max. memory limit in MB for `bowtie2` (default: 20000).
* `"atac.filter.mem_mb"` : (optional) Max. memory limit in MB for `filter` (default: 20000).
* `"atac.bam2ta.mem_mb"` : (optional) Max. memory limit in MB for `bam2ta` (default: 10000).
* `"atac.spr.mem_mb"` : (optional) Max. memory limit in MB for `spr` (default: 12000).
* `"atac.xcor.mem_mb"` : (optional) Max. memory limit in MB for `xcor` (default: 10000).
* `"atac.macs2_mem_mb"` : (optional) Max. memory limit in MB for `macs2` (default: 16000).
* `"atac.ataqc.mem_mb"` : (optional) Max. memory limit in MB for `ATAQC` (default: 16000).

Disks (`disks`) is used for Cloud platforms (Google Cloud Platforms, AWS, ...).

* `"atac.trim_adapter.disks"` : (optional) Disks for `trim_adapter` (default: "local-disk 100 HDD").
* `"atac.bowtie2.disks"` : (optional) Disks for `bowtie2` (default: "local-disk 100 HDD").
* `"atac.filter.disks"` : (optional) Disks for `filter` (default: "local-disk 100 HDD").
* `"atac.bam2ta.disks"` : (optional) Disks for `bam2ta` (default: "local-disk 100 HDD").
* `"atac.xcor.disks"` : (optional) Disks for `xcor` (default: "local-disk 100 HDD").
* `"atac.macs2_disks"` : (optional) Disks for `macs2` (default: "local-disk 100 HDD").

Walltime (`time`) settings (for SGE and SLURM only).

* `"atac.trim_adapter.time_hr"` : (optional) Walltime for `trim_adapter` (default: 24).
* `"atac.bowtie2.time_hr"` : (optional) Walltime for `bowtie2` (default: 48).
* `"atac.filter.time_hr"` : (optional) Walltime for `filter` (default: 24).
* `"atac.bam2ta.time_hr"` : (optional) Walltime for `bam2ta` (default: 6).
* `"atac.xcor.time_hr"` : (optional) Walltime for `xcor` (default: 6).
* `"atac.macs2_time_hr"` : (optional) Walltime for `macs2` (default: 24).
* `"atac.ataqc.time_hr"` : (optional) Walltime for `ATAQC` (default: 24).

