# Input JSON

An input JSON file includes all genomic data files, parameters and metadata for running pipelines. Our pipeline will use default values if they are not defined in an input JSON file. We provide a set of template JSON files: [minimum](../example_input_json/template.json) and [full](../example_input_json/template.full.json). We recommend to use a minimum template instead of full one. A full template includes all parameters of the pipeline with default values defined.

Please read through the following step-by-step instruction to compose a input JSON file.

>**IMPORTANT**: ALWAYS USE ABSOLUTE PATHS.

## Pipeline metadata

Parameter|Description
---------|-----------
`atac.title`| Title for experiment which will be shown in a final HTML report
`atac.description`| Description for experiment which will be shown in a final HTML report

## Pipeline parameters

Parameter|Default|Description
---------|-------|-----------
`atac.pipeline_type`| `atac` | `atac` for ATAC-seq or `dnase` for DNase-seq
`atac.align_only`| false | Peak calling and its downstream analyses will be disabled. Useful if you just want to map your FASTQs into filtered BAMs/TAG-ALIGNs and don't want to call peaks on them.
`atac.true_rep_only` | false | Disable pseudo replicate generation and all related analyses

## Reference genome

All reference genome specific reference files/parameters can be defined in a single TSV file `atac.genome_tsv`. However, you can also individally define each file/parameter instead of a TSV file. If both a TSV file and individual parameters are defined, then individual parameters will override those defined in a TSV file. For example, if you define both `atac.genome_tsv` and `atac.blacklist`, then `atac.blacklist` will override that is defined in `atac.genome_tsv`. This is useful when you want to use your own for a specific parameter while keeping all the other parameters same as original.

Parameter|Type|Description
---------|-------|-----------
`atac.genome_tsv`| File | Choose one of the TSV files listed below or build your own
`atac.genome_name`| String | Name of genome (e.g. hg38, hg19, ...)
`atac.ref_fa`| File | Reference FASTA file
`atac.ref_mito_fa`| File | Mito-only reference FASTA file
`atac.bowtie2_idx_tar`| File | Bowtie2 index TAR file (uncompressed) built from FASTA file
`atac.bowtie2_mito_idx_tar`| File | Mito-only Bowtie2 index TAR file (uncompressed) built from FASTA file
`atac.chrsz`| File | 2-col chromosome sizes file built from FASTA file with `faidx`
`atac.blacklist`| File | BED file. Peaks overlapping these regions will be filtered out
`atac.blacklist2`| File | Second blacklist. Two blacklist files (`atac.blacklist` and `atac.blacklist2`) will be merged.
`atac.gensz`| String | MACS2's genome sizes (hs for human, mm for mouse or sum of 2nd col in chrsz)
`atac.mito_chr_name`| String | Name of mitochondrial chromosome (e.g. chrM)
`atac.regex_bfilt_peak_chr_name`| String | Perl style reg-ex to keep peaks on selected chromosomes only matching with this pattern (default: `chr[\dXY]+`. This will keep chr1, chr2, ... chrX and chrY in `.bfilt.` peaks file. chrM is not included here)

Additional annotated genome data:

Parameter|Type|Description
---------|-------|-----------
`atac.tss` | File | TSS file
`atac.dnase` | File | Open chromatin region file
`atac.prom` | File | Promoter region file
`atac.enh` | File | Enhancer region file
`atac.reg2map` | File | File with cell type signals
`atac.reg2map_bed` | File | File of regions used to generate reg2map signals
`atac.roadmap_meta` | File | Roadmap metadata

We assume that users run pipeline with [Caper](https://github.com/ENCODE-DCC/caper/tree/master/caper). These TSVs work with Caper only since they have URLs instead of local paths or cloud bucket URIs. Caper will automatically download those URLs to a local temporary directory (`caper run ... --tmp-dir`).

We currently provide TSV files for 4 genomes as shown in the below table. You can [download/build](build_genome_database.md) them on your local computer. You can also [build a genome database for your own genome](build_genome_database.md).

Genome|URL
-|-
hg38|`https://storage.googleapis.com/encode-pipeline-genome-data/genome_tsv/v3/hg38.tsv`
mm10|`https://storage.googleapis.com/encode-pipeline-genome-data/genome_tsv/v3/mm10.tsv`
hg19|`https://storage.googleapis.com/encode-pipeline-genome-data/genome_tsv/v1/hg19_caper.tsv`
mm9|`https://storage.googleapis.com/encode-pipeline-genome-data/genome_tsv/v1/mm9_caper.tsv`

For DNAnexus CLI (AWS project):
Genome|DX URI
-|-
hg38|`dx://project-BKpvFg00VBPV975PgJ6Q03v6:pipeline-genome-data/genome_tsv/v3/hg38.dx.tsv`
mm10|`dx://project-BKpvFg00VBPV975PgJ6Q03v6:pipeline-genome-data/genome_tsv/v3/mm10.dx.tsv`

For DNAnexus CLI (Azure project): 
Genome|DX URI
-|-
hg38|`dx://project-F6K911Q9xyfgJ36JFzv03Z5J:pipeline-genome-data/genome_tsv/v3/hg38.dx_azure.tsv`
mm10|`dx://project-F6K911Q9xyfgJ36JFzv03Z5J:pipeline-genome-data/genome_tsv/v3/mm10.dx_azure.tsv`

For DNAnexus Web UI (AWS project): Choose one of the following TSV file on `https://platform.DNAnexus.com/projects/BKpvFg00VBPV975PgJ6Q03v6/data/pipeline-genome-data/genome_tsv/v3`.
Genome|File name
-|-
hg38|`hg38.dx.tsv`
mm10|`mm10.dx.tsv`

For DNAnexus Web UI (Azure project): Choose one of the following TSV file on `https://platform.DNAnexus.com/projects/F6K911Q9xyfgJ36JFzv03Z5J/data/pipeline-genome-data/genome_tsv/v3`.
Genome|File name
-|-
hg38|`hg38.dx_azure.tsv`
mm10|`mm10.dx_azure.tsv`

Additional information about each genome:

|Genome|Source|built from|
|-|-|-|
|hg38|ENCODE|[GRCh38_no_alt_analysis_set_GCA_000001405](https://www.encodeproject.org/files/GRCh38_no_alt_analysis_set_GCA_000001405.15/@@download/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta.gz)|
|mm10|ENCODE|[mm10_no_alt_analysis_set_ENCODE](https://www.encodeproject.org/files/mm10_no_alt_analysis_set_ENCODE/@@download/mm10_no_alt_analysis_set_ENCODE.fasta.gz)|
|hg19|UCSC|[GRCh37/hg19](http://hgdownload.cse.ucsc.edu/goldenpath/hg19/encodeDCC/referenceSequences/male.hg19.fa.gz)|
|mm9|UCSC|[mm9, NCBI Build 37](<http://hgdownload.cse.ucsc.edu/goldenPath/mm9/bigZips/mm9.2bit>)|

## How to download reference genome

1. Choose `GENOME` from `hg19`, `hg38`, `mm9` and `mm10` and specify a destination directory.
    ```bash
    $ bash genome/download_genome_data.sh [GENOME] [DESTINATION_DIR]
    ```
2. Find a TSV file on the destination directory and use it for `"atac.genome_tsv"` in your input JSON.

## Input genomic data

Choose endedness of your dataset first.

Parameter|Description
---------|-----------
`atac.paired_end`| Boolean to define endedness for ALL replicates. This will override per-replicate definition in `atac.paired_ends`
`atac.paired_ends`| Array of Boolean to define endedness for each replicate

Define `atac.paired_end` if all replicates in your dataset has the same endedness. You can also individually define endedness for each replicate. For example, rep1, rep2 are PE and rep3 is SE.

```javascript
{
    "atac.paired_ends" : [true, true, false]
}
```

Pipeline can start from any of the following data type (FASTQ, BAM, NODUP_BAM and TAG-ALIGN).

Parameter|Description
---------|-----------
`atac.fastqs_repX_R1`| Array of R1 FASTQ files for replicate X. These files will be merged into one FASTQ file for rep X.
`atac.fastqs_repX_R2`| Array of R2 FASTQ files for replicate X. These files will be merged into one FASTQ file for rep X. Do not define for single ended dataset. 
`atac.bams`| Array of BAM file for each replicate. (e.g. `["rep1.bam", "rep2.bam", ...]`)
`atac.nodup_bams`| Array of filtered/deduped BAM file for each replicate.
`atac.tas`| Array of TAG-ALIGN file for each replicate.

You can mix up different data types for individual replicate. For example, pipeline can start from FASTQs for rep1 and rep3, BAMs for rep2, NODUP_BAMs for rep4 and TAG-ALIGNs for rep5.

```javascript
{
    "atac.fastqs_rep1_R1" : ["rep1.fastq.gz"],
    "atac.fastqs_rep3_R1" : ["rep3.fastq.gz"],
    "atac.bams" : [null, "rep2.bam", null, null, null],
    "atac.nodup_bams" : [null, null, null, "rep4.nodup.bam", null],
    "atac.tas" : [null, null, null, null, "rep5.tagAlign.gz"]
}
```

## Adapter-trimming for FASTQs

If you choose to use auto-detection for adapters, then remove adapter arrays from input JSON. Otherwise define adapters for each FASTQ.

> **WARNING**: Individually defined adapters arrays should have the same dimension as FASTQs.

Parameter|Description
---------|-----------
`atac.adapter` | You can define an adapter sequence for ALL fastqs. If defined, this will override below adapter sequence definition for individual fastqs
`atac.adapters_repX_R1` | Array of adapter sequences for R1 FASTQs of replicate X
`atac.adapters_repX_R2` | Array of adapter sequences for R1 FASTQs of replicate X. Do not define it for singled-ended dataset

## Optional adapter-trimming parameters

Parameter|Default|Description
---------|-------|-----------
`atac.auto_detect_adapter` | false | You can use auto-detection for adapters. List of adapters can be detected: AGATCGGAAGAGC (Illumina), CTGTCTCTTATA (Nextera) and TGGAATTCTCGG (smallRNA)
`atac.cutadapt_param` | `-e 0.1 -m 5` | cutadapt (trim_adapter) parameters (default: min_trim_len=5, err_rate=0.1)

## Optional mapping parameters

Parameter|Type | Default|Description
---------|----|---|-----------
`atac.multimapping` | Int | 4 | Multimapping reads

## Optional filtering parameters

Parameter|Default|Description
---------|-------|-----------
`atac.mapq_thresh` | 30 | Threshold for mapped reads quality (samtools view -q). If not defined, automatically determined according to aligner.
`atac.dup_marker` | `picard` | Choose a dup marker between `picard` and `sambamba`. `picard` is recommended, use `sambamba` only when picard fails.
`atac.no_dup_removal` | false | Skip dup removal in a BAM filtering stage.

## Optional subsampling parameters

Parameter|Default|Description
---------|-------|-----------
`atac.subsample_reads` | 0 | Subsample reads (0: no subsampling). For PE dataset, this is not a number of read pairs but number of reads. Subsampled reads will be used for all downsteam analyses including peak-calling
`atac.xcor_subsample_reads` | 15000000 | Subsample reads for cross-corr. analysis only (0: no subsampling). Subsampled reads will be used for cross-corr. analysis only

## Optional peak-calling parameters

Parameter|Default|Description
---------|-------|-----------
`atac.cap_num_peak` | 500000 | Cap number of peaks called from a peak-caller (MACS2)
`atac.pval_thresh` | 0.01 | P-value threshold for MACS2 (macs2 callpeak -p).
`atac.smooth_win` | 150 | Size of smoothing window for MACS2 (macs2 callpeak --shift [-smooth_win/2] --extsize [smooth_win]).
`atac.enable_idr` | true | Enable IDR (irreproducible discovery rate)
`atac.idr_thresh` | 0.05 | Threshold for IDR

## Optional pipeline flags

Parameter|Default|Description
---------|-------|-----------
`atac.enable_xcor` | false | Enable cross-correlation analysis
`atac.enable_count_signal_track` | false | Enable count signal track generation
`atac.enable_preseq` | false | Enable preseq, which performs a yield prediction for reads
`atac.enable_jsd` | true | Enable deeptools fingerprint (JS distance)
`atac.enable_gc_bias` | true | Enable GC bias computation
`atac.enable_tss_enrich` | true | Enable TSS enrichment computation
`atac.enable_annot_enrich` | true | Enable Annotated region enrichment computation
`atac.enable_compare_to_roadmap` | false | Enable comparing signals to epigenome roadmap

## Optional parameter for TSS enrichment

Our pipeline automatically estimates read length from FASTQs, but `atac.read_len` will override those estimated ones. You need to define `atac.read_len` if you start from BAMs and want to get a TSS enrichment plot.

Parameter|Type | Description
---------|-----|-----------
`atac.read_len` | `Array[Int]` | Read length for each replicate. 

## Other optional parameters

Parameter|Default|Description
---------|-------|-----------
`atac.filter_chrs` | `["chrM", "MT"]` | Array of chromosome names to be filtered out from a final (filtered/nodup) BAM. Mitochondrial chromosomes are filtered out by default.

> **WARNING**: If your custom genome's mitochondrial chromosome name is different from `chrM` or `MT`, then define it correctly here. This parameter has nothing to do with a mito-chromosome name parameter `atac.mito_chr_name`. Changing `atac.mito_chr_name` does not affect this parameter.

## Resource parameters

> **WARNING**: It is recommened not to change the following parameters unless you get resource-related errors for a certain task and you want to increase resources for such task. The following parameters are provided for users who want to run our pipeline with Caper's `local` on HPCs and 2).

Resources defined here are PER REPLICATE. Therefore, total number of cores will be approximately `atac.align_cpu` x `NUMBER_OF_REPLICATES` because `align` is a bottlenecking task of the pipeline. Use this total number of cores if you manually `qsub` or `sbatch` your job (using local mode of Caper). `disk_factor` is used for Google Cloud and DNAnexus only.

For example, if sum of your FASTQs are 20GB then 4GB (base) + `atac.align_mem_factor` x 20GB = 5GB will be used for `align` task's instance memory.

If sum of your TAG-ALIGN BEDs (intermediate outputs) are 5GB then 4GB (base) + `atac.macs2_signal_track_mem_factor` x 5GB = 34GB will be used for `macs2_signal_track` task's instance memory.

Base memory/disk is 4GB/20GB for most tasks.

Parameter|Default|Description
---------|-------|-----------
`atac.align_cpu` | 6 |
`atac.align_mem_factor` | 0.15 | Multiplied to size of FASTQs to determine required memory
`atac.align_time_hr` | 48 | Walltime (HPCs only)
`atac.align_disk_factor` | 8.0 | Multiplied to size of FASTQs to determine required disk

Parameter|Default|Description
---------|-------|-----------
`atac.filter_cpu` | 4 |
`atac.filter_mem_factor` | 0.4 | Multiplied to size of BAM to determine required memory
`atac.filter_time_hr` | 24 | Walltime (HPCs only)
`atac.filter_disk_factor` | 4.0 | Multiplied to size of BAM to determine required disk

Parameter|Default|Description
---------|-------|-----------
`atac.bam2ta_cpu` | 2 |
`atac.bam2ta_mem_factor` | 0.3 | Multiplied to size of filtered BAM to determine required memory
`atac.bam2ta_time_hr` | 6 | Walltime (HPCs only)
`atac.bam2ta_disk_factor` | 4.0 | Multiplied to size of filtered BAM to determine required disk

Parameter|Default|Description
---------|-------|-----------
`atac.spr_mem_factor` | 4.5 | Multiplied to size of filtered BAM to determine required memory
`atac.spr_disk_factor` | 6.0 | Multiplied to size of filtered BAM to determine required disk

Parameter|Default|Description
---------|-------|-----------
`atac.jsd_cpu` | 4 |
`atac.jsd_mem_factor` | 0.1 | Multiplied to size of filtered BAM to determine required memory
`atac.jsd_time_hr` | 6 | Walltime (HPCs only)
`atac.jsd_disk_factor` | 2.0 | Multiplied to size of filtered BAM to determine required disk

Parameter|Default|Description
---------|-------|-----------
`atac.xcor_cpu` | 2 |
`atac.xcor_mem_factor` | 1.0 | Multiplied to size of TAG-ALIGN BED to determine required memory
`atac.xcor_time_hr` | 6 | Walltime (HPCs only)
`atac.xcor_disk_factor` | 4.5 | Multiplied to size of TAG-ALIGN BED to determine required disk

Parameter|Default|Description
---------|-------|-----------
`atac.call_peak_cpu` | 2 | MACS2 is single-threaded. More than 2 is not required.
`atac.call_peak_mem_factor` | 2.0 | Multiplied to size of TAG-ALIGN BED to determine required memory
`atac.call_peak_time_hr` | 24 | Walltime (HPCs only)
`atac.call_peak_disk_factor` | 15.0 | Multiplied to size of TAG-ALIGN BED to determine required disk

Parameter|Default|Description
---------|-------|-----------
`atac.macs2_signal_track_mem_factor` | 6.0 | Multiplied to size of TAG-ALIGN BED to determine required memory
`atac.macs2_signal_track_time_hr` | 24 | Walltime (HPCs only)
`atac.macs2_signal_track_disk_factor` | 40.0 | Multiplied to size of TAG-ALIGN BED to determine required disk

Parameter|Default|Description
---------|-------|-----------
`atac.preseq_mem_factor` | 0.5 | Multiplied to size of BAM to determine required memory
`atac.preseq_disk_factor` | 5.0 | Multiplied to size of BAM to determine required disk

If your system/cluster does not allow large memory allocation for Java applications, check the following resource parameters to manually define Java memory. It is **NOT RECOMMENDED** for most users to change these parameters since pipeline automatically takes 90% of task's memory for Java apps.

There are special parameters to control maximum Java heap memory (e.g. `java -Xmx4G`) for Java applications (e.g. Picard tools). They are strings including size units. Such string will be directly appended to Java's parameter `-Xmx`. If these parameters are not defined then pipeline uses 90% of each task's memory.

Parameter|Default
---------|-------
`atac.filter_picard_java_heap` | 90% of memory for `atac.filter` (dynamic)
`atac.preseq_picard_java_heap` | 90% of memory for `atac.preseq` (dynamic)
`atac.fraglen_stat_picard_java_heap` | 90% of memory for `atac.fraglen_stat_pe` (8GB)
`atac.gc_bias_picard_java_heap` | 90% of memory for `atac.gc_bias` (8GB)
