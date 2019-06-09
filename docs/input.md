# Input JSON

An input JSON file includes all genomic data files, parameters and metadata for running pipelines. Our pipeline will use default values if they are not defined in an input JSON file. We provide a set of template JSON files: [minimum](../examples/template.json) and [full](../examples/template.full.json). We recommend to use a minimum template instead of full one. A full template includes all parameters of the pipeline with default values defined.

Please read through the following step-by-step instruction to compose a input JSON file.

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
`atac.ref_fa`| File | Reference FASTA file
`atac.bowtie2_idx_tar`| File | Bowtie2 index TAR file (uncompressed) built from FASTA file
`atac.chrsz`| File | 2-col chromosome sizes file built from FASTA file with `faidx`
`atac.blacklist`| File | 3-col BED file. Peaks overlapping these regions will be filtered out
`atac.gensz`| String | MACS2's genome sizes (hs for human, mm for mouse or sum of 2nd col in chrsz)

Additional reference genome data for ATAqC:

Parameter|Type|Description
---------|-------|-----------
`atac.tss` | File | TSS file
`atac.dnase` | File | Open chromatin region file
`atac.prom` | File | Promoter region file
`atac.enh` | File | Enhancer region file
`atac.reg2map` | File | File with cell type signals
`atac.reg2map_bed` | File | File of regions used to generate reg2map signals
`atac.roadmap_meta` | File | Roadmap metadata

We currently provide TSV files for 4 genomes as shown in the below table. `GENOME` should be `hg38`, `mm10`, `hg19` or `mm9`. You can [download/build](build_genome_database.md) it on your local computer. You can also [build a genome database for your own genome](build_genome_database.md).

Platform|Path/URI
-|-
Google Cloud Platform|`gs://encode-pipeline-genome-data/[GENOME]_google.tsv`
Stanford Sherlock|`/home/groups/cherry/encode/pipeline_genome_data/[GENOME]_sherlock.tsv`
Stanford SCG|`/reference/ENCODE/pipeline_genome_data/[GENOME]_scg.tsv`
Local/SLURM/SGE|You need to [build](build_genome_database.md) or [download]() a genome database]. 
DNAnexus (CLI)|`dx://project-BKpvFg00VBPV975PgJ6Q03v6:pipeline-genome-data/[GENOME]_dx.tsv`
DNAnexus (CLI, Azure)|`dx://project-F6K911Q9xyfgJ36JFzv03Z5J:pipeline-genome-data/[GENOME]_dx_azure.tsv`
DNAnexus (Web)|Choose `[GENOME]_dx.tsv` from [here](https://platform.DNAnexus.com/projects/BKpvFg00VBPV975PgJ6Q03v6/data/pipeline-genome-data)
DNAnexus (Web, Azure)|Choose `[GENOME]_dx.tsv` from [here](https://platform.DNAnexus.com/projects/F6K911Q9xyfgJ36JFzv03Z5J/data/pipeline-genome-data)

Additional information about each genome:

|Genome|Source|built from|
|-|-|-|
|hg38|ENCODE|[GRCh38_no_alt_analysis_set_GCA_000001405](https://www.encodeproject.org/files/GRCh38_no_alt_analysis_set_GCA_000001405.15/@@download/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta.gz)|
|mm10|ENCODE|[mm10_no_alt_analysis_set_ENCODE](https://www.encodeproject.org/files/mm10_no_alt_analysis_set_ENCODE/@@download/mm10_no_alt_analysis_set_ENCODE.fasta.gz)|
|hg19|UCSC|[GRCh37/hg19](http://hgdownload.cse.ucsc.edu/goldenpath/hg19/encodeDCC/referenceSequences/male.hg19.fa.gz)|
|mm9|UCSC|[mm9, NCBI Build 37](<http://hgdownload.cse.ucsc.edu/goldenPath/mm9/bigZips/mm9.2bit>)|

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

You can mix up different data types for individual replicate/control replicate. For example, pipeline can start from FASTQs for rep1 and rep3, BAMs for rep2, NODUP_BAMs for rep4 and TAG-ALIGNs for rep5. You can define similarly for control replicates.

```javascript
{
    "atac.fastqs_rep1_R1" : ["rep1.fastq.gz"],
    "atac.fastqs_rep3_R1" : ["rep3.fastq.gz"],
    "atac.bams" : [null, "rep2.bam", null, null, null],
    "atac.nodup_bams" : [null, "rep2.bam", null, "rep4.nodup.bam", null],
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

Parameter|Default|Description
---------|-------|-----------
`atac.multimapping` | 4 | Multimapping reads
`atac.bowtie2_param_pe` | `-X2000 --mm --local` | bowtie2 parameters for each read-endedness (paired-end). See bowtie2 --help for details.
`atac.bowtie2_param_se` | `--local` | bowtie2 parameters for each read-endedness (single-ended). See bowtie2 --help for details.

## Optional filtering parameters

Parameter|Default|Description
---------|-------|-----------
`atac.mapq_thresh` | 30 | Threshold for mapped reads quality (samtools view -q)
`atac.dup_marker` | `picard` | Choose a dup marker between `picard` and `sambamba`. `picard` is recommended, use `sambamba` only when picard fails.
`atac.no_dup_removal` | false | Skip dup removal in a BAM filtering stage.

## Optional subsampling parameters

Parameter|Default|Description
---------|-------|-----------
`atac.subsample_reads` | 0 | Subsample reads (0: no subsampling). Subsampled reads will be used for all downsteam analyses including peak-calling
`atac.xcor_subsample_reads` | 15000000 | Subsample reads for cross-corr. analysis only (0: no subsampling). Subsampled reads will be used for cross-corr. analysis only

## Optional peak-calling parameters

Parameter|Default|Description
---------|-------|-----------
`atac.cap_num_peak` | 500000 | Cap number of peaks called from a peak-caller (MACS2)
`atac.pval_thresh` | 0.01 | P-value threshold for MACS2 (macs2 callpeak -p)
`atac.enable_idr` | true | Enable IDR (irreproducible discovery rate)
`atac.idr_thresh` | 0.05 | Threshold for IDR

## Optional pipeline flags

Parameter|Default|Description
---------|-------|-----------
`atac.enable_xcor` | true | Enable cross-correlation analysis
`atac.enable_count_signal_track` | false | Enable count signal track generation
`atac.keep_irregular_chr_in_bfilt_peak` | false | Keep irregular chromosome names. Use this for custom genomes without canonical chromosome names (chr1, chrX, ...)
`atac.disable_ataqc` | false | Disable ATAqC (including all annotation-based analyses in it)

## Other optional parameters

Parameter|Default|Description
---------|-------|-----------
`atac.mito_chr_name` | `chrM` | Name of mito chromosome. THIS IS NOT A REG-EX! you can define only one chromosome name for mito.
`atac.regex_filter_reads` | `chrM` | Regular expression to filter out reads with given chromosome name (1st column of BED/TAG-ALIGN). Any read with chr name that matches with this reg-ex pattern will be removed from outputs If your have changed the above parameter `atac.mito_chr_name` and still want to filter out mito reads then make sure that `atac.mito_chr_name` and `atac.regex_filter_reads` are the same

## Resource parameters

> **WARNING**: It is recommened not to change the following parameters unless you get resource-related errors for a certain task and you want to increase resources for such task. The following parameters are provided for users who want to run our pipeline with Caper's `local` on HPCs and 2).

Resources defined here are PER REPLICATE. Therefore, total number of cores will be approximately `atac.bowtie2_cpu` x `NUMBER_OF_REPLICATES because bowtie2 is a bottlenecking task of the pipeline. Use this total number of cores if you manually `qsub` or `sbatch` your job (using local mode of Caper). `disks` is used for Google Cloud and DNAnexus only.

Parameter|Default
---------|-------
`atac.trim_adapter_cpu` | 2
`atac.trim_adapter_mem_mb` | 12000
`atac.trim_adapter_time_hr` | 24
`atac.trim_adapter_disks` | `local-disk 100 HDD`

Parameter|Default
---------|-------
`atac.bowtie2_cpu` | 4
`atac.bowtie2_mem_mb` | 20000
`atac.bowtie2_time_hr` | 48
`atac.bowtie2_disks` | `local-disk 100 HDD`

Parameter|Default
---------|-------
`atac.filter_cpu` | 2
`atac.filter_mem_mb` | 20000
`atac.filter_time_hr` | 24
`atac.filter_disks` | `local-disk 100 HDD`

Parameter|Default
---------|-------
`atac.bam2ta_cpu` | 2
`atac.bam2ta_mem_mb` | 10000
`atac.bam2ta_time_hr` | 6
`atac.bam2ta_disks` | `local-disk 100 HDD`

Parameter|Default
---------|-------
`atac.spr_mem_mb` | 16000

Parameter|Default
---------|-------
`atac.fingerprint_cpu` | 2
`atac.fingerprint_mem_mb` | 12000
`atac.fingerprint_time_hr` | 6
`atac.fingerprint_disks` | `local-disk 100 HDD`

Parameter|Default
---------|-------
`atac.xcor_cpu` | 2
`atac.xcor_mem_mb` | 16000
`atac.xcor_time_hr` | 6
`atac.xcor_disks` | `local-disk 100 HDD`

Parameter|Default
---------|-------
`atac.macs2_mem_mb` | 16000
`atac.macs2_time_hr` | 24
`atac.macs2_disks` | `local-disk 100 HDD`

Make sure that ataqc_mem_java_mb < ataqc_mem_mb

Parameter|Default
---------|-------
`atac.ataqc_mem_mb` | 16000
`atac.ataqc_mem_java_mb` | 15000
`atac.ataqc_time_hr` | 24
`atac.ataqc_disks` | `local-disk 200 HDD`

