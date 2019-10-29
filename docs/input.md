# Input JSON

An input JSON file includes all genomic data files, parameters and metadata for running pipelines. Our pipeline will use default values if they are not defined in an input JSON file. We provide a set of template JSON files: [minimum](../example_input_json/template.json) and [full](../example_input_json/template.full.json). We recommend to use a minimum template instead of full one. A full template includes all parameters of the pipeline with default values defined.

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
`atac.genome_name`| String | Name of genome (e.g. hg38, hg19, ...)
`atac.ref_fa`| File | Reference FASTA file
`atac.ref_mito_fa`| File | Mito-only reference FASTA file
`atac.bowtie2_idx_tar`| File | Bowtie2 index TAR file (uncompressed) built from FASTA file
`atac.bowtie2_mito_idx_tar`| File | Mito-only Bowtie2 index TAR file (uncompressed) built from FASTA file
`atac.custom_aligner_idx_tar` | File | Index TAR file (uncompressed) for your own aligner. See details about [how to use a custom aligner](#how-to-use-a-custom-aligner)
`atac.custom_aligner_mito_idx_tar` | File | Mito-only index TAR file (uncompressed) for your own aligner. See details about [how to use a custom aligner](#how-to-use-a-custom-aligner)
`atac.chrsz`| File | 2-col chromosome sizes file built from FASTA file with `faidx`
`atac.blacklist`| File | 3-col BED file. Peaks overlapping these regions will be filtered out
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

We currently provide TSV files for 4 genomes as shown in the below table. `GENOME` should be `hg38`, `mm10`, `hg19` or `mm9`. You can [download/build](build_genome_database.md) it on your local computer. You can also [build a genome database for your own genome](build_genome_database.md).

Platform|Path/URI
-|-
Google Cloud Platform|`gs://encode-pipeline-genome-data/genome_tsv/v1/[GENOME]_gcp.tsv`
Stanford Sherlock|`/home/groups/cherry/encode/pipeline_genome_data/genome_tsv/v1/[GENOME]_sherlock.tsv`
Stanford SCG|`/reference/ENCODE/pipeline_genome_data/genome_tsv/v1/[GENOME]_scg.tsv`
Local/SLURM/SGE|You need to [build](build_genome_database.md) or [download]() a genome database]. 
DNAnexus (CLI)|`dx://project-BKpvFg00VBPV975PgJ6Q03v6:pipeline-genome-data/genome_tsv/v1/[GENOME]_dx.tsv`
DNAnexus (CLI, Azure)|`dx://project-F6K911Q9xyfgJ36JFzv03Z5J:pipeline-genome-data/genome_tsv/v1/[GENOME]_dx_azure.tsv`
DNAnexus (Web)|Choose `[GENOME]_dx.tsv` from [here](https://platform.DNAnexus.com/projects/BKpvFg00VBPV975PgJ6Q03v6/data/pipeline-genome-data/genome_tsv/v1)
DNAnexus (Web, Azure)|Choose `[GENOME]_dx.tsv` from [here](https://platform.DNAnexus.com/projects/F6K911Q9xyfgJ36JFzv03Z5J/data/pipeline-genome-data/genome_tsv/v1)

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
`atac.custom_align_py` | File | | Python script for your custom aligner. See details about [how to use a custom aligner](#how-to-use-a-custom-aligner)

## Optional filtering parameters

Parameter|Default|Description
---------|-------|-----------
`atac.mapq_thresh` | 30 | Threshold for mapped reads quality (samtools view -q). If not defined, automatically determined according to aligner.
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
`atac.custom_call_peak_py` | File | Python script for your custom peak caller. See details about [how to use a custom peak caller](#how-to-use-a-peak-caller)

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

Resources defined here are PER REPLICATE. Therefore, total number of cores will be approximately `atac.align_cpu` x `NUMBER_OF_REPLICATES because `align` is a bottlenecking task of the pipeline. Use this total number of cores if you manually `qsub` or `sbatch` your job (using local mode of Caper). `disks` is used for Google Cloud and DNAnexus only.

Parameter|Default
---------|-------
`atac.align_cpu` | 4
`atac.align_mem_mb` | 20000
`atac.align_time_hr` | 48
`atac.align_disks` | `local-disk 400 HDD`

Parameter|Default
---------|-------
`atac.filter_cpu` | 2
`atac.filter_mem_mb` | 20000
`atac.filter_time_hr` | 24
`atac.filter_disks` | `local-disk 400 HDD`

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
`atac.jsd_cpu` | 2
`atac.jsd_mem_mb` | 12000
`atac.jsd_time_hr` | 6
`atac.jsd_disks` | `local-disk 200 HDD`

Parameter|Default
---------|-------
`atac.xcor_cpu` | 2
`atac.xcor_mem_mb` | 16000
`atac.xcor_time_hr` | 6
`atac.xcor_disks` | `local-disk 100 HDD`

Parameter|Default
---------|-------
`atac.call_peak_cpu` | 1
`atac.call_peak_mem_mb` | 16000
`atac.call_peak_time_hr` | 24
`atac.call_peak_disks` | `local-disk 200 HDD`

Parameter|Default
---------|-------
`atac.macs2_signal_track_mem_mb` | 16000
`atac.macs2_signal_track_time_hr` | 24
`atac.macs2_signal_track_disks` | `local-disk 200 HDD`

Parameter|Default
---------|-------
`atac.preseq_mem_mb` | 16000

> **IMPORTANT**: If you see memory Java errors, check the following resource parameters.

There are special parameters to control maximum Java heap memory (e.g. `java -Xmx4G`) for Picard tools. They are strings including size units. Such string will be directly appended to Java's parameter `-Xmx`.

Parameter|Default
---------|-------
`atac.filter_picard_java_heap` | `4G`
`atac.preseq_picard_java_heap` | `6G`
`atac.fraglen_stat_picard_java_heap` | `6G`
`atac.gc_bias_picard_java_heap` | `6G`


## How to use a custom aligner

ENCODE ATAC-Seq pipeline currently supports `bowtie2` only. In order to use your own aligner you need to define the following parameters first. You can define `custom_aligner_idx_tar` either in your input JSON file or in your genome TSV file. Such index TAR file should be an uncompressed TAR file without any directory structured.

Parameter|Type|Description
---------|-------|-----------
`atac.custom_aligner_idx_tar` | File | Index TAR file (uncompressed) for your own aligner
`atac.custom_align_py` | File | Python script for your custom aligner

Here is a template for `custom_align.py`:

```python
#!/usr/bin/env python

import os
import argparse

def parse_arguments():
    parser = argparse.ArgumentParser(prog='ENCODE template aligner')
    parser.add_argument('index_prefix_or_tar', type=str,
                        help='Path for prefix (or a tarball .tar) \
                            for reference aligner index. \
                            Tar ball must be packed without compression \
                            and directory by using command line \
                            "tar cvf [TAR] [TAR_PREFIX].*')
    parser.add_argument('fastqs', nargs='+', type=str,
                        help='List of FASTQs (R1 and R2). \
                            FASTQs must be compressed with gzip (with .gz).')
    parser.add_argument('--paired-end', action="store_true",
                        help='Paired-end FASTQs.')
    parser.add_argument('--multimapping', default=4, type=int,
                        help='Multimapping reads')
    parser.add_argument('--nth', type=int, default=1,
                        help='Number of threads to parallelize.')
    parser.add_argument('--out-dir', default='', type=str,
                            help='Output directory.')
    args = parser.parse_args()

    # check if fastqs have correct dimension
    if args.paired_end and len(args.fastqs)!=2:
        raise argparse.ArgumentTypeError('Need 2 fastqs for paired end.')
    if not args.paired_end and len(args.fastqs)!=1:
        raise argparse.ArgumentTypeError('Need 1 fastq for single end.')

    return args

def align(fastq_R1, fastq_R2, ref_index_prefix, multimapping, nth, out_dir):
    basename = os.path.basename(os.path.splitext(fastq_R1)[0])    
    prefix = os.path.join(out_dir, basename)
    bam = '{}.bam'.format(prefix)

    # map your fastqs somehow
    os.system('touch {}'.format(bam))

    return bam

def main():
    # read params
    args = parse_arguments()
   
    # unpack index somehow on CWD
    os.system('tar xvf {}'.format(args.index_prefix_or_tar))

    bam = align(args.fastqs[0],
                args.fastqs[1] if args.paired_end else None,
                args.index_prefix_or_tar,
                args.multimapping,
                args.nth,
                args.out_dir)

if __name__=='__main__':
    main()

```

> **IMPORTANT**: Your custom python script should generate ONLY one `*.bam` file. For example, if there are two `.bam` files then pipeline will just pick the first one in an alphatical order.

## How to use a custom peak caller

Parameter|Type|Default|Description
---------|-------|-----------
`atac.peak_type` | String | `narrowPeak` | Only ENCODE peak types are supported: `narrowPeak`, `broadPeak` and `gappedPeak`
`atac.custom_call_peak_py` | File | | Python script for your custom peak caller

The file extension of your output peak file must be consitent with the `peak_type` you chose. For example, if you have chosen `narrowPeak` as `peak_type` then your output peak file should be `*.narrowPeak.gz`.

Here is a template for `custom_call_peak.py`:

```python
#!/usr/bin/env python

import os
import argparse

def parse_arguments():
    parser = argparse.ArgumentParser(prog='ENCODE template call_peak')
    parser.add_argument('ta', type=str,
                        help='Path for TAGALIGN file')
    parser.add_argument('--fraglen', type=int, required=True,
                        help='Fragment length.')
    parser.add_argument('--shift', type=int, default=0,
                        help='macs2 callpeak --shift.')
    parser.add_argument('--chrsz', type=str,
                        help='2-col chromosome sizes file.')
    parser.add_argument('--gensz', type=str,
                        help='Genome size (sum of entries in 2nd column of \
                            chr. sizes file, or hs for human, ms for mouse).')
    parser.add_argument('--pval-thresh', default=0.01, type=float,
                        help='P-Value threshold.')
    parser.add_argument('--cap-num-peak', default=500000, type=int,
                        help='Capping number of peaks by taking top N peaks.')
    parser.add_argument('--out-dir', default='', type=str,
                        help='Output directory.')
    args = parser.parse_args()
    return args

def call_peak(ta, chrsz, gensz, pval_thresh, shift, fraglen, cap_num_peak, out_dir):
    basename_ta = os.path.basename(os.path.splitext(ta)[0])
    basename_prefix = basename_ta

    prefix = os.path.join(out_dir, basename_prefix)
    npeak = '{}.narrowPeak.gz'.format(prefix)

    os.system('touch {}'.format(npeak))

    return npeak

def main():
    # read params
    args = parse_arguments()

    npeak = call_peak(
        args.ta, args.chrsz, args.gensz, args.pval_thresh,
        args.shift, args.fraglen, args.cap_num_peak, args.out_dir)

if __name__=='__main__':
    main()
```

> **IMPORTANT**: Your custom python script should generate ONLY one `*.*Peak.gz` file. For example, if there are `*.narrowPeak.gz` and `*.broadPeak.gz` files pipeline will just pick the first one in an alphatical order.
