# Input JSON

An input JSON file includes all genomic data files, parameters and metadata for running pipelines. Our pipeline will use default values if they are not defined in an input JSON file. We provide a set of template JSON files: [minimum](../example_input_json/template.json) and [full](../example_input_json/template.full.json). We recommend to use a minimum template instead of full one. A full template includes all parameters of the pipeline with default values defined.

>**IMPORTANT**: ALWAYS USE ABSOLUTE PATHS.

# Checklist

Mandatory parameters:

1) Pipeline type
    * `atac.pipeline_type`: `atac` for ATAC-seq or `dnase` for DNase-seq.

2) Experiment title/description
    * `atac.title`: experiment title for a final HTML report.
    * `atac.description`: experiment description for a final HTML report.

3) Read endedness
    * `atac.paired_end`: `true` if ALL replicates are paired-ended.
    * (Optional) `atac.paired_ends`: For samples with mixed read ends, you can define read endedness for each biological replicate (e.g. `[true, false]` means paired-ended biorep-1 and single-ended biorep-2).

4) Reference genome
    * `atac.genome_tsv`: Choose one from the following genome TSVs. `v3` was a standard for >=ENCODE4 and <=v2.0.2. It's updated to `v4` for >=v2.1.0 (on 12/01/2021). See the version history of genome TSV files [here](https://storage.googleapis.com/encode-pipeline-genome-data/genome_tsv/history.txt).
        Genome|URL
        -|-
        hg38|`https://storage.googleapis.com/encode-pipeline-genome-data/genome_tsv/v4/hg38.tsv`
        mm10|`https://storage.googleapis.com/encode-pipeline-genome-data/genome_tsv/v4/mm10.tsv`

        For DNAnexus CLI (AWS project):
        Genome|DX URI
        -|-
        hg38|`dx://project-BKpvFg00VBPV975PgJ6Q03v6:pipeline-genome-data/genome_tsv/v4/hg38.dx.tsv`
        mm10|`dx://project-BKpvFg00VBPV975PgJ6Q03v6:pipeline-genome-data/genome_tsv/v4/mm10.dx.tsv`

        For DNAnexus CLI (Azure project): 
        Genome|DX URI
        -|-
        hg38|`dx://project-F6K911Q9xyfgJ36JFzv03Z5J:pipeline-genome-data/genome_tsv/v4/hg38.dx_azure.tsv`
        mm10|`dx://project-F6K911Q9xyfgJ36JFzv03Z5J:pipeline-genome-data/genome_tsv/v4/mm10.dx_azure.tsv`

        For DNAnexus Web UI (AWS project): Choose one of the following TSV file on `https://platform.DNAnexus.com/projects/BKpvFg00VBPV975PgJ6Q03v6/data/pipeline-genome-data/genome_tsv/v4`.
        Genome|File name
        -|-
        hg38|`hg38.dx.tsv`
        mm10|`mm10.dx.tsv`

        For DNAnexus Web UI (Azure project): Choose one of the following TSV file on `https://platform.DNAnexus.com/projects/F6K911Q9xyfgJ36JFzv03Z5J/data/pipeline-genome-data/genome_tsv/v4`.
        Genome|File name
        -|-
        hg38|`hg38.dx_azure.tsv`
        mm10|`mm10.dx_azure.tsv`

    * To build a new TSV file from use your own FASTA (`.fa` and `.2bit`) see [this](build_genome_database.md).
    * You can also define genome specific parameters (defined in a genome TSV file) in an input JSON file. Parameters defined in an input JSON file will override those defined in a genome TSV file. For example, you can simply replace a blacklist file while keeping all other parameters in a genome TSV file for `hg38`.
        ```javascript
        {
            "atac.genome_tsv" : "/somewhere/hg38.tsv",
            "atac.blacklist": "/new/genome/data/new_blacklist.bed.gz"
        }
        ```

5) [Input files](#input-files) and [adapters](#adapters)
    * See [this](#input-files) for how to define FASTQ/BAM/TAG-ALIGNs for your sample.
    * See [this](#adapters) for how to define adapters to be trimmed.

6) Important parameters
    * `atac.auto_detect_adapter`: Automatically detect [some adapters](#adapters) and trim them in FASTQs.
    * `atac.multimapping`: Set it as 0 if you don't have multimapping reads. It's 4 by default.
    * `atac.read_len`: Array of Integers. Read length for each bio replicate. If you start from FASTQs simply forget about this parameter. Otherwise you want to get a TSS enrichment plot, then you must define it (e.g. `[75, 50]` means 75 for rep1 and 50 for rep2).
    * `atac.pval_thresh`: P-value threshold for MACS2 (macs2 callpeak -p).
    * `atac.smooth_win`: Size of smoothing window for MACS2 (macs2 callpeak --shift [-smooth_win/2] --extsize [smooth_win]).

7) [Resources](#resources)
    * It is recommened not to change the following parameters unless you get resource-related errors for a certain task and you want to increase resources for such task.

Optional parameters:

8) Useful parameters
    * `atac.subsample_reads`: Subsample reads. For PE dataset, this is not a number of read pairs but number of reads. This will affect all downsteam analyses including peak-calling. It's 0 by default, which means no subsampling.

9) Flags
    * `atac.align_only`: Peak calling and its downstream analyses will be disabled. Useful if you just want to align your FASTQs into filtered BAMs/TAG-ALIGNs and don't want to call peaks on them.
    * `atac.true_rep_only`: Disable pseudo replicate generation and all related analyses


## Input files

> **IMPORTANT**: Our pipeline considers a replicate (`rep`) as a biological replicate. You can still define technical replicates for each bio replicate. Tech replicates will be merged together to make a single FASTQ for each bio replicate.

> **IMPORTANT**: Our pipeline supports up to 10 bio replicates.

> **IMPORTANT**: Our pipeline has cross-validation analyses (IDR/overlap) comparing every pair of all replicates. Number of tasks for such analyses will be like <sub>n</sub>C<sub>2</sub>. This number will be 45 for 10 bio replicates. It's recommended to keep number of replicates <= 4.

Pipeline can start from any of the following data types (FASTQ, BAM, NODUP_BAM and TAG-ALIGN). 

1) Starting from FASTQs
    * Technical replicates for each bio-rep will be **MERGED** in the very early stage of the pipeline. Each read end R1 and R2 have separate arrays `atac.fastqs_repX_R1` and `atac.fastqs_repX_R2`. Do not define R2 array for single-ended replicates.
    * Example of 3 paired-ended biological replicates and 2 technical replicates for each bio rep. Two technical replicates `BIOREPX_TECHREP1.R1.fq.gz` and `BIOREPX_TECHREP2.R1.fq.gz` for each bio replicate will be merged.

        ```javascript
        {
            "atac.paired_end" : true,
            "atac.fastqs_rep1_R1" : ["BIOREP1_TECHREP1.R1.fq.gz", "BIOREP1_TECHREP2.R1.fq.gz"],
            "atac.fastqs_rep1_R2" : ["BIOREP1_TECHREP1.R2.fq.gz", "BIOREP1_TECHREP2.R2.fq.gz"],
            "atac.fastqs_rep2_R1" : ["BIOREP2_TECHREP1.R1.fq.gz", "BIOREP2_TECHREP2.R1.fq.gz"],
            "atac.fastqs_rep2_R2" : ["BIOREP2_TECHREP1.R2.fq.gz", "BIOREP2_TECHREP2.R2.fq.gz"],
            "atac.fastqs_rep3_R1" : ["BIOREP3_TECHREP1.R1.fq.gz", "BIOREP3_TECHREP2.R1.fq.gz"],
            "atac.fastqs_rep3_R2" : ["BIOREP3_TECHREP1.R2.fq.gz", "BIOREP3_TECHREP2.R2.fq.gz"]
        }
        ```

2) Starting from BAMs
    * Define a BAM for each replicate. Our pipeline does not determine read endedness from a BAM file. You need to explicitly define read endedness.
    * Example of 3 singled-ended replicates.
        ```javascript
        {
            "atac.paired_end" : false,
            "atac.bams" : ["rep1.bam", "rep2.bam", "rep3.bam"]
        }
        ```

3) Starting from filtered/deduped BAMs
    * Define a filtered/deduped BAM for each replicate. Our pipeline does not determine read endedness from a BAM file. You need to explicitly define read endedness. These BAMs should not have unmapped reads or duplicates.
    * Example of 2 singled-ended replicates.
        ```javascript
        {
            "atac.paired_end" : false,
            "atac.nodup_bams" : ["rep1.nodup.bam", "rep2.nodup.bam"]
        }
        ```

4) Starting from TAG-ALIGN BEDs
    * Define a TAG-ALIGN for each replicate. Our pipeline does not determine read endedness from a TAG-ALIGN file. You need to explicitly define read endedness.
    * Example of 4 paired-ended replicates.

        ```javascript
        {
            "atac.paired_end" : true,
            "atac.tas" : ["rep1.tagAlign.gz", "rep2.tagAlign.gz", "rep3.tagAlign.gz", "rep3.tagAlign.gz"]
        }
        ```

You can also mix up different data types for individual bio replicate. For example, pipeline can start from FASTQs for rep1 (SE) and rep3 (PE), BAMs for rep2 (SE), NODUP_BAMs for rep4 (SE) and TAG-ALIGNs for rep5 (PE).

```javascript
{
    "atac.paired_ends" : [false, false, true, false, true],
    "atac.fastqs_rep1_R1" : ["rep1.fastq.gz"],
    "atac.fastqs_rep3_R1" : ["rep3.R1.fastq.gz"],
    "atac.fastqs_rep3_R2" : ["rep3.R2.fastq.gz"],
    "atac.bams" : [null, "rep2.bam", null, null, null],
    "atac.nodup_bams" : [null, null, null, "rep4.nodup.bam", null],
    "atac.tas" : [null, null, null, null, "rep5.tagAlign.gz"]
}
```

## Adapters

1) Using automatic adapter detection
    * `atac.auto_detect_adapter`: Our pipeline can detect the following adpaters automatically from FASTQs.
        * Illumina: `AGATCGGAAGAGC`
        * Nextera: `CTGTCTCTTATA`
        * smallRNA: `TGGAATTCTCGG`

2) Define an adapter for **ALL FASTQs**
    * `atac.adapter`: This string will be used for **ALL FASTQs** unless you define an adapter individually for each FASTQ. Automatic adapter detection will be disabled.

3) Define an adapter individually for each FASTQ.
    * `atac.adapters_repX_RY`: You need to keep the same structure between two arrays (adapters and FASTQs).
    * The following is the same 3-bio-rep-2-tech-rep example used in the [previous section](#input-files). Since adapters for bio-rep 2 are blank. Pipeline will not trim adapters for these bio-rep 2 FASTQs unless you turn on the auto-detect-adapter flag (`atac.auto_detect_adapter`).

        ```json
            "atac.adapters_rep1_R1" : ["AAAAAAAAAA", "CCCCCCCCCC"],
            "atac.adapters_rep1_R2" : ["TTTTTTTTTT", "GGGGGGGGGG"],
            "atac.adapters_rep2_R1" : ["", ""],
            "atac.adapters_rep2_R2" : ["", ""],
            "atac.adapters_rep3_R1" : ["AAAAAAAAAA", "CCCCCCCCCC"],
            "atac.adapters_rep3_R2" : ["TTTTTTTTTT", "GGGGGGGGGG"],
        ```

## Resources

> **WARNING**: It is recommened not to change the following parameters unless you get resource-related errors for a certain task and you want to increase resources for such task. The following parameters are provided for users who want to run our pipeline with Caper's `local` on HPCs and 2).

Resources defined here are **PER BIO REPLICATE**. Therefore, total number of cores will be approximately `atac.align_cpu` x `NUMBER_OF_BIO_REPLICATES` because `align` is a bottlenecking task of the pipeline. This total number of cores will be useful **ONLY** when you use a `local` backend of Caper and manually `qsub` or `sbatch` your job. `disk_factor` is used for GCP/AWS/DNAnexus only.

See [this](input.md#resource-parameters) for details.
