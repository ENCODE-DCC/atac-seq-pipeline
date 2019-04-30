# resumer

## Introduction

This python script parses a metadata JSON file from a previous failed workflow and generates a new input JSON file to start a pipeline from where it left off.

## How to use it

Before running this script, you should have a metadata JSON file for a previous failed workflow. You can get it by adding a parameter `-m metadata.json` to the cromwell Java command line. If you stop a workflow (CTRL+C or kill) metadata then JSON file will not be generated.
```bash
$ java -jar ... cromwell-38.jar run chip.wdl -i original_input.json ... -m metadata.json
```

Unfortunately your workflow failed for some reasons but you can fix the problem and want to resume it from where it left off.
```bash
$ python resumer.py metadata.json
```

You will get a new input JSON file `resume.FAILED_WORKFLOW_ID.json` and run cromwell with it instead of the original one `original_input.json`.
```bash
$ java -jar ... cromwell-38.jar run chip.wdl -i resume.FAILED_WORKFLOW_ID.json ...
```

## Usage

```bash
usage: Resumer for ENCODE ATAC/Chip-Seq pipelines [-h]
                                                  [--output-def-json-file OUTPUT_DEF_JSON_FILE]
                                                  metadata_json_file

Parse cromwell's metadata JSON file and generate a new input JSON file to
resume a pipeline from where it left off.

positional arguments:
  metadata_json_file    Cromwell metadata JSON file from a previous failed
                        run.

optional arguments:
  -h, --help            show this help message and exit
  --output-def-json-file OUTPUT_DEF_JSON_FILE
                        Output definition JSON file for your pipeline. If not
                        specified, it will look for a valid JSON file on
                        script's directory. You can use your own JSON file for
                        your pipeline. Entries in "Array[Object]" is for
                        Array[Object] in an input JSON. This is useful to take
                        outputs from a scatter block. For example, the 1st
                        entry of "Array[Object]" in chip.json is "chip.bwa" :
                        {"bam" : "chip.bams", "flagstat_qc" :
                        "chip.flagstat_qcs"}. chip.flagstat_qcs : [...(taken
                        from an output of chip.bwa.flagstat_qc)...] will be
                        added to your new input JSON. For example, the 1st
                        entry of "Object" in chip.json is "chip.pool_ta" :
                        {"ta_pooled" : "chip.ta_pooled"}. chip.ta_pooled :
                        "(taken from an output of chip.pool_ta.ta_pooled)"
                        will be added to your new input JSON.
```

## Examples

```bash
$ python resumer.py metadata.json
```

## How it works (for developers)

In order to use this script, your pipeline should be able to start from any type of inputs (e.g. FASTQ, BAM, ...) and inputs to the previous task (e.g. map_fastq) should be ignored if next step (e.g. filter_bam)'s input is already given in the input JSON file.

```
# example toy_chip workflow that processes through FASTQ->BAM->FILT_BAM->PEAK->REPORT
# this pipeline can start from any types of input FASTQ, BAM, FILT_BAM, PEAK
# key idea of resuming workflow is to skip previous step 
# if next step's input is already given in the input JSON file
# this is controlled by `Boolean` variables (`need_to_process_XXX`).

workflow toy_chip {
	# input definition
	Array[File] fastqs = [] # per replicate
	Array[File] bams = [] # per replicate
	Array[File] filt_bams = [] # per replicate
	Array[File] peaks = [] # per replicate

	Boolean need_to_process_peak = true # trivial
	Boolean need_to_process_filt_bam = need_to_process_peak && length(peaks)==0
	Boolean need_to_process_bam = need_to_process_filt_bam && length(filt_bams)==0
	Boolean need_to_process_fastq = need_to_process_bam && length(bams)==0

	scatter(fastq in if need_to_process_fastq then fastqs else []) {
		call map_fastq { input: fastq = fastq }
	}

	# temporary array to deal with outputs from either previous step or from an input JSON file
	Array[File] bams_ = flatten([map_fastq.bam, bams]) 
	scatter(bam in if need_to_process_bam then bams_ else []) {
		call filter_bam { input: bam = bam }
	}
	
	Array[File] filt_bams_ = flatten([filter_bam.filt_bam, filt_bams]) # temporary array again
	scatter(filt_bam in if need_to_process_filt_bam then filt_bams_ else []) {
		call call_peak { input: filt_bam = filt_bam }
	}
	
	Array[File] peaks_ = flatten([call_peak.peak, peaks]) # temporary array again
	if (need_to_process_peak) {
		call generate_report { input: peaks = peaks_ }
	}
}
```

Output definition JSON file `toy_chip.json` for the above example workflow should look like:
```javascript
{
    "Array[Object]" : {
        "toy_chip.map_fastq" : {
        	"bam" : "toy_chip.bams"
    	},
        "toy_chip.filter_bam" : {
        	"filt_bam" : "toy_chip.filt_bams"
    	}
        "toy_chip.call_peak" : {
        	"peak" : "toy_chip.peaks"
    	}
    }
}
```

An original input JSON file to start from fastqs.
```javscript
{
	"toy_chip.fastqs" : ["rep1.fastq.gz", "rep1.fastq.gz"]
}
```

Run a pipeline with this original input JSON.
```bash
$ java -jar cromwell-38.jar run toy_chip.wdl -i org_input.json -m metadata.json
```

Pipeline fails due to some errors in `call_peak` task. Run `resumer.py` to make a new input JSON file to resume.
```bash
$ python resumer.py metadata.json --output-def-json-file toy_chip.json
```

Then `result.WORKFLOW_ID.json` will be generated.
```javscript
{
	"toy_chip.fastqs" : ["rep1.fastq.gz", "rep1.fastq.gz"]
	"toy_chip.bams" : ["rep1.bam", "rep1.bam"]
	"toy_chip.filt_bams" : ["rep1.filt.bam", "rep1.filt.bam"]
}
```

You feed it to the cromwell java command line after fixing the problem. Then pipeline will start from ``scatter` block for `call_peak` tasks.
```bash
$ java -jar cromwell-38.jar run toy_chip.wdl -i resume.WORKFLOW_ID.json
```

## Output definition JSON file (for developers)

An output definition JSON file must have at least one object from `"Array[Object]"` and `"Object"`. It can have both. The following JSON is a simplified version of an output definition JSON file for ChIP-Seq pipeline (`chip.json`).
```javascript
{
    "Array[Object]" : {
        "chip.bwa" : {
            "bam" : "chip.bams",
            "flagstat_qc" : "chip.flagstat_qcs"
        }
    },

    "Object" : {
        "chip.pool_ta" : {
            "ta_pooled" : "chip.ta_pooled"
        }
    }
}
```

`"Array[Object]"` is useful to take an array of outputs from a `scatter` block and `"Object"` is good for taking a single value from any tasks.

Using this JSON file for `resumer.py` will add the following extra input data definitions to the original input JSON file.
```javascript
{
	"chip.bams" : [...(an array of values taken from chip.bwa.bam)...],
	"chip.flagstat_qcs" : [...(an array of values taken from chip.bwa.flagstat_qc)...],
	"chip.ta_pooled" : "...(a value taken from chip.pool_ta.ta_pooled)..."
}


