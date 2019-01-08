# qc_jsons_to_tsv

## Introduction

This python script finds `qc.json` (or can be specified) recursively, parses all of them and make a TSV spreadsheet file.

## Usage

```
usage: qc.json parser for ENCODE ATAC/Chip-Seq pipelines [-h]
                                                         [--search-dir SEARCH_DIR]
                                                         [--out-file OUT_FILE]
                                                         [--criteria-def-json-file CRITERIA_DEF_JSON_FILE]
                                                         [--qc-json-file-basename QC_JSON_FILE_BASENAME]
                                                         [--tsv-mapping-qc-json-path-to-title TSV_MAPPING_QC_JSON_PATH_TO_TITLE]
                                                         [--tsv-mapping-workflow-id-to-title TSV_MAPPING_WORKFLOW_ID_TO_TITLE]

Recursively find qc.json, parse it and make a TSV spreadsheet of all quality
metrics.

optional arguments:
  -h, --help            show this help message and exit
  --search-dir SEARCH_DIR
                        Root directory to be recursively searched for qc.json
                        (or --qc-json-file-basename).
  --out-file OUT_FILE   (Optional) Output TSV filename. Prints to STDOUT if
                        not defined.
  --criteria-def-json-file CRITERIA_DEF_JSON_FILE
                        (Optional but important) Specify criteria definition
                        JSON file. "criteria" category will be added to the
                        output file/STDOUT.
  --qc-json-file-basename QC_JSON_FILE_BASENAME
                        (Optional) Specify QC JSON file basename to be parsed.
                        Files with this name will be recursively found and
                        parsed.
  --tsv-mapping-qc-json-path-to-title TSV_MAPPING_QC_JSON_PATH_TO_TITLE
                        (Optional) Two-column TSV (ABSOLUTE path for qc.json
                        (or --qc-json-file-basename) [TAB] title).This is
                        useful when you forgot to define titles for your
                        pipeline runs.
  --tsv-mapping-workflow-id-to-title TSV_MAPPING_WORKFLOW_ID_TO_TITLE
                        (Optional) Two-column TSV (cromwell workflow ID [TAB]
                        title).This is useful when you forgot to define titles
                        for your pipeline runs.
```

## Examples

```
python qc_jsons_to_tsv.py --search-dir test/v1.1.4 --criteria-def-json-file criteria.default.json > test_v1.1.4.tsv
python qc_jsons_to_tsv.py --search-dir test/v1.1.5 --criteria-def-json-file criteria.default.json > test_v1.1.5.tsv
```