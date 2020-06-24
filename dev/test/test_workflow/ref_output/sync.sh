#!/bin/bash

gsutil -m rsync -r -d . gs://encode-pipeline-test-samples/encode-atac-seq-pipeline/ref_output
