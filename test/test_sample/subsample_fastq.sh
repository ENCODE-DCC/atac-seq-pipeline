#!/bin/bash

RATIO=0.005

# SE
cd /mnt/data/pipeline_test_samples/encode-atac-seq-pipeline/ENCSR889WQX/fastq/rep1
seqtk sample -s100 ENCFF325FCQ.fastq.gz ${RATIO} | gzip -nc > ENCFF325FCQ.subsampled.fastq.gz
seqtk sample -s100 ENCFF439VSY.fastq.gz ${RATIO} | gzip -nc > ENCFF439VSY.subsampled.fastq.gz
seqtk sample -s100 ENCFF683IQS.fastq.gz ${RATIO} | gzip -nc > ENCFF683IQS.subsampled.fastq.gz
seqtk sample -s100 ENCFF744CHW.fastq.gz ${RATIO} | gzip -nc > ENCFF744CHW.subsampled.fastq.gz
gsutil cp *.subsampled.fastq.gz gs://encode-pipeline-test-samples/encode-atac-seq-pipeline/ENCSR889WQX/fastq_subsampled/rep1

cd /mnt/data/pipeline_test_samples/encode-atac-seq-pipeline/ENCSR889WQX/fastq/rep2
seqtk sample -s100 ENCFF463QCX.fastq.gz ${RATIO} | gzip -nc > ENCFF463QCX.subsampled.fastq.gz
seqtk sample -s100 ENCFF992TSA.fastq.gz ${RATIO} | gzip -nc > ENCFF992TSA.subsampled.fastq.gz
gsutil cp *.subsampled.fastq.gz gs://encode-pipeline-test-samples/encode-atac-seq-pipeline/ENCSR889WQX/fastq_subsampled/rep2

# PE
cd /mnt/data/pipeline_test_samples/encode-atac-seq-pipeline/ENCSR356KRQ/fastq/rep1/pair1/
seqtk sample -s100 ENCFF106QGY.fastq.gz ${RATIO} | gzip -nc > ENCFF106QGY.subsampled.fastq.gz
seqtk sample -s100 ENCFF341MYG.fastq.gz ${RATIO} | gzip -nc > ENCFF341MYG.subsampled.fastq.gz
gsutil cp *.subsampled.fastq.gz gs://encode-pipeline-test-samples/encode-atac-seq-pipeline/ENCSR356KRQ/fastq_subsampled/rep1/pair1
cd /mnt/data/pipeline_test_samples/encode-atac-seq-pipeline/ENCSR356KRQ/fastq/rep1/pair2/
seqtk sample -s100 ENCFF248EJF.fastq.gz ${RATIO} | gzip -nc > ENCFF248EJF.subsampled.fastq.gz
seqtk sample -s100 ENCFF368TYI.fastq.gz ${RATIO} | gzip -nc > ENCFF368TYI.subsampled.fastq.gz
gsutil cp *.subsampled.fastq.gz gs://encode-pipeline-test-samples/encode-atac-seq-pipeline/ENCSR356KRQ/fastq_subsampled/rep1/pair2

cd /mnt/data/pipeline_test_samples/encode-atac-seq-pipeline/ENCSR356KRQ/fastq/rep2/pair1/
seqtk sample -s100 ENCFF193RRC.fastq.gz ${RATIO} | gzip -nc > ENCFF193RRC.subsampled.fastq.gz
seqtk sample -s100 ENCFF366DFI.fastq.gz ${RATIO} | gzip -nc > ENCFF366DFI.subsampled.fastq.gz
seqtk sample -s100 ENCFF641SFZ.fastq.gz ${RATIO} | gzip -nc > ENCFF641SFZ.subsampled.fastq.gz
seqtk sample -s100 ENCFF751XTV.fastq.gz ${RATIO} | gzip -nc > ENCFF751XTV.subsampled.fastq.gz
seqtk sample -s100 ENCFF859BDM.fastq.gz ${RATIO} | gzip -nc > ENCFF859BDM.subsampled.fastq.gz
seqtk sample -s100 ENCFF927LSG.fastq.gz ${RATIO} | gzip -nc > ENCFF927LSG.subsampled.fastq.gz
gsutil cp *.subsampled.fastq.gz gs://encode-pipeline-test-samples/encode-atac-seq-pipeline/ENCSR356KRQ/fastq_subsampled/rep2/pair1
cd /mnt/data/pipeline_test_samples/encode-atac-seq-pipeline/ENCSR356KRQ/fastq/rep2/pair2/
seqtk sample -s100 ENCFF007USV.fastq.gz ${RATIO} | gzip -nc > ENCFF007USV.subsampled.fastq.gz
seqtk sample -s100 ENCFF031ARQ.fastq.gz ${RATIO} | gzip -nc > ENCFF031ARQ.subsampled.fastq.gz
seqtk sample -s100 ENCFF573UXK.fastq.gz ${RATIO} | gzip -nc > ENCFF573UXK.subsampled.fastq.gz
seqtk sample -s100 ENCFF590SYZ.fastq.gz ${RATIO} | gzip -nc > ENCFF590SYZ.subsampled.fastq.gz
seqtk sample -s100 ENCFF734PEQ.fastq.gz ${RATIO} | gzip -nc > ENCFF734PEQ.subsampled.fastq.gz
seqtk sample -s100 ENCFF886FSC.fastq.gz ${RATIO} | gzip -nc > ENCFF886FSC.subsampled.fastq.gz
gsutil cp *.subsampled.fastq.gz gs://encode-pipeline-test-samples/encode-atac-seq-pipeline/ENCSR356KRQ/fastq_subsampled/rep2/pair2


