source activate chip-seq-pipeline

WORKDIR=/mnt/data/pipeline_test_samples/atac/ENCSR889WQX_chr19
SAMPLE=ENCSR889WQX
CHR=chr19
BAM1=/srv/scratch/shared/kadru/leepc12/run/atac-seq-pipeline/bds/ENCSR889WQX/out/align/rep1/ENCFF439VSY.trim_pooled.bam
BAM2=/srv/scratch/shared/kadru/leepc12/run/atac-seq-pipeline/bds/ENCSR889WQX/out/align/rep2/ENCFF463QCX.trim_ENCFF992TSA.trim.bam
SMALL_BAM1=$SAMPLE.$CHR.rep1.bam
SMALL_BAM2=$SAMPLE.$CHR.rep2.bam
SMALL_FQ1=$SAMPLE.$CHR.rep1.fastq
SMALL_FQ2=$SAMPLE.$CHR.rep2.fastq
mkdir -p $WORKDIR && cd $WORKDIR
samtools view -b $BAM1 "$CHR" > $SMALL_BAM1
samtools view -b $BAM2 "$CHR" > $SMALL_BAM2
bamToFastq -i $SMALL_BAM1 -fq $SMALL_FQ1
bamToFastq -i $SMALL_BAM2 -fq $SMALL_FQ2
gzip -f $SMALL_FQ1 $SMALL_FQ2

WORKDIR=/mnt/data/pipeline_test_samples/atac/ENCSR356KRQ/fastq_chr19
SAMPLE=ENCSR356KRQ
CHR=chr19
BAM1=/srv/scratch/shared/kadru/leepc12/run/atac-seq-pipeline/bds/ENCSR356KRQ/out/align/rep1/ENCFF341MYG.trim_ENCFF106QGY.trim.PE2SE.bam
BAM2=/srv/scratch/shared/kadru/leepc12/run/atac-seq-pipeline/bds/ENCSR356KRQ/out/align/rep2/ENCFF641SFZ.trim_pooled.PE2SE.bam
SMALL_BAM1=$SAMPLE.$CHR.rep1.bam
SMALL_BAM2=$SAMPLE.$CHR.rep2.bam
SMALL_NMSRT_BAM1=$SAMPLE.$CHR.rep1.nmsrt.bam
SMALL_NMSRT_BAM2=$SAMPLE.$CHR.rep2.nmsrt.bam
SMALL_FQ1_R1=$SAMPLE.$CHR.rep1-R1.fastq
SMALL_FQ1_R2=$SAMPLE.$CHR.rep1-R2.fastq
SMALL_FQ2_R1=$SAMPLE.$CHR.rep2-R1.fastq
SMALL_FQ2_R2=$SAMPLE.$CHR.rep2-R2.fastq
mkdir -p $WORKDIR && cd $WORKDIR
samtools view -b $BAM1 "$CHR" > $SMALL_BAM1
samtools view -b $BAM2 "$CHR" > $SMALL_BAM2
samtools sort -n $SMALL_BAM1 $SMALL_NMSRT_BAM1
samtools sort -n $SMALL_BAM2 $SMALL_NMSRT_BAM2
samtools bam2fq $SMALL_NMSRT_BAM1.bam > tmp.fq1
samtools bam2fq $SMALL_NMSRT_BAM2.bam > tmp.fq2
cat tmp.fq1 | grep '^@.*/1$' -A 3 --no-group-separator > $SMALL_FQ1_R1
cat tmp.fq1 | grep '^@.*/2$' -A 3 --no-group-separator > $SMALL_FQ1_R2
cat tmp.fq2 | grep '^@.*/1$' -A 3 --no-group-separator > $SMALL_FQ2_R1
cat tmp.fq2 | grep '^@.*/2$' -A 3 --no-group-separator > $SMALL_FQ2_R2
#bamToFastq -i $SMALL_NMSRT_BAM1.bam -fq $SMALL_FQ1_R1 -fq2 $SMALL_FQ1_R2
#bamToFastq -i $SMALL_NMSRT_BAM2.bam -fq $SMALL_FQ2_R1 -fq2 $SMALL_FQ2_R2
gzip -f $SMALL_FQ1_R1 $SMALL_FQ1_R2 $SMALL_FQ2_R1 $SMALL_FQ2_R2
rm -f $SMALL_BAM1 $SMALL_BAM2 $SMALL_NMSRT_BAM1.bam $SMALL_NMSRT_BAM2.bam tmp.fq1 tmp.fq2



