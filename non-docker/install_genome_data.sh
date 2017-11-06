#!/bin/bash
# Stop on error
set -e

if [[ "$#" -lt 2 ]]; then
  echo
  echo "This script installs data for genome [GENOME] on a directory [DEST_DIR]."
  echo "A TSV file [DEST_DIR]/[GENOME]_local.conf will be generated. Use it for pipeline."
  echo
  echo "Supported genomes: hg19, mm9, hg38 and mm10"
  echo
  echo "Usage: ./install_genome_data.sh [GENOME] [DEST_DIR]"
  echo "  Example: ./install_genome_data.sh hg38 /your/genome/data/path/hg38"
  echo
  exit 1
fi

# pipeline specific params
# CONDA_ENV="atac-seq-pipeline"
CONDA_ENV="bds_atac"
BUILD_BWT2_IDX=1
BUILD_BWA_IDX=0

GENOME=$1
#DEST_DIR=$(readlink -f $2)
DEST_DIR=$(cd $(dirname $2) && pwd -P)/$(basename $2)
TSV=${DEST_DIR}/${GENOME}_local.tsv

echo "=== Creating destination directory and TSV file..."
mkdir -p ${DEST_DIR}
cd ${DEST_DIR}

if [[ $GENOME == "hg19" ]]; then
  REF_FA="http://hgdownload.cse.ucsc.edu/goldenpath/hg19/encodeDCC/referenceSequences/male.hg19.fa.gz"
  BLACKLIST="http://hgdownload.cse.ucsc.edu/goldenpath/hg19/encodeDCC/wgEncodeMapability/wgEncodeDacMapabilityConsensusExcludable.bed.gz"

elif [[ $GENOME == "mm9" ]]; then
  REF_FA="http://hgdownload.cse.ucsc.edu/goldenPath/mm9/bigZips/mm9.2bit"
  BLACKLIST="http://mitra.stanford.edu/kundaje/genome_data/mm9/mm9-blacklist.bed.gz"

elif [[ $GENOME == "hg38" ]]; then
  REF_FA="https://www.encodeproject.org/files/GRCh38_no_alt_analysis_set_GCA_000001405.15/@@download/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta.gz"
  BLACKLIST="http://mitra.stanford.edu/kundaje/genome_data/hg38/hg38.blacklist.bed.gz"

elif [[ $GENOME == "mm10" ]]; then
  REF_FA="https://www.encodeproject.org/files/mm10_no_alt_analysis_set_ENCODE/@@download/mm10_no_alt_analysis_set_ENCODE.fasta.gz"
  BLACKLIST="http://mitra.stanford.edu/kundaje/genome_data/mm10/mm10.blacklist.bed.gz"
fi

if [[ ${REF_FA} == "" ]]; then
  echo "Error: unsupported genome $GENOME"
  exit 1
fi

echo "=== Downloading files..."
wget -c -O $(basename ${REF_FA}) ${REF_FA}
if [[ $BLACKLIST != "" ]]; then wget -N -c $BLACKLIST; fi

source activate ${CONDA_ENV}

echo "=== Processing reference fasta file..."
if [[ ${REF_FA} == *.gz ]]; then 
  REF_FA_PREFIX=$(basename ${REF_FA} .gz)
  gzip -d -f -c ${REF_FA_PREFIX}.gz > ${REF_FA_PREFIX}
elif [[ ${REF_FA} == *.bz2 ]]; then
  REF_FA_PREFIX=$(basename ${REF_FA} .bz2)
  bzip2 -d -f -c ${REF_FA_PREFIX}.bz2 > ${REF_FA_PREFIX}
elif [[ ${REF_FA} == *.2bit ]]; then
  REF_FA_PREFIX=$(basename ${REF_FA} .2bit).fa
  twoBitToFa $(basename ${REF_FA}) ${REF_FA_PREFIX}
else
  REF_FA_PREFIX=$(basename ${REF_FA})  
fi

echo "=== Generating fasta index and chrom.sizes file..."
cd ${DEST_DIR}
mkdir -p seq
cd seq
rm -f ${REF_FA_PREFIX}
ln -s ../${REF_FA_PREFIX} ${REF_FA_PREFIX}
faidx -x ${REF_FA_PREFIX}
cp -f *.fai ../
CHRSZ=$GENOME.chrom.sizes
cut -f1,2 ${REF_FA_PREFIX}.fai > ../$CHRSZ

echo "=== Determinig gensz..."
cd ${DEST_DIR}
GENSZ=$(cat $CHRSZ | awk '{sum+=$2} END{print sum}')
if [[ $GENOME == hg* ]]; then GENSZ=hs; fi
if [[ $GENOME == mm* ]]; then GENSZ=mm; fi

if [[ ${BUILD_BWT2_IDX} == 1 ]]; then
  echo "=== Building bowtie2 index..."
  mkdir -p ${DEST_DIR}/bowtie2_index
  cd ${DEST_DIR}/bowtie2_index
  rm -f ${REF_FA_PREFIX}
  ln -s ../${REF_FA_PREFIX} ${REF_FA_PREFIX}
  if [[ ! -f ${REF_FA_PREFIX}.rev.1.bt2 ]]; then
    bowtie2-build ${REF_FA_PREFIX} ${REF_FA_PREFIX}
    tar cvf ${REF_FA_PREFIX}.tar ${REF_FA_PREFIX}.*.bt2
  fi
fi

if [[ ${BUILD_BWA_IDX} == 1 ]]; then
  echo "=== Building bwa index..."
  mkdir -p ${DEST_DIR}/bwa_index
  cd ${DEST_DIR}/bwa_index
  rm -f ${REF_FA_PREFIX}
  ln -s ../${REF_FA_PREFIX} ${REF_FA_PREFIX}
  if [[ ! -f ${REF_FA_PREFIX}.sa ]]; then
    bwa index ${REF_FA_PREFIX}
    tar cvf ${REF_FA_PREFIX}.tar ${REF_FA_PREFIX}.*
  fi
fi

echo "=== Creating TSV file... (${TSV})"
cd ${DEST_DIR}
rm -f ${TSV}
touch ${TSV}

if [[ $BLACKLIST != "" ]]; then
  echo -e "blacklist\t${DEST_DIR}/$(basename $BLACKLIST)" >> ${TSV};
fi
echo -e "chrsz\t${DEST_DIR}/$(basename $CHRSZ)" >> ${TSV}
echo -e "gensz\t$GENSZ" >> ${TSV}
if [[ ${BUILD_BWT2_IDX} == 1 ]]; then
  echo -e "bowtie2_idx_tar\t${DEST_DIR}/bowtie2_index/${REF_FA_PREFIX}" >> ${TSV}
fi
if [[ ${BUILD_BWA_IDX} == 1 ]]; then
  echo -e "bwa_idx_tar\t${DEST_DIR}/bwa_index/${REF_FA_PREFIX}" >> ${TSV}
fi
echo -e "ref_fa\t${DEST_DIR}/${REF_FA_PREFIX}" >> ${TSV}

echo "=== All done."
