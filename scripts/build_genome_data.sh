#!/bin/bash
# Stop on error
set -e

if [[ "$#" -lt 2 ]]; then
  echo
  echo "ACTIVATE PIPELINE'S CONDA ENVIRONMENT BEFORE RUNNING THIS SCRIPT!"
  echo
  echo "This script installs data for genome [GENOME] on a directory [DEST_DIR]."
  echo "A TSV file [DEST_DIR]/[GENOME].tsv will be generated. Use it for pipeline."
  echo
  echo "Supported genomes: hg19, mm9, hg38 and mm10"
  echo
  echo "Usage: ./build_genome_data.sh [GENOME] [DEST_DIR]"
  echo "  Example: ./build_genome_data.sh hg38 /your/genome/data/path/hg38"
  echo
  exit 2
fi

# pipeline specific params
BUILD_BWT2_IDX=1
BUILD_BWA_IDX=0

GENOME=$1
DEST_DIR=$(cd $(dirname $2) && pwd -P)/$(basename $2)
TSV=${DEST_DIR}/${GENOME}.tsv

echo "=== Creating destination directory and TSV file..."
mkdir -p ${DEST_DIR}
cd ${DEST_DIR}

if [[ "${GENOME}" == "hg19" ]]; then
  MITO_CHR_NAME="chrM"
  REF_FA="http://hgdownload.cse.ucsc.edu/goldenpath/hg19/encodeDCC/referenceSequences/male.hg19.fa.gz"
  BLACKLIST="http://hgdownload.cse.ucsc.edu/goldenpath/hg19/encodeDCC/wgEncodeMapability/wgEncodeDacMapabilityConsensusExcludable.bed.gz"
  # optional data
  TSS="https://storage.googleapis.com/encode-pipeline-genome-data/hg19/ataqc/hg19_gencode_tss_unique.bed.gz"
  DNASE="https://storage.googleapis.com/encode-pipeline-genome-data/hg19/ataqc/reg2map_honeybadger2_dnase_all_p10_ucsc.bed.gz"
  PROM="https://storage.googleapis.com/encode-pipeline-genome-data/hg19/ataqc/reg2map_honeybadger2_dnase_prom_p2.bed.gz"
  ENH="https://storage.googleapis.com/encode-pipeline-genome-data/hg19/ataqc/reg2map_honeybadger2_dnase_enh_p2.bed.gz"
  REG2MAP="https://storage.googleapis.com/encode-pipeline-genome-data/hg19/ataqc/dnase_avgs_reg2map_p10_merged_named.pvals.gz"
  ROADMAP_META="https://storage.googleapis.com/encode-pipeline-genome-data/hg19/ataqc/eid_to_mnemonic.txt"

elif [[ "${GENOME}" == "mm9" ]]; then
  MITO_CHR_NAME="chrM"
  REF_FA="http://hgdownload.cse.ucsc.edu/goldenPath/mm9/bigZips/mm9.2bit"
  BLACKLIST="https://storage.googleapis.com/encode-pipeline-genome-data/mm9/mm9-blacklist.bed.gz"
  # optional data
  TSS="https://storage.googleapis.com/encode-pipeline-genome-data/mm9/ataqc/mm9_gencode_tss_unique.bed.gz"
  DNASE="https://storage.googleapis.com/encode-pipeline-genome-data/mm9/ataqc/mm9_univ_dhs_ucsc.from_mm10.bed.gz"
  PROM="https://storage.googleapis.com/encode-pipeline-genome-data/mm9/ataqc/tss_mm9_master.from_mm10.bed.gz"
  ENH="https://storage.googleapis.com/encode-pipeline-genome-data/mm9/ataqc/mm9_enh_dhs_ucsc.from_mm10.bed.gz"
  REG2MAP_BED="https://storage.googleapis.com/encode-pipeline-genome-data/mm9/ataqc/mm9_dhs_universal_ucsc_v1.bed.gz"
  REG2MAP="https://storage.googleapis.com/encode-pipeline-genome-data/mm9/ataqc/dnase_avgs_merged_named.fseq.vals.gz"
  ROADMAP_META="https://storage.googleapis.com/encode-pipeline-genome-data/mm9/ataqc/accession_to_name.txt"

elif [[ "${GENOME}" == "hg38" ]]; then
  MITO_CHR_NAME="chrM"
  REF_FA="https://www.encodeproject.org/files/GRCh38_no_alt_analysis_set_GCA_000001405.15/@@download/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta.gz"
  BLACKLIST="https://storage.googleapis.com/encode-pipeline-genome-data/hg38/hg38.blacklist.bed.gz"
  # optional data
  TSS="https://storage.googleapis.com/encode-pipeline-genome-data/hg38/ataqc/hg38_gencode_tss_unique.bed.gz"
  DNASE="https://storage.googleapis.com/encode-pipeline-genome-data/hg38/ataqc/reg2map_honeybadger2_dnase_all_p10_ucsc.hg19_to_hg38.bed.gz"
  PROM="https://storage.googleapis.com/encode-pipeline-genome-data/hg38/ataqc/reg2map_honeybadger2_dnase_prom_p2.hg19_to_hg38.bed.gz"
  ENH="https://storage.googleapis.com/encode-pipeline-genome-data/hg38/ataqc/reg2map_honeybadger2_dnase_enh_p2.hg19_to_hg38.bed.gz"
  REG2MAP_BED="https://storage.googleapis.com/encode-pipeline-genome-data/hg38/ataqc/hg38_celltype_compare_subsample.bed.gz"
  REG2MAP="https://storage.googleapis.com/encode-pipeline-genome-data/hg38/ataqc/hg38_dnase_avg_fseq_signal_formatted.txt.gz"
  ROADMAP_META="https://storage.googleapis.com/encode-pipeline-genome-data/hg38/ataqc/hg38_dnase_avg_fseq_signal_metadata.txt"

elif [[ "${GENOME}" == "mm10" ]]; then
  MITO_CHR_NAME="chrM"
  REF_FA="https://www.encodeproject.org/files/mm10_no_alt_analysis_set_ENCODE/@@download/mm10_no_alt_analysis_set_ENCODE.fasta.gz"
  BLACKLIST="https://storage.googleapis.com/encode-pipeline-genome-data/mm10/mm10.blacklist.bed.gz"
  # optional data
  TSS="https://storage.googleapis.com/encode-pipeline-genome-data/mm10/ataqc/mm10_gencode_tss_unique.bed.gz"
  DNASE="https://storage.googleapis.com/encode-pipeline-genome-data/mm10/ataqc/mm10_univ_dhs_ucsc.bed.gz"
  PROM="https://storage.googleapis.com/encode-pipeline-genome-data/mm10/ataqc/tss_mm10_master.bed.gz"
  ENH="https://storage.googleapis.com/encode-pipeline-genome-data/mm10/ataqc/mm10_enh_dhs_ucsc.bed.gz"
  REG2MAP_BED="https://storage.googleapis.com/encode-pipeline-genome-data/mm10/ataqc/mm10_celltype_compare_subsample.bed.gz"
  REG2MAP="https://storage.googleapis.com/encode-pipeline-genome-data/mm10/ataqc/mm10_dnase_avg_fseq_signal_formatted.txt.gz"
  ROADMAP_META="https://storage.googleapis.com/encode-pipeline-genome-data/mm10/ataqc/mm10_dnase_avg_fseq_signal_metadata.txt"

elif [[ "${GENOME}" == "hg38_chr19_chrM" ]]; then
  MITO_CHR_NAME="chrM"
  REF_FA="https://storage.googleapis.com/encode-pipeline-genome-data/hg38_chr19_chrM/GRCh38_no_alt_analysis_set_GCA_000001405.15.chr19_chrM.fasta.gz"
  BLACKLIST="https://storage.googleapis.com/encode-pipeline-genome-data/hg38/hg38.blacklist.bed.gz"
  # optional data
  TSS="https://storage.googleapis.com/encode-pipeline-genome-data/hg38/ataqc/hg38_gencode_tss_unique.bed.gz"
  DNASE="https://storage.googleapis.com/encode-pipeline-genome-data/hg38/ataqc/reg2map_honeybadger2_dnase_all_p10_ucsc.hg19_to_hg38.bed.gz"
  PROM="https://storage.googleapis.com/encode-pipeline-genome-data/hg38/ataqc/reg2map_honeybadger2_dnase_prom_p2.hg19_to_hg38.bed.gz"
  ENH="https://storage.googleapis.com/encode-pipeline-genome-data/hg38/ataqc/reg2map_honeybadger2_dnase_enh_p2.hg19_to_hg38.bed.gz"
  REG2MAP_BED="https://storage.googleapis.com/encode-pipeline-genome-data/hg38/ataqc/hg38_celltype_compare_subsample.bed.gz"
  REG2MAP="https://storage.googleapis.com/encode-pipeline-genome-data/hg38/ataqc/hg38_dnase_avg_fseq_signal_formatted.txt.gz"
  ROADMAP_META="https://storage.googleapis.com/encode-pipeline-genome-data/hg38/ataqc/hg38_dnase_avg_fseq_signal_metadata.txt"

elif [[ "${GENOME}" == "mm10_chr19_chrM" ]]; then
  MITO_CHR_NAME="chrM"
  REF_FA="https://storage.googleapis.com/encode-pipeline-genome-data/mm10_chr19_chrM/mm10_no_alt_analysis_set_ENCODE.chr19_chrM.fasta.gz"
  BLACKLIST="https://storage.googleapis.com/encode-pipeline-genome-data/mm10/mm10.blacklist.bed.gz"
  # optional data
  TSS="https://storage.googleapis.com/encode-pipeline-genome-data/mm10/ataqc/mm10_gencode_tss_unique.bed.gz"
  DNASE="https://storage.googleapis.com/encode-pipeline-genome-data/mm10/ataqc/mm10_univ_dhs_ucsc.bed.gz"
  PROM="https://storage.googleapis.com/encode-pipeline-genome-data/mm10/ataqc/tss_mm10_master.bed.gz"
  ENH="https://storage.googleapis.com/encode-pipeline-genome-data/mm10/ataqc/mm10_enh_dhs_ucsc.bed.gz"
  REG2MAP_BED="https://storage.googleapis.com/encode-pipeline-genome-data/mm10/ataqc/mm10_celltype_compare_subsample.bed.gz"
  REG2MAP="https://storage.googleapis.com/encode-pipeline-genome-data/mm10/ataqc/mm10_dnase_avg_fseq_signal_formatted.txt.gz"
  ROADMAP_META="https://storage.googleapis.com/encode-pipeline-genome-data/mm10/ataqc/mm10_dnase_avg_fseq_signal_metadata.txt"

elif [[ "${GENOME}" == "YOUR_OWN_GENOME" ]]; then
  MITO_CHR_NAME="chrM"
  REF_FA="URL_FOR_YOUR_FASTA_OR_2BIT"
  BLACKLIST= # leave it empty if you don't have it
fi

if [[ -z "${REF_FA}" ]]; then
  echo "Error: unsupported genome $GENOME"
  exit 1
fi
if [[ -z "${MITO_CHR_NAME}" ]]; then
  echo "Error: Mitochondrial chromosome name must be defined"
  exit 1
fi

echo "=== Downloading files..."
wget -c -O $(basename ${REF_FA}) ${REF_FA}
if [[ ! -z "${BLACKLIST}" ]]; then wget -N -c ${BLACKLIST}; fi

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
gzip -nc ${REF_FA_PREFIX} > ${REF_FA_PREFIX}.gz

echo "=== Generating fasta index and chrom.sizes file..."
cd ${DEST_DIR}
samtools faidx ${REF_FA_PREFIX}
CHRSZ=$GENOME.chrom.sizes
cut -f1,2 ${REF_FA_PREFIX}.fai > ${CHRSZ}

echo "=== Extracting mito chromosome from fasta"
REF_FA_PREFIX_WO_EXT=${REF_FA_PREFIX%.*}
REF_MITO_FA_PREFIX=${REF_FA_PREFIX_WO_EXT}.${MITO_CHR_NAME}.fa
samtools faidx ${REF_FA_PREFIX} ${MITO_CHR_NAME} > ${REF_MITO_FA_PREFIX}
gzip -nc ${REF_MITO_FA_PREFIX} > ${REF_MITO_FA_PREFIX}.gz

echo "=== Determinig gensz..."
cd ${DEST_DIR}
GENSZ=$(cat $CHRSZ | awk '{sum+=$2} END{print sum}')
if [[ "${GENOME}" == hg* ]]; then GENSZ=hs; fi
if [[ "${GENOME}" == mm* ]]; then GENSZ=mm; fi

if [[ ${BUILD_BWT2_IDX} == 1 ]]; then
  echo "=== Building bowtie2 index..."
  mkdir -p ${DEST_DIR}/bowtie2_index
  cd ${DEST_DIR}/bowtie2_index

  # whole chr
  rm -f ${REF_FA_PREFIX}
  ln -s ../${REF_FA_PREFIX} ${REF_FA_PREFIX}
  if [[ ! -f ${REF_FA_PREFIX}.rev.1.bt2 ]]; then
    bowtie2-build ${REF_FA_PREFIX} ${REF_FA_PREFIX}
    tar cvf ${REF_FA_PREFIX}.tar ${REF_FA_PREFIX}.*.bt2
  fi
  # mito chr only
  rm -f ${REF_MITO_FA_PREFIX}
  ln -s ../${REF_MITO_FA_PREFIX} ${REF_MITO_FA_PREFIX}
  if [[ ! -f ${REF_MITO_FA_PREFIX}.rev.1.bt2 ]]; then
    bowtie2-build ${REF_MITO_FA_PREFIX} ${REF_MITO_FA_PREFIX}
    tar cvf ${REF_MITO_FA_PREFIX}.tar ${REF_MITO_FA_PREFIX}.*.bt2
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
  ln -s ../${REF_MITO_FA_PREFIX} ${REF_MITO_FA_PREFIX}
  if [[ ! -f ${REF_MITO_FA_PREFIX}.sa ]]; then
    bwa index ${REF_MITO_FA_PREFIX}
    tar cvf ${REF_MITO_FA_PREFIX}.tar ${REF_MITO_FA_PREFIX}.*
  fi
fi

echo "=== Creating TSV file... (${TSV})"
cd ${DEST_DIR}
rm -f ${TSV}
touch ${TSV}

echo -e "genome_name\t${GENOME}" >> ${TSV}
echo -e "ref_fa\t${DEST_DIR}/${REF_FA_PREFIX}.gz" >> ${TSV}
echo -e "ref_mito_fa\t${DEST_DIR}/${REF_MITO_FA_PREFIX}.gz" >> ${TSV}
if [[ ! -z "${MITO_CHR_NAME}" ]]; then
  echo -e "mito_chr_name\t${MITO_CHR_NAME}" >> ${TSV}
fi
if [[ ! -z "${BLACKLIST}" ]]; then
  echo -e "blacklist\t${DEST_DIR}/$(basename ${BLACKLIST})" >> ${TSV};
fi
echo -e "chrsz\t${DEST_DIR}/$(basename ${CHRSZ})" >> ${TSV}
echo -e "gensz\t${GENSZ}" >> ${TSV}
if [[ ${BUILD_BWT2_IDX} == 1 ]]; then
  echo -e "bowtie2_idx_tar\t${DEST_DIR}/bowtie2_index/${REF_FA_PREFIX}.tar" >> ${TSV}
  echo -e "bowtie2_mito_idx_tar\t${DEST_DIR}/bowtie2_index/${REF_MITO_FA_PREFIX}.tar" >> ${TSV}
fi
if [[ ${BUILD_BWA_IDX} == 1 ]]; then
  echo -e "bwa_idx_tar\t${DEST_DIR}/bwa_index/${REF_FA_PREFIX}.tar" >> ${TSV}
  echo -e "bwa_mito_idx_tar\t${DEST_DIR}/bwa_index/${REF_MITO_FA_PREFIX}.tar" >> ${TSV}
fi

echo "=== Downloading ATAQC file... (${TSV})"
cd ${DEST_DIR}
mkdir -p ataqc
cd ataqc

if [[ ! -z "${TSS}" ]]; then
  wget -N -c ${TSS}
  echo -e "tss\t${DEST_DIR}/ataqc/$(basename ${TSS})" >> ${TSV}
fi
if [[ ! -z "${DNASE}" ]]; then
  wget -N -c ${DNASE}
  echo -e "dnase\t${DEST_DIR}/ataqc/$(basename ${DNASE})" >> ${TSV}
fi
if [[ ! -z "${PROM}" ]]; then
  wget -N -c ${PROM}
  echo -e "prom\t${DEST_DIR}/ataqc/$(basename ${PROM})" >> ${TSV}
fi
if [[ ! -z "${ENH}" ]]; then
  wget -N -c ${ENH}
  echo -e "enh\t${DEST_DIR}/ataqc/$(basename ${ENH})" >> ${TSV}
fi
if [[ ! -z "${REG2MAP}" ]]; then
  wget -N -c ${REG2MAP}
  echo -e "reg2map\t${DEST_DIR}/ataqc/$(basename ${REG2MAP})" >> ${TSV}
fi
if [[ ! -z "${REG2MAP_BED}" ]]; then
  wget -N -c ${REG2MAP_BED}
  echo -e "reg2map_bed\t${DEST_DIR}/ataqc/$(basename ${REG2MAP_BED})" >> ${TSV}
fi
if [[ ! -z "${ROADMAP_META}" ]]; then
  wget -N -c ${ROADMAP_META}
  echo -e "roadmap_meta\t${DEST_DIR}/ataqc/$(basename ${ROADMAP_META})" >> ${TSV}
fi

echo "=== All done."

