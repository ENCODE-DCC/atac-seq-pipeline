#!/bin/bash
# Stop on error
set -e

if [[ "$#" -lt 2 ]]; then
  echo
  echo "This script downloads/installs data for genome [GENOME] on a directory [DEST_DIR]."
  echo "A TSV file [DEST_DIR]/[GENOME].tsv will be generated. Use it for pipeline."
  echo
  echo "Supported genomes: hg19, mm9, hg38 and mm10"
  echo
  echo "Usage: ./download_genome_data.sh [GENOME] [DEST_DIR]"
  echo "  Example: ./download_genome_data.sh hg38 /your/genome/data/path/hg38"
  echo
  exit 2
fi

# parameters for genome database version (v1: <ENCODE4, v3: >=ENCODE4)
VER=v3

GENOME=$1
DEST_DIR=$(cd $(dirname $2) && pwd -P)/$(basename $2)
TSV=${DEST_DIR}/${GENOME}.tsv

echo "=== Creating destination directory and TSV file..."
mkdir -p ${DEST_DIR}
cd ${DEST_DIR}

if [[ "${GENOME}" == "hg19" ]]; then
  REGEX_BFILT_PEAK_CHR_NAME="chr[\dXY]+"
  MITO_CHR_NAME="chrM"
  REF_FA="http://hgdownload.cse.ucsc.edu/goldenpath/hg19/encodeDCC/referenceSequences/male.hg19.fa.gz"
  REF_MITO_FA="https://storage.googleapis.com/encode-pipeline-genome-data/hg19/male.hg19.chrM.fa.gz"
  CHRSZ="https://storage.googleapis.com/encode-pipeline-genome-data/hg19/hg19.chrom.sizes"
  BWT2_IDX="https://storage.googleapis.com/encode-pipeline-genome-data/hg19/bowtie2_index/male.hg19.fa.tar"
  BWT2_MITO_IDX="https://storage.googleapis.com/encode-pipeline-genome-data/hg19/bowtie2_index/male.hg19.chrM.fa.tar"
  BWA_IDX="https://storage.googleapis.com/encode-pipeline-genome-data/hg19/bwa_index/male.hg19.fa.tar"
  BWA_MITO_IDX="https://storage.googleapis.com/encode-pipeline-genome-data/hg19/bwa_index/male.hg19.chrM.fa.tar"
  BLACKLIST="https://storage.googleapis.com/encode-pipeline-genome-data/hg19/wgEncodeDacMapabilityConsensusExcludable.bed.gz"

  TSS="https://storage.googleapis.com/encode-pipeline-genome-data/hg19/ataqc/hg19_gencode_tss_unique.bed.gz"
  DNASE="https://storage.googleapis.com/encode-pipeline-genome-data/hg19/ataqc/reg2map_honeybadger2_dnase_all_p10_ucsc.bed.gz"
  PROM="https://storage.googleapis.com/encode-pipeline-genome-data/hg19/ataqc/reg2map_honeybadger2_dnase_prom_p2.bed.gz"
  ENH="https://storage.googleapis.com/encode-pipeline-genome-data/hg19/ataqc/reg2map_honeybadger2_dnase_enh_p2.bed.gz"
  REG2MAP="https://storage.googleapis.com/encode-pipeline-genome-data/hg19/ataqc/dnase_avgs_reg2map_p10_merged_named.pvals.gz"
  ROADMAP_META="https://storage.googleapis.com/encode-pipeline-genome-data/hg19/ataqc/eid_to_mnemonic.txt"

elif [[ "${GENOME}" == "mm9" ]]; then
  REGEX_BFILT_PEAK_CHR_NAME="chr[\dXY]+"
  MITO_CHR_NAME="chrM"
  REF_FA="https://storage.googleapis.com/encode-pipeline-genome-data/mm9/mm9.fa.gz"
  REF_MITO_FA="https://storage.googleapis.com/encode-pipeline-genome-data/mm9/mm9.chrM.fa.gz"
  CHRSZ="https://storage.googleapis.com/encode-pipeline-genome-data/mm9/mm9.chrom.sizes"
  BWT2_IDX="https://storage.googleapis.com/encode-pipeline-genome-data/mm9/bowtie2_index/mm9.fa.tar"
  BWT2_MITO_IDX="https://storage.googleapis.com/encode-pipeline-genome-data/mm9/bowtie2_index/mm9.chrM.fa.tar"
  BWA_IDX="https://storage.googleapis.com/encode-pipeline-genome-data/mm9/bwa_index/mm9.fa.tar"
  BWA_MITO_IDX="https://storage.googleapis.com/encode-pipeline-genome-data/mm9/bwa_index/mm9.chrM.fa.tar"
  BLACKLIST="https://storage.googleapis.com/encode-pipeline-genome-data/mm9/mm9-blacklist.bed.gz"

  TSS="https://storage.googleapis.com/encode-pipeline-genome-data/mm9/ataqc/mm9_gencode_tss_unique.bed.gz"
  DNASE="https://storage.googleapis.com/encode-pipeline-genome-data/mm9/ataqc/mm9_univ_dhs_ucsc.from_mm10.bed.gz"
  PROM="https://storage.googleapis.com/encode-pipeline-genome-data/mm9/ataqc/tss_mm9_master.from_mm10.bed.gz"
  ENH="https://storage.googleapis.com/encode-pipeline-genome-data/mm9/ataqc/mm9_enh_dhs_ucsc.from_mm10.bed.gz"
  REG2MAP_BED="https://storage.googleapis.com/encode-pipeline-genome-data/mm9/ataqc/mm9_dhs_universal_ucsc_v1.bed.gz"
  REG2MAP="https://storage.googleapis.com/encode-pipeline-genome-data/mm9/ataqc/dnase_avgs_merged_named.fseq.vals.gz"
  ROADMAP_META="https://storage.googleapis.com/encode-pipeline-genome-data/mm9/ataqc/accession_to_name.txt"

elif [[ "${GENOME}" == "hg38" ]]; then
  REGEX_BFILT_PEAK_CHR_NAME="chr[\dXY]+"
  MITO_CHR_NAME="chrM"
  REF_FA="https://www.encodeproject.org/files/GRCh38_no_alt_analysis_set_GCA_000001405.15/@@download/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta.gz"
  REF_MITO_FA="https://storage.googleapis.com/encode-pipeline-genome-data/hg38/GRCh38_no_alt_analysis_set_GCA_000001405.15.chrM.fa.gz"
  CHRSZ="https://storage.googleapis.com/encode-pipeline-genome-data/hg38/hg38.chrom.sizes"
  BWT2_IDX="https://storage.googleapis.com/encode-pipeline-genome-data/hg38/bowtie2_index/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta.tar"
  BWT2_MITO_IDX="https://storage.googleapis.com/encode-pipeline-genome-data/hg38/bowtie2_index/GRCh38_no_alt_analysis_set_GCA_000001405.15.chrM.fa.tar"
  BWA_IDX="https://storage.googleapis.com/encode-pipeline-genome-data/hg38/bwa_index/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta.tar"
  BWA_MITO_IDX="https://storage.googleapis.com/encode-pipeline-genome-data/hg38/bwa_index/GRCh38_no_alt_analysis_set_GCA_000001405.15.chrM.fa.tar"
  BLACKLIST="https://storage.googleapis.com/encode-pipeline-genome-data/hg38/hg38.blacklist.bed.gz"
  TSS="https://storage.googleapis.com/encode-pipeline-genome-data/hg38/ataqc/hg38_gencode_tss_unique.bed.gz"
  DNASE="https://storage.googleapis.com/encode-pipeline-genome-data/hg38/ataqc/reg2map_honeybadger2_dnase_all_p10_ucsc.hg19_to_hg38.bed.gz"
  PROM="https://storage.googleapis.com/encode-pipeline-genome-data/hg38/ataqc/reg2map_honeybadger2_dnase_prom_p2.hg19_to_hg38.bed.gz"
  ENH="https://storage.googleapis.com/encode-pipeline-genome-data/hg38/ataqc/reg2map_honeybadger2_dnase_enh_p2.hg19_to_hg38.bed.gz"
  REG2MAP_BED="https://storage.googleapis.com/encode-pipeline-genome-data/hg38/ataqc/hg38_celltype_compare_subsample.bed.gz"
  REG2MAP="https://storage.googleapis.com/encode-pipeline-genome-data/hg38/ataqc/hg38_dnase_avg_fseq_signal_formatted.txt.gz"
  ROADMAP_META="https://storage.googleapis.com/encode-pipeline-genome-data/hg38/ataqc/hg38_dnase_avg_fseq_signal_metadata.txt"
  if [[ "${VER}" == "v2" ]]; then
    BWT2_IDX="https://storage.googleapis.com/encode-pipeline-genome-data/hg38/bowtie2_index/ENCFF110MCL.tar.gz"
    BWA_IDX="https://storage.googleapis.com/encode-pipeline-genome-data/hg38/bwa_index/ENCFF643CGH.tar.gz"
    BLACKLIST="https://storage.googleapis.com/encode-pipeline-genome-data/hg38/ENCFF419RSJ.bed.gz"
    TSS="https://storage.googleapis.com/encode-pipeline-genome-data/hg38/ataqc/tss.pc.gencode.v29.bed.gz"
  elif [[ "${VER}" == "v3" ]]; then
    REF_MITO_FA="https://www.encodeproject.org/files/GRCh38_no_alt_analysis_set_GCA_000001405.15_mito_only/@@download/GRCh38_no_alt_analysis_set_GCA_000001405.15_mito_only.fasta.gz"
    BWT2_IDX="https://www.encodeproject.org/files/ENCFF110MCL/@@download/ENCFF110MCL.tar.gz"
    BWT2_MITO_IDX="https://www.encodeproject.org/files/GRCh38_no_alt_analysis_set_GCA_000001405.15_mito_only_bowtie2_index/@@download/GRCh38_no_alt_analysis_set_GCA_000001405.15_mito_only_bowtie2_index.tar.gz"
    BWA_IDX="https://www.encodeproject.org/files/ENCFF643CGH/@@download/ENCFF643CGH.tar.gz"
    BWA_MITO_IDX="https://www.encodeproject.org/files/GRCh38_no_alt_analysis_set_GCA_000001405.15_mito_only_bwa_index/@@download/GRCh38_no_alt_analysis_set_GCA_000001405.15_mito_only_bwa_index.tar.gz"
    CHRSZ="https://www.encodeproject.org/files/GRCh38_EBV.chrom.sizes/@@download/GRCh38_EBV.chrom.sizes.tsv"
    BLACKLIST="https://www.encodeproject.org/files/ENCFF356LFX/@@download/ENCFF356LFX.bed.gz"
    TSS="https://www.encodeproject.org/files/ENCFF493CCB/@@download/ENCFF493CCB.bed.gz"
    DNASE="https://www.encodeproject.org/files/ENCFF304XEX/@@download/ENCFF304XEX.bed.gz"
    PROM="https://www.encodeproject.org/files/ENCFF140XLU/@@download/ENCFF140XLU.bed.gz"
    ENH="https://www.encodeproject.org/files/ENCFF212UAV/@@download/ENCFF212UAV.bed.gz"
  fi

elif [[ "${GENOME}" == "mm10" ]]; then
  REGEX_BFILT_PEAK_CHR_NAME="chr[\dXY]+"
  MITO_CHR_NAME="chrM"
  REF_FA="https://www.encodeproject.org/files/mm10_no_alt_analysis_set_ENCODE/@@download/mm10_no_alt_analysis_set_ENCODE.fasta.gz"
  REF_MITO_FA="https://storage.googleapis.com/encode-pipeline-genome-data/mm10/mm10_no_alt_analysis_set_ENCODE.chrM.fa.gz"
  CHRSZ="https://storage.googleapis.com/encode-pipeline-genome-data/mm10/mm10.chrom.sizes"
  BWT2_IDX="https://storage.googleapis.com/encode-pipeline-genome-data/mm10/bowtie2_index/mm10_no_alt_analysis_set_ENCODE.fasta.tar"
  BWT2_MITO_IDX="https://storage.googleapis.com/encode-pipeline-genome-data/mm10/bowtie2_index/mm10_no_alt_analysis_set_ENCODE.chrM.fa.tar"
  BWA_IDX="https://storage.googleapis.com/encode-pipeline-genome-data/mm10/bwa_index/mm10_no_alt_analysis_set_ENCODE.fasta.tar"
  BWA_MITO_IDX="https://storage.googleapis.com/encode-pipeline-genome-data/mm10/bwa_index/mm10_no_alt_analysis_set_ENCODE.chrM.fa.tar"
  BLACKLIST="https://storage.googleapis.com/encode-pipeline-genome-data/mm10/mm10.blacklist.bed.gz"
  TSS="https://storage.googleapis.com/encode-pipeline-genome-data/mm10/ataqc/mm10_gencode_tss_unique.bed.gz"
  DNASE="https://storage.googleapis.com/encode-pipeline-genome-data/mm10/ataqc/mm10_univ_dhs_ucsc.bed.gz"
  PROM="https://storage.googleapis.com/encode-pipeline-genome-data/mm10/ataqc/tss_mm10_master.bed.gz"
  ENH="https://storage.googleapis.com/encode-pipeline-genome-data/mm10/ataqc/mm10_enh_dhs_ucsc.bed.gz"
  REG2MAP_BED="https://storage.googleapis.com/encode-pipeline-genome-data/mm10/ataqc/mm10_celltype_compare_subsample.bed.gz"
  REG2MAP="https://storage.googleapis.com/encode-pipeline-genome-data/mm10/ataqc/mm10_dnase_avg_fseq_signal_formatted.txt.gz"
  ROADMAP_META="https://storage.googleapis.com/encode-pipeline-genome-data/mm10/ataqc/mm10_dnase_avg_fseq_signal_metadata.txt"
  if [[ "${VER}" == "v2" ]]; then
    BWT2_IDX="https://storage.googleapis.com/encode-pipeline-genome-data/mm10/bowtie2_index/ENCFF309GLL.tar.gz"
    BWA_IDX="https://storage.googleapis.com/encode-pipeline-genome-data/mm10/bwa_index/ENCFF018NEO.tar.gz"
    BLACKLIST="https://storage.googleapis.com/encode-pipeline-genome-data/mm10/ENCFF547MET.bed.gz"
  elif [[ "${VER}" == "v3" ]]; then
    REF_MITO_FA="https://www.encodeproject.org/files/mm10_no_alt_analysis_set_ENCODE_mito_only/@@download/mm10_no_alt_analysis_set_ENCODE_mito_only.fasta.gz"
    BWT2_IDX="https://www.encodeproject.org/files/ENCFF309GLL/@@download/ENCFF309GLL.tar.gz"
    BWT2_MITO_IDX="https://www.encodeproject.org/files/mm10_no_alt_analysis_set_ENCODE_mito_only_bowtie2_index/@@download/mm10_no_alt_analysis_set_ENCODE_mito_only_bowtie2_index.tar.gz"
    BWA_IDX="https://www.encodeproject.org/files/ENCFF018NEO/@@download/ENCFF018NEO.tar.gz"
    BWA_MITO_IDX="https://www.encodeproject.org/files/mm10_no_alt_analysis_set_ENCODE_mito_only_bwa_index/@@download/mm10_no_alt_analysis_set_ENCODE_mito_only_bwa_index.tar.gz"
    CHRSZ="https://www.encodeproject.org/files/mm10_no_alt.chrom.sizes/@@download/mm10_no_alt.chrom.sizes.tsv"
    BLACKLIST="https://www.encodeproject.org/files/ENCFF547MET/@@download/ENCFF547MET.bed.gz"
    TSS="https://www.encodeproject.org/files/ENCFF498BEJ/@@download/ENCFF498BEJ.bed.gz"
    DNASE="https://www.encodeproject.org/files/ENCFF015KVI/@@download/ENCFF015KVI.bed.gz"
    PROM="https://www.encodeproject.org/files/ENCFF206BQS/@@download/ENCFF206BQS.bed.gz"
    ENH="https://www.encodeproject.org/files/ENCFF580RGZ/@@download/ENCFF580RGZ.bed.gz"
  fi

elif [[ "${GENOME}" == "hg38_chr19_chrM" ]]; then
  REGEX_BFILT_PEAK_CHR_NAME="chr[\dXY]+"
  MITO_CHR_NAME="chrM"
  REF_FA="https://storage.googleapis.com/encode-pipeline-genome-data/hg38_chr19_chrM/GRCh38_no_alt_analysis_set_GCA_000001405.15.chr19_chrM.fasta.gz"
  REF_MITO_FA="https://storage.googleapis.com/encode-pipeline-genome-data/hg38_chr19_chrM/GRCh38_no_alt_analysis_set_GCA_000001405.15.chr19_chrM.chrM.fa.gz"
  CHRSZ="https://storage.googleapis.com/encode-pipeline-genome-data/hg38_chr19_chrM/hg38_chr19_chrM.chrom.sizes"
  BWT2_IDX="https://storage.googleapis.com/encode-pipeline-genome-data/hg38_chr19_chrM/bowtie2_index/GRCh38_no_alt_analysis_set_GCA_000001405.15.chr19_chrM.fasta.tar"
  BWT2_MITO_IDX="https://storage.googleapis.com/encode-pipeline-genome-data/hg38_chr19_chrM/bowtie2_index/GRCh38_no_alt_analysis_set_GCA_000001405.15.chr19_chrM.chrM.fa.tar"
  BWA_IDX="https://storage.googleapis.com/encode-pipeline-genome-data/hg38_chr19_chrM/bwa_index/GRCh38_no_alt_analysis_set_GCA_000001405.15.chr19_chrM.fasta.tar"
  BWA_MITO_IDX="https://storage.googleapis.com/encode-pipeline-genome-data/hg38_chr19_chrM/bwa_index/GRCh38_no_alt_analysis_set_GCA_000001405.15.chr19_chrM.chrM.fa.tar"
  BLACKLIST="https://storage.googleapis.com/encode-pipeline-genome-data/hg38/hg38.blacklist.bed.gz"
  TSS="https://storage.googleapis.com/encode-pipeline-genome-data/hg38/ataqc/hg38_gencode_tss_unique.bed.gz"
  DNASE="https://storage.googleapis.com/encode-pipeline-genome-data/hg38/ataqc/reg2map_honeybadger2_dnase_all_p10_ucsc.hg19_to_hg38.bed.gz"
  PROM="https://storage.googleapis.com/encode-pipeline-genome-data/hg38/ataqc/reg2map_honeybadger2_dnase_prom_p2.hg19_to_hg38.bed.gz"
  ENH="https://storage.googleapis.com/encode-pipeline-genome-data/hg38/ataqc/reg2map_honeybadger2_dnase_enh_p2.hg19_to_hg38.bed.gz"
  REG2MAP_BED="https://storage.googleapis.com/encode-pipeline-genome-data/hg38/ataqc/hg38_celltype_compare_subsample.bed.gz"
  REG2MAP="https://storage.googleapis.com/encode-pipeline-genome-data/hg38/ataqc/hg38_dnase_avg_fseq_signal_formatted.txt.gz"
  ROADMAP_META="https://storage.googleapis.com/encode-pipeline-genome-data/hg38/ataqc/hg38_dnase_avg_fseq_signal_metadata.txt"
  if [[ "${VER}" == "v2" ]]; then
    BLACKLIST="https://storage.googleapis.com/encode-pipeline-genome-data/hg38/ENCFF419RSJ.bed.gz"
    TSS="https://storage.googleapis.com/encode-pipeline-genome-data/hg38/ataqc/tss.pc.gencode.v29.bed.gz"
  elif [[ "${VER}" == "v3" ]]; then
    BLACKLIST="https://www.encodeproject.org/files/ENCFF356LFX/@@download/ENCFF356LFX.bed.gz"
    TSS="https://www.encodeproject.org/files/ENCFF493CCB/@@download/ENCFF493CCB.bed.gz"
    DNASE="https://www.encodeproject.org/files/ENCFF304XEX/@@download/ENCFF304XEX.bed.gz"
    PROM="https://www.encodeproject.org/files/ENCFF140XLU/@@download/ENCFF140XLU.bed.gz"
    ENH="https://www.encodeproject.org/files/ENCFF212UAV/@@download/ENCFF212UAV.bed.gz"
  fi

elif [[ "${GENOME}" == "mm10_chr19_chrM" ]]; then
  REGEX_BFILT_PEAK_CHR_NAME="chr[\dXY]+"
  MITO_CHR_NAME="chrM"
  REF_FA="https://storage.googleapis.com/encode-pipeline-genome-data/mm10_chr19_chrM/mm10_no_alt_analysis_set_ENCODE.chr19_chrM.fasta.gz"
  REF_MITO_FA="https://storage.googleapis.com/encode-pipeline-genome-data/mm10_chr19_chrM/mm10_no_alt_analysis_set_ENCODE.chr19_chrM.chrM.fa.gz"
  CHRSZ="https://storage.googleapis.com/encode-pipeline-genome-data/mm10_chr19_chrM/mm10_chr19_chrM.chrom.sizes"
  BWT2_IDX="https://storage.googleapis.com/encode-pipeline-genome-data/mm10_chr19_chrM/bowtie2_index/mm10_no_alt_analysis_set_ENCODE.chr19_chrM.fasta.tar"
  BWT2_MITO_IDX="https://storage.googleapis.com/encode-pipeline-genome-data/mm10_chr19_chrM/bowtie2_index/mm10_no_alt_analysis_set_ENCODE.chr19_chrM.chrM.fa.tar"
  BWA_IDX="https://storage.googleapis.com/encode-pipeline-genome-data/mm10_chr19_chrM/bwa_index/mm10_no_alt_analysis_set_ENCODE.chr19_chrM.fasta.tar"
  BWA_MITO_IDX="https://storage.googleapis.com/encode-pipeline-genome-data/mm10_chr19_chrM/bwa_index/mm10_no_alt_analysis_set_ENCODE.chr19_chrM.chrM.fa.tar"
  BLACKLIST="https://storage.googleapis.com/encode-pipeline-genome-data/mm10/mm10.blacklist.bed.gz"
  TSS="https://storage.googleapis.com/encode-pipeline-genome-data/mm10/ataqc/mm10_gencode_tss_unique.bed.gz"
  DNASE="https://storage.googleapis.com/encode-pipeline-genome-data/mm10/ataqc/mm10_univ_dhs_ucsc.bed.gz"
  PROM="https://storage.googleapis.com/encode-pipeline-genome-data/mm10/ataqc/tss_mm10_master.bed.gz"
  ENH="https://storage.googleapis.com/encode-pipeline-genome-data/mm10/ataqc/mm10_enh_dhs_ucsc.bed.gz"
  REG2MAP_BED="https://storage.googleapis.com/encode-pipeline-genome-data/mm10/ataqc/mm10_celltype_compare_subsample.bed.gz"
  REG2MAP="https://storage.googleapis.com/encode-pipeline-genome-data/mm10/ataqc/mm10_dnase_avg_fseq_signal_formatted.txt.gz"
  ROADMAP_META="https://storage.googleapis.com/encode-pipeline-genome-data/mm10/ataqc/mm10_dnase_avg_fseq_signal_metadata.txt"
  if [[ "${VER}" == "v2" ]]; then
    BLACKLIST="https://storage.googleapis.com/encode-pipeline-genome-data/mm10/ENCFF547MET.bed.gz"
  elif [[ "${VER}" == "v3" ]]; then
    BLACKLIST="https://www.encodeproject.org/files/ENCFF547MET/@@download/ENCFF547MET.bed.gz"
    TSS="https://www.encodeproject.org/files/ENCFF498BEJ/@@download/ENCFF498BEJ.bed.gz"
    DNASE="https://www.encodeproject.org/files/ENCFF015KVI/@@download/ENCFF015KVI.bed.gz"
    PROM="https://www.encodeproject.org/files/ENCFF206BQS/@@download/ENCFF206BQS.bed.gz"
    ENH="https://www.encodeproject.org/files/ENCFF580RGZ/@@download/ENCFF580RGZ.bed.gz"
  fi
fi

if [[ -z "${REF_FA}" ]]; then
  echo "Error: unsupported genome $GENOME"
  exit 1
fi
if [[ -z "${MITO_CHR_NAME}" ]]; then
  echo "Error: Mitochondrial chromosome name must be defined"
  exit 1
fi
if [[ -z "${REGEX_BFILT_PEAK_CHR_NAME}" ]]; then
  echo "Error: Perl style reg-ex for filtering peaks must be defined"
  exit 1
fi


echo "=== Downloading files..."
wget -c -O $(basename ${REF_FA}) ${REF_FA}
wget -c -O $(basename ${REF_MITO_FA}) ${REF_MITO_FA}
if [[ ! -z "${BLACKLIST}" ]]; then wget -N -c ${BLACKLIST}; fi

echo "=== Processing reference fasta file..."
REF_FA_PREFIX=$(basename ${REF_FA})
REF_MITO_FA_PREFIX=$(basename ${REF_MITO_FA})

echo "=== Generating fasta index and chrom.sizes file..."
cd ${DEST_DIR}
wget -c -O $(basename ${CHRSZ}) ${CHRSZ}
CHRSZ=$(basename ${CHRSZ})

echo "=== Determinig gensz..."
cd ${DEST_DIR}
GENSZ=$(cat ${CHRSZ} | awk '{sum+=$2} END{print sum}')
if [[ "${GENOME}" == hg* ]]; then GENSZ=hs; fi
if [[ "${GENOME}" == mm* ]]; then GENSZ=mm; fi

echo "=== Downloading bowtie2 index..."
mkdir -p ${DEST_DIR}/bowtie2_index
cd ${DEST_DIR}/bowtie2_index
wget -c ${BWT2_IDX}
wget -c ${BWT2_MITO_IDX}

echo "=== Downloading bwa index..."
mkdir -p ${DEST_DIR}/bwa_index
cd ${DEST_DIR}/bwa_index
wget -c ${BWA_IDX}
wget -c ${BWA_MITO_IDX}

echo "=== Creating TSV file... (${TSV})"
cd ${DEST_DIR}
rm -f ${TSV}
touch ${TSV}

echo -e "genome_name\t${GENOME}" >> ${TSV}
echo -e "ref_fa\t${DEST_DIR}/$(basename $REF_FA_PREFIX)" >> ${TSV}
echo -e "ref_mito_fa\t${DEST_DIR}/$(basename $REF_MITO_FA_PREFIX)" >> ${TSV}
echo -e "mito_chr_name\t${MITO_CHR_NAME}" >> ${TSV}
printf "regex_bfilt_peak_chr_name\t%s\n" "${REGEX_BFILT_PEAK_CHR_NAME}" >> ${TSV}
echo -e "blacklist\t${DEST_DIR}/$(basename ${BLACKLIST})" >> ${TSV};
echo -e "chrsz\t${DEST_DIR}/$(basename ${CHRSZ})" >> ${TSV}
echo -e "gensz\t${GENSZ}" >> ${TSV}
echo -e "bowtie2_idx_tar\t${DEST_DIR}/bowtie2_index/$(basename $BWT2_IDX)" >> ${TSV}
echo -e "bowtie2_mito_idx_tar\t${DEST_DIR}/bowtie2_index/$(basename $BWT2_MITO_IDX)" >> ${TSV}
echo -e "bwa_idx_tar\t${DEST_DIR}/bwa_index/$(basename $BWA_IDX)" >> ${TSV}
echo -e "bwa_mito_idx_tar\t${DEST_DIR}/bwa_index/$(basename $BWA_MITO_IDX)" >> ${TSV}

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
