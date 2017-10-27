#!/usr/bin/env python

# ENCODE DCC filter python script
# Author: Jin Lee (leepc12@gmail.com)

import sys
import os
import re
import json
import collections
import argparse
import multiprocessing
import subprocess
import copy
import encode_dcc_common

def parse_arguments():
    parser = argparse.ArgumentParser(prog='ENCODE DCC filter python script',
                                        description='')
    parser.add_argument('bam', type=str,
                        help='Path for BAM file.')
    parser.add_argument('--dup-marker', type=str, choices=['picard','sambamba'],
                        default='picard',
                        help='Dup marker for filtering mapped reads in BAM.')    
    parser.add_argument('--nth', type=int, default=1,
                        help='No. threads to parallelize.')
    parser.add_argument('--paired-end', action="store_true",
                        help='Paired-end BAM')
    parser.add_argument('--out-dir', default='.', type=str,
                            help='Output directory. Prefix will be taken from BAM.')
    parser.add_argument('--multimapping', default=0, type=int,
                        help='Multimapping.')
    parser.add_argument('--mapq-thresh', default=30, type=int,
                        help='Threshold for low MAPQ reads removal.')
    parser.add_argument('--no-dup-removal', action="store_true",
                        help='No dupe reads removal when filtering raw BAM')
    args = parser.parse_args()
    return args

def rm_unmapped_lowq_reads_se(bam, nth, out_dir):
    prefix = os.path.join(out_dir,
        os.path.basename(strip_ext_bam(bam)))
    filt_bam = '{}.filt.bam'.format(prefix)

    if [[ $multimapping > 0 ]]; then \
        sambamba sort -t $nth_dedup $bam -n -o $qname_sort_bam; \
        samtools view -h $qname_sort_bam | $(which assign_multimappers.py) -k $multimapping | \
        samtools view -F 1804 -Su /dev/stdin | \
        sambamba sort -t $nth_dedup /dev/stdin -o $filt_bam; \
        rm -f $qname_sort_bam; \
        else \
            samtools view -F 1804 -q $mapq_thresh -u $bam | \
            sambamba sort -t $nth_dedup /dev/stdin -o $filt_bam; \
        fi

def rm_unmapped_lowq_reads_pe(bam, nth, out_dir):
    prefix = os.path.join(out_dir,
        os.path.basename(strip_ext_bam(bam)))
    filt_bam = '{}.filt.bam'.format(prefix)

    if [[ $multimapping > 0 ]]; then \
        samtools view -F 524 -f 2 -u $bam | \
        sambamba sort -t $nth_dedup -n /dev/stdin -o $tmp_filt_bam; \
        samtools view -h $tmp_filt_bam | \
        $(which assign_multimappers.py) -k $multimapping --paired-end | \
        samtools fixmate -r /dev/stdin $tmp_filt_bam.fixmate.bam; \
        else \
            samtools view -F 1804 -f 2 -q $mapq_thresh -u $bam | \
            sambamba sort -t $nth_dedup -n /dev/stdin -o $tmp_filt_bam; \
            samtools fixmate -r $tmp_filt_bam $tmp_filt_bam.fixmate.bam; \
        fi
    sys samtools view -F 1804 -f 2 -u $tmp_filt_bam.fixmate.bam | sambamba sort -t $nth_dedup /dev/stdin -o $filt_bam

    sys rm -f $tmp_filt_bam.fixmate.bam
    sys rm -f $tmp_filt_bam

def mark_dup_picard(filt_bam, nth, out_dir): # shared by both se and pe
    if [ -f "${MARKDUP}" ]; then \
            java -Xmx4G -jar ${MARKDUP} \
                INPUT="$bam" OUTPUT="$dupmark_bam" \
                METRICS_FILE="$dup_qc" VALIDATION_STRINGENCY=LENIENT \
                ASSUME_SORTED=true REMOVE_DUPLICATES=false; \
            else \
            picard MarkDuplicates \
                INPUT="$bam" OUTPUT="$dupmark_bam" \
                METRICS_FILE="$dup_qc" VALIDATION_STRINGENCY=LENIENT \
                ASSUME_SORTED=true REMOVE_DUPLICATES=false; \
            fi
def mark_dup_sambamba(filt_bam, nth, out_dir): # shared by both se and pe
    sambamba markdup -t $nth_markdup --hash-table-size=17592186044416 \
            --overflow-list-size=20000000 --io-buffer-size=256 $bam $dupmark_bam \
            2> $dup_qc

def rm_dup_se(markdup_bam, out_dir):
    sys samtools view -F 1804 -b $dupmark_bam > $nodup_bam

                //# Index Final BAM file
                sys sambamba index -t $nth_dedup $nodup_bam

                sys sambamba flagstat -t $nth_dedup $nodup_bam > $map_qc

                //# =============================
                //# Compute library complexity
                //# =============================
                //# sort by position and strand
                //# Obtain unique count statistics

                //# PBC File output
                //# TotalReadPairs [tab] DistinctReadPairs [tab] OneReadPair [tab] TwoReadPairs [tab] NRF=Distinct/Total [tab] PBC1=OnePair/Distinct [tab] PBC2=OnePair/TwoPair
                sys bedtools bamtobed -i $dupmark_bam | \
                    awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$6}' | \
                    grep -v 'chrM' | sort | uniq -c | \
                    awk 'BEGIN{mt=0;m0=0;m1=0;m2=0} ($1==1){m1=m1+1} ($1==2){m2=m2+1} {m0=m0+1} {mt=mt+$1} END{m1_m2=-1.0; if(m2>0) m1_m2=m1/m2; printf "%d\t%d\t%d\t%d\t%f\t%f\t%f\n",mt,m0,m1,m2,m0/mt,m1/m0,m1_m2}' > $pbc_qc
                sys $shcmd_finalize

def rm_dup_pe(markdup_bam, out_dir):
    sys samtools view -F 1804 -f 2 -b $dupmark_bam > $nodup_bam

    sys sambamba index -t $nth_dedup $nodup_bam

    sys sambamba flagstat -t $nth_dedup $nodup_bam > $map_qc

    //# =============================
    //# Compute library complexity
    //# =============================
    //# Sort by name
    //# convert to bedPE and obtain fragment coordinates
    //# sort by position and strand
    //# Obtain unique count statistics

    //# TotalReadPairs [tab] DistinctReadPairs [tab] OneReadPair [tab] TwoReadPairs [tab] NRF=Distinct/Total [tab] PBC1=OnePair/Distinct [tab] PBC2=OnePair/TwoPair
    sys sambamba sort -t $nth_dedup -n $dupmark_bam -o $dupmark_bam.tmp.bam

    sys bedtools bamtobed -bedpe -i $dupmark_bam.tmp.bam | \
        awk 'BEGIN{OFS="\t"}{print $1,$2,$4,$6,$9,$10}' | \
        grep -v 'chrM' | sort | uniq -c | \
        awk 'BEGIN{mt=0;m0=0;m1=0;m2=0} ($1==1){m1=m1+1} ($1==2){m2=m2+1} {m0=m0+1} {mt=mt+$1} END{m1_m2=-1.0; if(m2>0) m1_m2=m1/m2; printf "%d\t%d\t%d\t%d\t%f\t%f\t%f\n",mt,m0,m1,m2,m0/mt,m1/m0,m1_m2}' > $pbc_qc
    sys rm -f $dupmark_bam.tmp.bam

def main():
    # read params
    args = parse_arguments()
    log.info('Initializing and making output directory...')

    # make out_dir
    mkdir_p(args.out_dir)

    # declare temp arrays
    temp_files = [] # files to deleted later at the end

    # generate read length file
    log.info('Removing unmapped/low-quality reads...')
    if args.paired_end:
        filt_bam = make_read_length_file(
                            args.fastqs[1], args.out_dir)
    else:
        filt_bam = 
    
    # if bowtie2 index is tarball then unpack it
    if args.bowtie2_index_prefix_or_tar.endswith('.tar'):
        log.info('Unpacking bowtie2 index tar...')
        tar = args.bowtie2_index_prefix_or_tar
        # untar
        untar(tar, args.out_dir)
        bowtie2_index_prefix = os.path.join(
            args.out_dir, os.path.basename(strip_ext_tar(tar)))
        temp_files.append('{}.*.bt2'.format(
            bowtie2_index_prefix))
    else:
        bowtie2_index_prefix = args.bowtie2_index_prefix_or_tar

    # check if bowties indices are unpacked on out_dir
    chk_bowtie2_index(bowtie2_index_prefix)

    # bowtie2
    log.info('Running bowtie2...')
    if args.paired_end:
        bam, align_log = bowtie2_pe(
            args.fastqs[0], args.fastqs[1], 
            bowtie2_index_prefix,
            args.multimapping, args.score_min, args.nth,
            args.out_dir)
    else:
        bam, align_log = bowtie2_se(
            args.fastqs[0], 
            bowtie2_index_prefix,
            args.multimapping, args.score_min, args.nth,
            args.out_dir)

    # initialize multithreading
    log.info('Initializing multithreading...')
    num_process = min(2,args.nth)
    log.info('Number of threads={}.'.format(num_process))
    pool = multiprocessing.Pool(num_process)
    
    # samtools index
    log.info('Running samtools index...')
    ret_val1 = pool.apply_async(
        samtools_index, (bam, args.out_dir))

    # samtools flagstat qc
    log.info('Running samtools flagstat...')
    ret_val2 = pool.apply_async(
        samtools_flagstat, (bam, args.out_dir))

    bai = ret_val1.get()
    flagstat_qc = ret_val2.get()

    # close multithreading
    pool.close()
    pool.join()

    # remove temporary/intermediate files
    log.info('Removing temporary files...')
    run_shell_cmd('rm -rf {}'.format(' '.join(temp_files)))

    log.info('All done.')

if __name__=='__main__':
    main()