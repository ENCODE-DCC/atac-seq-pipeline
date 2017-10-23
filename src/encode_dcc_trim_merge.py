#!/usr/bin/env python

# ENCODE DCC adapter trimmer and fastq merger
# Author: Jin Lee (leepc12@gmail.com)

import sys
import os
import re
import json
import logging
import gzip
import collections
import argparse
import multiprocessing
import subprocess
import copy

logger = logging.getLogger(__name__)

def parse_arguments():
    parser = argparse.ArgumentParser(prog='ENCODE DCC adapter trimmer and fastq merger',
                                        description='')
    parser.add_argument('fastqs', nargs='+', type=str,
                        help='TSV file path or \
                            list of FASTQ files delimited by space. Use \
                            TSV for multiple techincal replicates).')
    parser.add_argument('--nth', type=int, default=1,
                        help='No. threads to parallelize bowtie2.')
    parser.add_argument('--paired-end', action="store_true",
                        help='FASTQs have paired-end')
    group_out = parser.add_mutually_exclusive_group()
    group_out.add_argument('--out-dir', default='.', type=str,
                            help='Output directory. Prefix will be taken from FASTQs.')
    group_out.add_argument('--out-prefix', type=str,
                            help='Output prefix path.')
    parser.add_argument('--adapters', nargs='+', type=str,
                        help='TSV file path or \
                            list of adapters delimited by space (\
                            for multiple techincal replicates use TSV).')
    parser.add_argument('--auto-detect-adapter', action='store_true',
                        help='Automatically detect/trim adapters \
                            (supported system: Illumina, Nextera and smallRNA).')
    parser.add_argument('--min-trim-len', type=int, default=5,
                        help='Minimum trim length for cutadapt -m \
                            (throwing away processed reads shorter than this).')
    parser.add_argument('--err-rate', type=float, default=0.1,
                        help='Maximum allowed adapter error rate for cutadapt -e \
                            (no. errors divided by the length \
                            of the matching adapter region).')
    args = parser.parse_args()
    # parse fastqs command line
    if args.fastqs[0].endswith('.gz'): # it's fastq
        args.fastqs = [args.fastqs]
    else: # it's TSV
        args.fastqs = read_tsv(args.fastqs[0])
    # parse --adapters command line
    if args.adapters:
        if os.path.exists(args.adapters[0]): # it's TSV
            args.adapters = read_tsv(args.adapters[0])
        else:
            args.adapters = [args.adapters]
        if args.adapters and len(args.adapters) != len(args.fastqs):
            raise ValueError(
                'fastqs and adapters must have the same dimension.')
    else: # fill empty string in adapter list
        args.adapters = copy.deepcopy(args.fastqs)
        for i in range(len(args.adapters)):
            for j in range(len(args.adapters[i])):
                args.adapters[i][j] = ''

    return args

def detect_adapter(fastq):
    adapter = ''
    return adapter

def trim_adapter_se(fastq, adapter, min_trim_len, err_rate, out_dir):
    prefix = os.path.join(out_dir,
        os.path.basename(strip_ext_fastq(fastq)))
    trimmed = '{}.trim.fastq.gz'.format(prefix)
    cmd = 'cutadapt {} -e {} -a {} {} | gzip -nc > {}'.format(
        '-m {}'.format(min_trim_len) if min_trim_len > 0 else '',
        err_rate,
        adapter,
        fastq,
        trimmed)     
    print(cmd)
    print(subprocess.check_output(cmd, shell=True))
    return trimmed

def trim_adapter_pe(fastq1, fastq2, adapter1, adapter2,
        min_trim_len, err_rate, out_dir):
    prefix1 = os.path.join(out_dir,
        os.path.basename(strip_ext_fastq(fastq1)))
    prefix2 = os.path.join(out_dir,
        os.path.basename(strip_ext_fastq(fastq2)))
    trimmed1 = '{}.trim.fastq.gz'.format(prefix1)
    trimmed2 = '{}.trim.fastq.gz'.format(prefix2)
    cmd = 'cutadapt {} -e {} -a {} -A {} {} {} -o {} -p {}'.format(
        '-m {}'.format(min_trim_len) if min_trim_len > 0 else '',
        err_rate,
        adapter1, adapter2,
        fastq1, fastq2,
        trimmed1, trimmed2)
    print(cmd)
    print(subprocess.check_output(cmd, shell=True))
    return [trimmed1, trimmed2]

def merge_fastqs(fastqs, out_dir):
    prefix1 = os.path.join(out_dir,
        os.path.basename(strip_ext_fastq(fastqs[0])))
    merged = '{}.merged.fastq.gz'.format(prefix)
    cmd = 'zcat {} | gzip -nc > {}'.format(
        ' '.join(fastqs),
        merged)
    print(cmd)
    print(subprocess.check_output(cmd, shell=True))
    return merged

def main():
    args = parse_arguments()
    temp_files = []
    R1_to_merge = []
    R2_to_merge = []
        
    # detect adapters (with multithreading)
    pool_detect = multiprocessing.Pool(args.nth) 
    ret_vals = {}
    for i in range(len(args.fastqs)):
        # for each technical replicate
        fastqs = args.fastqs[i] # R1 and R2
        adapters = args.adapters[i]
        if args.paired_end:
            if args.auto_detect_adapter and \
                not (adapters[0] and adapters[1]):
                ret_val1 = pool_detect.apply_async(detect_adapter,
                    (fastqs[0]))
                ret_val2 = pool_detect.apply_async(detect_adapter,
                    (fastqs[1]))
                ret_vals[i]=[ret_val1,ret_val2]
        else:
            if args.auto_detect_adapter and \
                not adapters[0]:
                ret_val1 = pool_detect.apply_async(detect_adapter,
                    (fastqs[0]))
                ret_vals[i]=[ret_val1]
    # update array with detected adapters
    for i in ret_vals:
        for j, ret_vals[i]:
            args.adapters[i][j] = ret_vals[i][j].get()
    pool_detect.close()
    pool_detect.join()

    # trim adapters
    pool_detect = multiprocessing.Pool(args.nth)
    adapters[tech_rep][0,1]
    
    # merge fastqs
    for i in range(len(args.fastqs)):
        # for each technical replicate
        fastqs = args.fastqs[i] # R1 and R2
        adapters = args.adapters[i] if args.adapters else [None,None]
        if args.paired_end:
            if not args.auto_detect_adapter and \
                not (adapters[0] and adapters[1]):
                trimmed_fastqs = fastqs
            else:
                if args.auto_detect_adapter and \
                    not (adapters[0] and adapters[1]): # detect adapters
                    adapters[0] = detect_adapter(fastqs[0])
                    adapters[1] = detect_adapter(fastqs[1])
                if adapters[0] and adapters[1]:
                    trimmed_fastqs = trim_adapter_pe(
                        fastqs[0], fastqs[1], 
                        adapters[1], adapter[2],
                        args.min_trim_len,
                        args.err_rate,
                        args.out_dir)
                    temp_files.append(trimmed_fastqs)
                else:
                    trimmed_fastqs = fastqs
            R1_to_merge.append(trimmed_fastqs[0])
            R2_to_merge.append(trimmed_fastqs[1])
        else:

            if not args.auto_detect_adapter and \
                not adapters[0]:
                trimmed_fastq = fastq[0]
            else:
                if args.auto_detect_adapter and not adapters[0]:
                    adapters[0] = detect_adapter(fastqs[0])
                if adapters[0]:
                    trimmed_fastqs = trim_adapter_se(
                        fastqs[0],
                        adapters[0],
                        args.min_trim_len,
                        args.err_rate,
                        args.out_dir)
                    temp_files.append(trimmed_fastqs)
                else:
                    trimmed_fastq = fastqs[0]
            R1_to_merge.append(trimmed_fastq)
    # merge fastqs
    R1_merged = merge_fastqs(R1_to_merge)
    R2_merged = merge_fastqs(R2_to_merge)
    # remove temporary/intermediate files
    subprocess.check_output(shlex.split(
        'rm -rf {}'.format(' '.join(temp_files))))    

if __name__=='__main__':
    main()