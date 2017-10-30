#!/usr/bin/env python

# ENCODE DCC adapter trimmer and fastq merger
# Author: Jin Lee (leepc12@gmail.com)

import sys
import os
import argparse
import multiprocessing
import copy
from detect_adapter import detect_most_likely_adapter
from encode_dcc_common import *

def parse_arguments(debug=False):
    parser = argparse.ArgumentParser(prog='ENCODE DCC adapter trimmer and fastq merger',
                                        description='')
    parser.add_argument('fastqs', nargs='+', type=str,
                        help='TSV file path or \
                            list of FASTQs delimited by space. \
                            FASTQs must be compressed with gzip (with .gz). \
                            Use TSV for multiple techincal replicates).')
    parser.add_argument('--nth', type=int, default=1,
                        help='Number of threads to parallelize.')
    parser.add_argument('--paired-end', action="store_true",
                        help='Paired-end FASTQs.')
    parser.add_argument('--out-dir', default='.', type=str,
                            help='Output directory.')
    parser.add_argument('--out-meta-tsv', default='meta.tsv', type=str,
                            help='Basename for metadata TSV file with all output files.')
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
    parser.add_argument('--log-level', default='INFO', 
                        choices=['NOTSET','DEBUG','INFO',
                            'WARNING','CRITICAL','ERROR','CRITICAL'],
                        help='Log level')
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
            if not args.adapters[0]: # if empty TSV
                args.adapters = None
        else:
            args.adapters = [args.adapters]        
    if not args.adapters: # fill empty string in adapter list
        args.adapters = copy.deepcopy(args.fastqs)
        for i, adapters in enumerate(args.adapters):
            for j, adapter in enumerate(adapters):
                args.adapters[i][j] = ''

    # check if fastqs, adapers have correct dimension
    if len(args.adapters)!=len(args.fastqs):
        raise ValueError('fastqs and adapters dimension mismatch.')
    for i, fastqs in enumerate(args.fastqs):
        if args.paired_end and len(fastqs)!=2:
            raise ValueError('Need 2 fastqs per replicate for paired end.')
        if not args.paired_end and len(fastqs)!=1:
            raise ValueError('Need 1 fastq per replicate for single end.')
        if len(fastqs)!=len(args.adapters[i]):
            raise ValueError('fastqs and adapters dimension mismatch.')
            
    log.setLevel(args.log_level)
    log.info(sys.argv)
    return args

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
    run_shell_cmd(cmd)
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
    run_shell_cmd(cmd)
    return [trimmed1, trimmed2]

def merge_fastqs(fastqs, out_dir):
    prefix = os.path.join(out_dir,
        os.path.basename(strip_ext_fastq(fastqs[0])))
    merged = '{}.merged.fastq.gz'.format(prefix)

    if len(fastqs)>1:
        cmd = 'zcat -f {} | gzip -nc > {}'.format(
            ' '.join(fastqs),
            merged)
    else:
        cmd = 'cp {} {}'.format(
            fastqs[0],
            merged)
    run_shell_cmd(cmd)
    return merged

def main():
    # read params
    args = parse_arguments()
    log.info('Initializing and making output directory...')

    # make out_dir
    mkdir_p(args.out_dir)

    # declare temp arrays
    temp_files = [] # files to deleted later at the end

    # initialize multithreading
    log.info('Initializing multithreading...')
    if args.paired_end:
        num_process = min(2*len(args.fastqs),args.nth)
    else:
        num_process = min(len(args.fastqs),args.nth)
    log.info('Number of threads={}.'.format(num_process))
    pool = multiprocessing.Pool(num_process)

    # detect adapters
    log.info('Detecting adapters...')
    ret_vals = {}
    for i in range(len(args.fastqs)):
        # for each technical replicate
        log.info('Detecting adapters for tech_rep{}...'.format(
                i+1))
        fastqs = args.fastqs[i] # R1 and R2
        adapters = args.adapters[i]
        if args.paired_end:
            if args.auto_detect_adapter and \
                not (adapters[0] and adapters[1]):                
                ret_val1 = pool.apply_async(
                    detect_most_likely_adapter,
                        (fastqs[0],))
                ret_val2 = pool.apply_async(
                    detect_most_likely_adapter,
                        (fastqs[1],))
                ret_vals[i]=[ret_val1,ret_val2]
        else:
            if args.auto_detect_adapter and \
                not adapters[0]:
                ret_val1 = pool.apply_async(
                    detect_most_likely_adapter,
                        (fastqs[0],))
                ret_vals[i]=[ret_val1]

    # update array with detected adapters
    for i in ret_vals:
        for j in range(len(ret_vals[i])):
            args.adapters[i][j] = str(ret_vals[i][j].get(BIG_INT))
            log.info('Detected adapters for tech_rep{}, R{}: {}'.format(
                    i+1, j+1, args.adapters[i][j]))

    # trim adapters
    log.info('Trimming adapters...')
    trimmed_fastqs = copy.deepcopy(args.fastqs)
    ret_vals = {}
    for i in range(len(args.fastqs)):
        # for each technical replicate
        fastqs = args.fastqs[i] # R1 and R2
        adapters = args.adapters[i]
        if args.paired_end:
            if adapters[0] and adapters[1]:
                ret_vals[i] = pool.apply_async(
                    trim_adapter_pe,(
                        fastqs[0], fastqs[1], 
                        adapters[0], adapters[1],
                        args.min_trim_len,
                        args.err_rate,
                        args.out_dir))
        else:
            if adapters[0]:
                ret_vals[i] = pool.apply_async(
                    trim_adapter_se,(
                        fastqs[0],
                        adapters[0],
                        args.min_trim_len,
                        args.err_rate,
                        args.out_dir))

    # update array with trimmed fastqs
    for i in ret_vals:
        if args.paired_end:
            trimmed_fastqs[i] = ret_vals[i].get(BIG_INT)
        else:
            trimmed_fastqs[i] = [ret_vals[i].get(BIG_INT)]
        temp_files.extend(trimmed_fastqs[i])

    # merge fastqs
    log.info('Merging fastqs...')
    R1_to_merge = []
    R2_to_merge = []
    for i in range(len(trimmed_fastqs)):
        R1_to_merge.append(trimmed_fastqs[i][0])
        if args.paired_end:
            R2_to_merge.append(trimmed_fastqs[i][1])
    log.info('R1: {}'.format(R1_to_merge))
    ret_val1 = pool.apply_async(merge_fastqs,
                    (R1_to_merge,args.out_dir,))
    R1_merged = ret_val1.get(BIG_INT)
    if args.paired_end:
        log.info('R2: {}'.format(R2_to_merge))
        ret_val2 = pool.apply_async(merge_fastqs,
                        (R2_to_merge,args.out_dir,))
        R2_merged = ret_val2.get(BIG_INT)

    # close multithreading
    pool.close()
    pool.join()

    # write metadata TSV for outputs (for WDL)
    tsv = os.path.join(args.out_dir, args.out_meta_tsv)
    if args.paired_end:
        arr = [[R1_merged, R2_merged]]
    else:
        arr = [[R1_merged]]
    write_tsv(tsv, arr)

    # remove temporary/intermediate files
    log.info('Removing temporary files...')
    rm_f(temp_files)

    log.info('All done.')

if __name__=='__main__':
    main()