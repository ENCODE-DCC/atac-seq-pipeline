#!/usr/bin/env python

# ENCODE DCC fastq merger wrapper
# Author: Jin Lee (leepc12@gmail.com)

import sys
import os
import argparse
import multiprocessing
from encode_dcc_common import *

def parse_arguments(debug=False):
    parser = argparse.ArgumentParser(prog='ENCODE DCC fastq merger.',
                                        description='')
    parser.add_argument('--fastqs-R1', nargs='+', type=str, required=True,
                        help='FASTQs to be merged for R1. \
                            FASTQs must be compressed with gzip (with .gz).')
    parser.add_argument('--fastqs-R2', nargs='*', type=str,
                        help='FASTQs to be merged for R2. \
                            FASTQs must be compressed with gzip (with .gz).')
    parser.add_argument('--paired-end', action="store_true",
                        help='Paired-end FASTQs.')
    parser.add_argument('--nth', type=int, default=1,
                        help='Number of threads to parallelize.')
    parser.add_argument('--out-dir', default='', type=str,
                            help='Output directory.')
    parser.add_argument('--out-txt', default='out.txt', type=str,
                            help='Single-column text file with all output filenames \
                            ([OUT_DIR]/[BASENAME]). This file will be generated \
                            under [OUT_DIR].')
    parser.add_argument('--log-level', default='INFO', 
                        choices=['NOTSET','DEBUG','INFO',
                            'WARNING','CRITICAL','ERROR','CRITICAL'],
                        help='Log level')
    args = parser.parse_args()

    if args.paired_end:
        if len(args.fastqs_R2)==0:
            raise argparse.ArgumentTypeError(
            'Define --fastqs-R2 for paired end dataset.')
        if len(args.fastqs_R1)!=len(args.fastqs_R2):
            raise argparse.ArgumentTypeError(
            'Dimension mismatch between --fastqs-R2 and --fastqs-R2.')
    else:
        if len(args.fastqs_R2)>0:
            raise argparse.ArgumentTypeError(
            'Do not define --fastqs-R2 for single ended dataset.')
    log.setLevel(args.log_level)
    log.info(sys.argv)
    return args

def merge_fastqs(fastqs, out_dir):
    if len(fastqs)>1:
        prefix = os.path.join(out_dir,
            os.path.basename(strip_ext_fastq(fastqs[0])))
        merged = '{}.merged.fastq.gz'.format(prefix)
        cmd = 'zcat -f {} | gzip -nc > {}'.format(
            ' '.join(fastqs),
            merged)
        run_shell_cmd(cmd)
        return merged
    else:
        return make_hard_link(fastqs[0], out_dir)

def main():
    # read params
    args = parse_arguments()
    log.info('Initializing and making output directory...')

    # make out_dir
    mkdir_p(args.out_dir)

    # initialize multithreading
    log.info('Initializing multithreading...')
    if args.paired_end:
        num_process = min(2*len(args.fastqs_R1),args.nth)
    else:
        num_process = min(len(args.fastqs_R1),args.nth)
    log.info('Number of threads={}.'.format(num_process))
    pool = multiprocessing.Pool(num_process)

    # merge fastqs
    log.info('Merging fastqs...')
    log.info('R1 to be merged: {}'.format(args.fastqs_R1))
    ret_val1 = pool.apply_async(merge_fastqs,
                    (args.fastqs_R1, args.out_dir,))
    R1_merged = ret_val1.get(BIG_INT)
    if args.paired_end:
        log.info('R2 to be merged: {}'.format(args.fastqs_R2))
        ret_val2 = pool.apply_async(merge_fastqs,
                        (args.fastqs_R2, args.out_dir,))
        R2_merged = ret_val2.get(BIG_INT)

    # close multithreading
    pool.close()
    pool.join()

    # write txt for output filenames (for WDL)
    log.info('Writinig output txt...')    
    txt = os.path.join(args.out_dir, args.out_txt)
    if args.paired_end:
        arr = [R1_merged, R2_merged]
    else:
        arr = [R1_merged]
    write_txt(txt, arr)

    log.info('List all files in output directory...')
    ls_l(args.out_dir)

    log.info('All done.')

if __name__=='__main__':
    main()