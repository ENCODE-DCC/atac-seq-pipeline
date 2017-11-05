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
    parser.add_argument('fastqs', nargs='+', type=str,
                        help='TSV file path or list of FASTQs. \
                            FASTQs must be compressed with gzip (with .gz). \
                            Use TSV for paired end fastqs or \
                            multiple techincal replicates \
                            row=tech_rep_id, col=pair_id). \
                            All rows will be merged.')
    parser.add_argument('--nth', type=int, default=1,
                        help='Number of threads to parallelize.')
    parser.add_argument('--out-dir', default='.', type=str,
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

    # parse fastqs command line
    if args.fastqs[0].endswith('.gz'): # it's single-end fastq
        args.fastqs = [[f] for f in args.fastqs] # make it a matrix
    else: # it's TSV
        args.fastqs = read_tsv(args.fastqs[0])

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

    paired_end = len(args.fastqs[0])>1
    # initialize multithreading
    log.info('Initializing multithreading...')
    if paired_end:
        num_process = min(2*len(args.fastqs),args.nth)
    else:
        num_process = min(len(args.fastqs),args.nth)
    log.info('Number of threads={}.'.format(num_process))
    pool = multiprocessing.Pool(num_process)

    # merge fastqs
    log.info('Merging fastqs...')
    R1_to_merge = []
    R2_to_merge = []
    for i in range(len(args.fastqs)):
        R1_to_merge.append(args.fastqs[i][0])
        if paired_end:
            R2_to_merge.append(args.fastqs[i][1])

    log.info('R1 to be merged: {}'.format(R1_to_merge))
    ret_val1 = pool.apply_async(merge_fastqs,
                    (R1_to_merge,args.out_dir,))
    R1_merged = ret_val1.get(BIG_INT)
    if paired_end:
        log.info('R2 to be merged: {}'.format(R2_to_merge))
        ret_val2 = pool.apply_async(merge_fastqs,
                        (R2_to_merge,args.out_dir,))
        R2_merged = ret_val2.get(BIG_INT)

    # close multithreading
    pool.close()
    pool.join()

    # write txt for output filenames (for WDL)
    log.info('Writinig output txt...')    
    txt = os.path.join(args.out_dir, args.out_txt)
    if paired_end:
        arr = [R1_merged, R2_merged]
    else:
        arr = [R1_merged]
    write_txt(txt, arr)

    log.info('All done.')

if __name__=='__main__':
    main()