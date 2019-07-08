#!/usr/bin/env python

# Author: Jin Lee (leepc12@gmail.com)

import sys
import os
import re
import argparse
from encode_lib_genomic import *

def parse_arguments():
    parser = argparse.ArgumentParser(prog='ENCODE post align',
                                        description='')
    parser.add_argument('fastq', type=str,
                        help='Path for FASTQ R1')
    parser.add_argument('bam', type=str,
                        help='Path for BAM')
    parser.add_argument('--out-dir', default='', type=str,
                            help='Output directory.')
    parser.add_argument('--log-level', default='INFO', 
                        choices=['NOTSET','DEBUG','INFO',
                            'WARNING','CRITICAL','ERROR','CRITICAL'],
                        help='Log level')
    args = parser.parse_args()

    log.setLevel(args.log_level)
    log.info(sys.argv)    
    return args

def make_read_length_file(fastq, out_dir):
    basename = os.path.basename(strip_ext_fastq(fastq))
    prefix = os.path.join(out_dir, basename)
    txt = '{}.read_length.txt'.format(prefix)
    read_length = get_read_length(fastq)
    with open(txt,'w') as fp:
        fp.write(str(read_length))
    return txt

def main():
    # read params
    args = parse_arguments()

    log.info('Initializing and making output directory...')
    mkdir_p(args.out_dir)

    # generate read length file
    log.info('Generating read length file...')
    R1_read_length_file = make_read_length_file(
                            args.fastq, args.out_dir)

    log.info('Running samtools index...')
    bai = samtools_index(args.bam, args.out_dir)

    log.info('SAMstat...')
    flagstat = samstat(bam, args.nth, args.out_dir)

    log.info('List all files in output directory...')
    ls_l(args.out_dir)

    log.info('All done.')

if __name__=='__main__':
    main()
