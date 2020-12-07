#!/usr/bin/env python

# Author: Jin Lee (leepc12@gmail.com)

import sys
import os
import argparse
from encode_lib_common import (
    log,
    ls_l,
    mkdir_p,
    rm_f,
)
from encode_lib_genomic import (
    bam_to_pbam,
)


def parse_arguments():
    parser = argparse.ArgumentParser(prog='ENCODE bam to pbam',
                                     description='')
    parser.add_argument('bam', type=str,
                        help='Path for BAM.')
    parser.add_argument('--ref-fa', type=str,
                        help='Path for reference fasta.')
    parser.add_argument('--delete-original-bam', action='store_true',
                        help='Delete original BAM after conversion.')
    parser.add_argument('--out-dir', default='', type=str,
                        help='Output directory.')
    parser.add_argument('--log-level', default='INFO',
                        choices=['NOTSET', 'DEBUG', 'INFO',
                                 'WARNING', 'CRITICAL', 'ERROR',
                                 'CRITICAL'],
                        help='Log level')
    args = parser.parse_args()

    log.setLevel(args.log_level)
    log.info(sys.argv)
    return args


def main():
    # read params
    args = parse_arguments()

    log.info('Initializing and making output directory...')
    mkdir_p(args.out_dir)

    # generate read length file
    log.info('Converting BAM into pBAM...')
    bam_to_pbam(args.bam, args.ref_fa, args.out_dir)

    if args.delete_original_bam:
        log.info('Deleting original BAM...')
        rm_f(args.bam)

    log.info('List all files in output directory...')
    ls_l(args.out_dir)

    log.info('All done.')


if __name__ == '__main__':
    main()
