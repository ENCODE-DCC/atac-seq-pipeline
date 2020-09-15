#!/usr/bin/env python
import sys
import os
import argparse
from encode_lib_common import (
    assert_file_not_empty, get_num_lines, log, ls_l, mkdir_p, rm_f,
    run_shell_cmd, strip_ext_ta)
from encode_lib_genomic import (
    subsample_ta_pe, subsample_ta_se)

def parse_arguments():
    parser = argparse.ArgumentParser(
        prog='ENCODE DCC control TAG-ALIGN subsampler.'
             'This script does not check if number of reads in TA is higher than '
             'subsampling number (--subsample). '
             'If number of reads in TA is lower than subsampling number then '
             'TA will be just shuffled.')
    parser.add_argument('ta', type=str,
                        help='Path for control TAGALIGN file.')
    parser.add_argument('--paired-end', action="store_true",
                        help='Paired-end TAGALIGN.')
    parser.add_argument('--subsample', default=0, type=int,
                        help='Number of reads to subsample.')
    parser.add_argument('--out-dir', default='', type=str,
                        help='Output directory.')
    parser.add_argument('--log-level', default='INFO',
                        choices=['NOTSET', 'DEBUG', 'INFO',
                                 'WARNING', 'CRITICAL', 'ERROR',
                                 'CRITICAL'],
                        help='Log level')
    args = parser.parse_args()
    if not args.subsample:
        raise ValueError('--subsample should be a positive integer.')

    log.setLevel(args.log_level)
    log.info(sys.argv)
    return args


def main():
    # read params
    args = parse_arguments()
    log.info('Initializing and making output directory...')
    mkdir_p(args.out_dir)

    if args.paired_end:
        subsampled_ta = subsample_ta_pe(
            args.ta, args.subsample,
            non_mito=False, mito_chr_name=None, r1_only=False,
            out_dir=args.out_dir)
    else:
        subsampled_ta = subsample_ta_se(
            args.ta, args.subsample,
            non_mito=False, mito_chr_name=None,
            out_dir=args.out_dir)
    log.info('Checking if output is empty...')
    assert_file_not_empty(subsampled_ta)

    log.info('List all files in output directory...')
    ls_l(args.out_dir)

    log.info('All done.')


if __name__ == '__main__':
    main()
