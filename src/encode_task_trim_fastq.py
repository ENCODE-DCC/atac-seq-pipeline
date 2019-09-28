#!/usr/bin/env python

# ENCODE DCC fastq merger wrapper
# Author: Jin Lee (leepc12@gmail.com)

import sys
import os
import argparse
from encode_lib_common import (
    assert_file_not_empty, copy_f_to_f, log, ls_l, mkdir_p,
    run_shell_cmd, strip_ext_fastq)


def parse_arguments(debug=False):
    parser = argparse.ArgumentParser(
        prog='ENCODE DCC fastq merger.')
    parser.add_argument('fastq', type=str,
                        help='FASTQ to be trimmed.')
    parser.add_argument('--trim-bp', type=int, default=50,
                        help='Number of basepair after trimming.')
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


def trim_fastq(fastq, trim_bp, out_dir):
    prefix = os.path.join(out_dir,
                          os.path.basename(strip_ext_fastq(fastq)))
    trimmed = '{}.trim_{}bp.fastq.gz'.format(prefix, trim_bp)

    cmd = 'python $(which trimfastq.py) {} {} | gzip -nc > {}'.format(
        fastq, trim_bp, trimmed)
    run_shell_cmd(cmd)

    # if shorter than trim_bp
    cmd2 = 'zcat -f {} | (grep \'sequences shorter than desired length\' '
    cmd2 += '|| true) | wc -l'
    cmd2 = cmd2.format(
        trimmed)
    if int(run_shell_cmd(cmd2)) > 0:
        copy_f_to_f(fastq, trimmed)

    return trimmed


def main():
    # read params
    args = parse_arguments()

    log.info('Initializing and making output directory...')
    mkdir_p(args.out_dir)

    log.info('Trimming fastqs ({} bp)...'.format(args.trim_bp))
    trimmed = trim_fastq(args.fastq, args.trim_bp, args.out_dir)
    assert_file_not_empty(trimmed)

    log.info('List all files in output directory...')
    ls_l(args.out_dir)

    log.info('All done.')


if __name__ == '__main__':
    main()
