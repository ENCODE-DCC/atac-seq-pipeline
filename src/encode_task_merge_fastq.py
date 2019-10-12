#!/usr/bin/env python

# ENCODE DCC fastq merger wrapper
# Author: Jin Lee (leepc12@gmail.com)

import sys
import os
import argparse
from encode_lib_common import (
    copy_f_to_f, log, ls_l, mkdir_p, read_tsv, run_shell_cmd,
    strip_ext_fastq)


def parse_arguments(debug=False):
    parser = argparse.ArgumentParser(prog='ENCODE DCC fastq merger.',
                                     description='')
    parser.add_argument(
        'fastqs', nargs='+', type=str,
        help='TSV file path or list of FASTQs. '
             'FASTQs must be compressed with gzip (with .gz). '
             'Use TSV for multiple fastqs to be merged later. '
             'row=merge_id, col=end_id).')
    parser.add_argument('--paired-end', action="store_true",
                        help='Paired-end FASTQs.')
    parser.add_argument('--nth', type=int, default=1,
                        help='Number of threads to parallelize.')
    parser.add_argument('--out-dir', default='', type=str,
                        help='Output directory.')
    parser.add_argument('--log-level', default='INFO',
                        choices=['NOTSET', 'DEBUG', 'INFO',
                                 'WARNING', 'CRITICAL', 'ERROR',
                                 'CRITICAL'],
                        help='Log level')
    args = parser.parse_args()

    # parse fastqs command line
    if args.fastqs[0].endswith('.gz') or args.fastqs[0].endswith('.fastq') or \
            args.fastqs[0].endswith('.fq'):  # it's fastq
        args.fastqs = [[f] for f in args.fastqs]  # make it a matrix
    else:  # it's TSV
        args.fastqs = read_tsv(args.fastqs[0])

    for i, fastqs in enumerate(args.fastqs):
        if args.paired_end and len(fastqs) != 2:
            raise argparse.ArgumentTypeError(
                'Need 2 fastqs per replicate for paired end.')
        if not args.paired_end and len(fastqs) != 1:
            raise argparse.ArgumentTypeError(
                'Need 1 fastq per replicate for single end.')

    log.setLevel(args.log_level)
    log.info(sys.argv)
    return args


def merge_fastqs(fastqs, end, out_dir):
    """make merged fastqs on $out_dir/R1, $out_dir/R2
    """
    out_dir = os.path.join(out_dir, end)
    mkdir_p(out_dir)
    prefix = os.path.join(out_dir,
                          os.path.basename(strip_ext_fastq(fastqs[0])))
    merged = '{}.merged.fastq.gz'.format(prefix)

    if len(fastqs) > 1:
        cmd = 'zcat -f {} | gzip -nc > {}'.format(
            ' '.join(fastqs),
            merged)
        run_shell_cmd(cmd)
        return merged
    else:
        return copy_f_to_f(fastqs[0], merged)


def main():
    # read params
    args = parse_arguments()

    log.info('Initializing and making output directory...')
    mkdir_p(args.out_dir)

    # update array with trimmed fastqs
    fastqs_R1 = []
    fastqs_R2 = []
    for fastqs in args.fastqs:
        fastqs_R1.append(fastqs[0])
        if args.paired_end:
            fastqs_R2.append(fastqs[1])

    log.info('Merging fastqs...')
    log.info('R1 to be merged: {}'.format(fastqs_R1))
    merge_fastqs(fastqs_R1, 'R1', args.out_dir)
    if args.paired_end:
        log.info('R2 to be merged: {}'.format(fastqs_R2))
        merge_fastqs(fastqs_R2, 'R2', args.out_dir)

    log.info('List all files in output directory...')
    ls_l(args.out_dir)

    log.info('All done.')


if __name__ == '__main__':
    main()
