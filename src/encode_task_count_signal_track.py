#!/usr/bin/env python

# ENCODE DCC Count signal track generation
# Author: Jin Lee (leepc12@gmail.com)

import sys
import os
import argparse
from encode_lib_common import (
    log, ls_l, mkdir_p, rm_f, run_shell_cmd, strip_ext_ta,
    get_gnu_sort_param,
)


def parse_arguments():
    parser = argparse.ArgumentParser(
        prog='ENCODE DCC Count signal track generation')
    parser.add_argument('ta', type=str,
                        help='Path for TAGALIGN file.')
    parser.add_argument('--chrsz', type=str,
                        help='2-col chromosome sizes file.')
    parser.add_argument('--mem-gb', type=float, default=4.0,
                        help='Max. memory for this job in GB. '
                        'This will be used to determine GNU sort -S (defaulting to 0.5 of this value). '
                        'It should be total memory for this task (not memory per thread).')
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


def count_signal_track(ta, chrsz, mem_gb, out_dir):
    prefix = os.path.join(out_dir, os.path.basename(strip_ext_ta(ta)))
    pos_bw = '{}.positive.bigwig'.format(prefix)
    neg_bw = '{}.negative.bigwig'.format(prefix)
    # temporary files
    pos_bedgraph = '{}.positive.bedgraph'.format(prefix)
    neg_bedgraph = '{}.negative.bedgraph'.format(prefix)

    temp_files = []

    run_shell_cmd(
        'zcat -f {ta} | sort -k1,1 -k2,2n {sort_param} | '
        'bedtools genomecov -5 -bg -strand + -g {chrsz} -i stdin > {pos_bedgraph}'.format(
            ta=ta,
            sort_param=get_gnu_sort_param(mem_gb * 1024 ** 3, ratio=0.5),
            chrsz=chrsz,
            pos_bedgraph=pos_bedgraph,
        )
    )

    run_shell_cmd(
        'zcat -f {ta} | sort -k1,1 -k2,2n {sort_param} | '
        'bedtools genomecov -5 -bg -strand - -g {chrsz} -i stdin > {neg_bedgraph}'.format(
            ta=ta,
            sort_param=get_gnu_sort_param(mem_gb * 1024 ** 3, ratio=0.5),
            chrsz=chrsz,
            neg_bedgraph=neg_bedgraph,
        )
    )

    run_shell_cmd(
        'bedGraphToBigWig {pos_bedgraph} {chrsz} {pos_bw}'.format(
            pos_bedgraph=pos_bedgraph,
            chrsz=chrsz,
            pos_bw=pos_bw,
        )
    )

    run_shell_cmd(
        'bedGraphToBigWig {neg_bedgraph} {chrsz} {neg_bw}'.format(
            neg_bedgraph=neg_bedgraph,
            chrsz=chrsz,
            neg_bw=neg_bw,
        )
    )

    # remove temporary files
    temp_files.append(pos_bedgraph)
    temp_files.append(neg_bedgraph)
    rm_f(temp_files)

    return pos_bw, neg_bw


def main():
    # read params
    args = parse_arguments()

    log.info('Initializing and making output directory...')
    mkdir_p(args.out_dir)

    log.info('Generating count signal tracks...')
    pos_bw, neg_bw = count_signal_track(
        args.ta,
        args.chrsz,
        args.mem_gb,
        args.out_dir
    )

    log.info('List all files in output directory...')
    ls_l(args.out_dir)

    log.info('All done.')


if __name__ == '__main__':
    main()
