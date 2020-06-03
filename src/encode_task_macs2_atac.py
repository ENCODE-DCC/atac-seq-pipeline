#!/usr/bin/env python

# ENCODE DCC MACS2 call peak wrapper
# Author: Jin Lee (leepc12@gmail.com)

import sys
import os
import argparse
from encode_lib_common import (
    assert_file_not_empty, human_readable_number,
    log, ls_l, mkdir_p, rm_f, run_shell_cmd, strip_ext_ta)
from encode_lib_genomic import bed_clip


def parse_arguments():
    parser = argparse.ArgumentParser(prog='ENCODE MACS2 callpeak',
                                     description='')
    parser.add_argument('ta', type=str,
                        help='Path for TAGALIGN file.')
    parser.add_argument('--chrsz', type=str,
                        help='2-col chromosome sizes file.')
    parser.add_argument('--gensz', type=str,
                        help='Genome size (sum of entries in 2nd column of \
                            chr. sizes file, or hs for human, ms for mouse).')
    parser.add_argument('--pval-thresh', default=0.01, type=float,
                        help='P-Value threshold.')
    parser.add_argument('--smooth-win', default=150, type=int,
                        help='Smoothing window size.')
    parser.add_argument('--cap-num-peak', default=300000, type=int,
                        help='Capping number of peaks by taking top N peaks.')
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


def macs2(ta, chrsz, gensz, pval_thresh, smooth_win, cap_num_peak,
          out_dir):
    prefix = os.path.join(out_dir,
                          os.path.basename(strip_ext_ta(ta)))
    npeak = '{}.{}.{}.narrowPeak.gz'.format(
        prefix,
        'pval{}'.format(pval_thresh),
        human_readable_number(cap_num_peak))
    # temporary files
    npeak_tmp = '{}.tmp'.format(npeak)
    npeak_tmp2 = '{}.tmp2'.format(npeak)

    shiftsize = -int(round(float(smooth_win)/2.0))
    temp_files = []

    cmd0 = 'macs2 callpeak '
    cmd0 += '-t {} -f BED -n {} -g {} -p {} '
    cmd0 += '--shift {} --extsize {} '
    cmd0 += '--nomodel -B --SPMR '
    cmd0 += '--keep-dup all --call-summits '
    cmd0 = cmd0.format(
        ta,
        prefix,
        gensz,
        pval_thresh,
        shiftsize,
        smooth_win)
    run_shell_cmd(cmd0)

    cmd1 = 'LC_COLLATE=C sort -k 8gr,8gr "{}"_peaks.narrowPeak | '
    cmd1 += 'awk \'BEGIN{{OFS="\\t"}}'
    cmd1 += '{{$4="Peak_"NR; if ($2<0) $2=0; if ($3<0) $3=0; if ($10==-1) '
    cmd1 += '$10=$2+int(($3-$2+1)/2.0); print $0}}\' > {}'
    cmd1 = cmd1.format(
        prefix,
        npeak_tmp)
    run_shell_cmd(cmd1)

    cmd2 = 'head -n {} {} > {}'.format(
        cap_num_peak,
        npeak_tmp,
        npeak_tmp2)
    run_shell_cmd(cmd2)

    # clip peaks between 0-chromSize.
    bed_clip(npeak_tmp2, chrsz, npeak)

    rm_f([npeak_tmp, npeak_tmp2])

    # remove temporary files
    temp_files.append("{}_*".format(prefix))
    rm_f(temp_files)

    return npeak


def main():
    # read params
    args = parse_arguments()

    log.info('Initializing and making output directory...')
    mkdir_p(args.out_dir)

    log.info('Calling peaks macs2...')
    npeak = macs2(
        args.ta, args.chrsz, args.gensz, args.pval_thresh,
        args.smooth_win, args.cap_num_peak,
        args.out_dir)

    log.info('Checking if output is empty...')
    assert_file_not_empty(npeak)

    log.info('List all files in output directory...')
    ls_l(args.out_dir)

    log.info('All done.')


if __name__ == '__main__':
    main()
