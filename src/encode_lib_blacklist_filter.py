#!/usr/bin/env python

# ENCODE DCC blacklist filter wrapper
# Author: Jin Lee (leepc12@gmail.com)

import sys
import os
import argparse
from encode_lib_common import (
    get_ext, get_num_lines, gunzip, log, mkdir_p,
    rm_f, run_shell_cmd, strip_ext, strip_ext_bam)


def parse_arguments():
    parser = argparse.ArgumentParser(prog='ENCODE DCC Blacklist filter.')
    parser.add_argument('peak', type=str, help='Peak file.')
    parser.add_argument('--blacklist', type=str,
                        help='Blacklist BED file.')
    parser.add_argument('--regex-bfilt-peak-chr-name',
                        help='Keep chromosomes matching this pattern only '
                             'in .bfilt. peak files.')
    parser.add_argument('--out-dir', default='', type=str,
                        help='Output directory.')
    parser.add_argument('--log-level', default='INFO',
                        choices=['NOTSET', 'DEBUG', 'INFO',
                                 'WARNING', 'CRITICAL', 'ERROR',
                                 'CRITICAL'],
                        help='Log level')
    args = parser.parse_args()
    if args.blacklist is None or args.blacklist.endswith('null'):
        args.blacklist = ''

    log.setLevel(args.log_level)
    log.info(sys.argv)
    return args


def blacklist_filter(peak, blacklist, regex_bfilt_peak_chr_name, out_dir):
    prefix = os.path.join(
        out_dir,
        os.path.basename(strip_ext(peak)))
    peak_ext = get_ext(peak)
    filtered = '{}.bfilt.{}.gz'.format(prefix, peak_ext)
    if regex_bfilt_peak_chr_name is None:
        regex_bfilt_peak_chr_name = ''

    if blacklist is None or blacklist == '' or get_num_lines(peak) == 0 \
            or get_num_lines(blacklist) == 0:
        cmd = 'zcat -f {} | '
        cmd += 'grep -P \'{}\\b\' | '
        cmd += 'gzip -nc > {}'
        cmd = cmd.format(
            peak,
            regex_bfilt_peak_chr_name,
            filtered)
        run_shell_cmd(cmd)
    else:
        # due to bedtools bug when .gz is given for -a and -b
        tmp1 = gunzip(peak, 'tmp1', out_dir)
        tmp2 = gunzip(blacklist, 'tmp2', out_dir)

        cmd = 'bedtools intersect -nonamecheck -v -a {} -b {} | '
        cmd += 'awk \'BEGIN{{OFS="\\t"}} '
        cmd += '{{if ($5>1000) $5=1000; print $0}}\' | '
        cmd += 'grep -P \'{}\\b\' | '
        cmd += 'gzip -nc > {}'
        cmd = cmd.format(
            tmp1,  # peak
            tmp2,  # blacklist
            regex_bfilt_peak_chr_name, # regex
            filtered)
        run_shell_cmd(cmd)
        rm_f([tmp1, tmp2])
    return filtered


def blacklist_filter_bam(bam, blacklist, out_dir):
    prefix = os.path.join(out_dir,
                          os.path.basename(strip_ext_bam(bam)))
    filtered = '{}.bfilt.bam'.format(prefix)

    if blacklist == '' or get_num_lines(blacklist) == 0:
        cmd = 'zcat -f {} | gzip -nc > {}'.format(bam, filtered)
        run_shell_cmd(cmd)
    else:
        # due to bedtools bug when .gz is given for -a and -b
        tmp2 = gunzip(blacklist, 'tmp2', out_dir)

        cmd = 'bedtools intersect -nonamecheck -v -abam {} -b {} > {}'
        cmd = cmd.format(
            bam,
            tmp2,  # blacklist
            filtered)
        run_shell_cmd(cmd)
        rm_f([tmp2])
    return filtered


def main():
    # read params
    args = parse_arguments()
    log.info('Initializing and making output directory...')

    # make out_dir (root of all outputs)
    mkdir_p(args.out_dir)

    # reproducibility QC
    log.info('Filtering peak with blacklist...')
    blacklist_filter(
        args.peak, args.blacklist,
        args.keep_irregular_chr, args.out_dir)

    log.info('All done.')


if __name__ == '__main__':
    main()
