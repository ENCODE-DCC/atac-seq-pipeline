#!/usr/bin/env python

# Author: Jin Lee (leepc12@gmail.com)

import sys
import argparse
from encode_lib_common import (
    assert_file_not_empty,
    log,
    ls_l,
    mkdir_p,
)
from encode_lib_genomic import (
    peak_to_bigbed,
    peak_to_hammock,
    get_region_size_metrics,
    get_num_peaks,
    peak_to_starch,
)
from encode_lib_blacklist_filter import blacklist_filter
from encode_lib_frip import frip


def parse_arguments():
    parser = argparse.ArgumentParser(prog='ENCODE post_call_peak (atac)',
                                     description='')
    parser.add_argument(
        'peak', type=str,
        help='Path for PEAK file. Peak filename should be "*.*Peak.gz". '
             'e.g. rep1.narrowPeak.gz')
    parser.add_argument('--ta', type=str,
                        help='TAG-ALIGN file.')
    parser.add_argument('--peak-type', type=str, required=True,
                        choices=['narrowPeak', 'regionPeak',
                                 'broadPeak', 'gappedPeak'],
                        help='Peak file type.')
    parser.add_argument('--chrsz', type=str,
                        help='2-col chromosome sizes file.')
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


def main():
    # read params
    args = parse_arguments()

    log.info('Initializing and making output directory...')
    mkdir_p(args.out_dir)

    log.info('Blacklist-filtering peaks...')
    bfilt_peak = blacklist_filter(
        args.peak, args.blacklist, args.regex_bfilt_peak_chr_name, args.out_dir)

    log.info('Checking if output is empty...')
    assert_file_not_empty(bfilt_peak)

    log.info('Converting peak to bigbed...')
    peak_to_bigbed(bfilt_peak, args.peak_type, args.chrsz,
                   args.out_dir)

    log.info('Converting peak to starch...')
    peak_to_starch(bfilt_peak, args.out_dir)

    log.info('Converting peak to hammock...')
    peak_to_hammock(bfilt_peak, args.out_dir)

    log.info('FRiP without fragment length...')
    frip(args.ta, bfilt_peak, args.out_dir)

    log.info('Calculating (blacklist-filtered) peak region size QC/plot...')
    get_region_size_metrics(bfilt_peak)

    log.info('Calculating number of peaks (blacklist-filtered)...')
    get_num_peaks(bfilt_peak)

    log.info('List all files in output directory...')
    ls_l(args.out_dir)

    log.info('All done.')


if __name__ == '__main__':
    main()
