#!/usr/bin/env python

# ENCODE DCC Naive overlap wrapper
# Author: Jin Lee (leepc12@gmail.com)

import sys
import os
import argparse
from encode_lib_common import (
    assert_file_not_empty,
    gunzip,
    log,
    ls_l,
    mkdir_p,
    rm_f,
    run_shell_cmd,
    get_gnu_sort_param,
)
from encode_lib_genomic import (
    peak_to_bigbed,
    peak_to_hammock,
    peak_to_starch,
)
from encode_lib_blacklist_filter import blacklist_filter
from encode_lib_frip import frip, frip_shifted


def parse_arguments():
    parser = argparse.ArgumentParser(
        prog='ENCODE DCC Naive overlap.',
        description='NarrowPeak or RegionPeak only.')
    parser.add_argument('peak1', type=str,
                        help='Peak 1.')
    parser.add_argument('peak2', type=str,
                        help='Peak 2.')
    parser.add_argument('peak_pooled', type=str,
                        help='Pooled peak.')
    parser.add_argument('--prefix', default='overlap', type=str,
                        help='Prefix basename for output overlap peak.')
    parser.add_argument('--peak-type', type=str, required=True,
                        choices=['narrowPeak', 'regionPeak',
                                 'broadPeak', 'gappedPeak'],
                        help='Peak file type.')
    parser.add_argument('--nonamecheck', action='store_true',
                        help='bedtools intersect -nonamecheck. \
                        use this if you get bedtools intersect \
                        naming convenction warnings/errors).')
    parser.add_argument('--blacklist', type=str,
                        help='Blacklist BED file.')
    parser.add_argument('--regex-bfilt-peak-chr-name',
                        help='Keep chromosomes matching this pattern only '
                             'in .bfilt. peak files.')
    parser.add_argument('--ta', type=str,
                        help='TAGALIGN file for FRiP.')
    parser.add_argument('--chrsz', type=str,
                        help='2-col chromosome sizes file.')
    parser.add_argument('--fraglen', type=int, default=0,
                        help='Fragment length for TAGALIGN file. \
                        If given, do shifted FRiP (for ChIP-Seq).')
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
    if args.blacklist is None or args.blacklist.endswith('null'):
        args.blacklist = ''

    log.setLevel(args.log_level)
    log.info(sys.argv)
    return args

# only for narrowPeak (or regionPeak) type


def naive_overlap(basename_prefix, peak1, peak2, peak_pooled, peak_type,
                  nonamecheck, mem_gb, out_dir):
    prefix = os.path.join(out_dir, basename_prefix)
    prefix += '.overlap'
    overlap_peak = '{}.{}.gz'.format(prefix, peak_type)

    nonamecheck_param = '-nonamecheck' if nonamecheck else ''
    if peak_type.lower() in ('narrowpeak', 'regionpeak'):
        awk_param = '{s1=$3-$2; s2=$13-$12; if (($21/s1 >= 0.5) || ($21/s2 >= 0.5)) {print $0}}'
        cut_param = '1-10'
    elif peak_type.lower() == 'broadpeak':
        awk_param = '{s1=$3-$2; s2=$12-$11; if (($19/s1 >= 0.5) || ($19/s2 >= 0.5)) {print $0}}'
        cut_param = '1-9'
    elif peak_type.lower() == 'gappedpeak':
        awk_param = '{s1=$3-$2; s2=$18-$17; if (($31/s1 >= 0.5) || ($31/s2 >= 0.5)) {print $0}}'
        cut_param = '1-15'
    else:
        raise ValueError('Unsupported peak_type.')

    # due to bedtools bug when .gz is given for -a and -b
    tmp1 = gunzip(peak1, 'tmp1', out_dir)
    tmp2 = gunzip(peak2, 'tmp2', out_dir)
    tmp_pooled = gunzip(peak_pooled, 'tmp_pooled', out_dir)

    # Find pooled peaks that overlap peak1 and peak2
    # where overlap is defined as the fractional overlap
    # wrt any one of the overlapping peak pairs >= 0.5
    run_shell_cmd(
        'intersectBed {nonamecheck_param} -wo '
        '-a {tmp_pooled} -b {tmp1} | '
        'awk \'BEGIN{{FS="\\t";OFS="\\t"}} {awk_param}\' | '
        'cut -f {cut_param} | sort {sort_param} | uniq | '
        'intersectBed {nonamecheck_param} -wo '
        '-a stdin -b {tmp2} | '
        'awk \'BEGIN{{FS="\\t";OFS="\\t"}} {awk_param}\' | '
        'cut -f {cut_param} | sort {sort_param} | uniq | gzip -nc > {overlap_peak}'.format(
            nonamecheck_param=nonamecheck_param,
            tmp_pooled=tmp_pooled, # peak_pooled
            tmp1=tmp1, # peak1
            awk_param=awk_param,
            cut_param=cut_param,
            sort_param=get_gnu_sort_param(mem_gb * 1024 ** 3, ratio=0.5),
            tmp2=tmp2, # peak2
            overlap_peak=overlap_peak,
        )
    )
    rm_f([tmp1, tmp2, tmp_pooled])
    return overlap_peak


def main():
    # read params
    args = parse_arguments()

    log.info('Initializing and making output directory...')
    mkdir_p(args.out_dir)

    log.info('Do naive overlap...')
    overlap_peak = naive_overlap(
        args.prefix, args.peak1, args.peak2, args.peak_pooled,
        args.peak_type, args.nonamecheck, args.mem_gb, args.out_dir,
    )

    log.info('Blacklist-filtering peaks...')
    bfilt_overlap_peak = blacklist_filter(
        overlap_peak, args.blacklist, args.regex_bfilt_peak_chr_name, args.out_dir)

    log.info('Checking if output is empty...')
    assert_file_not_empty(bfilt_overlap_peak)

    log.info('Converting peak to bigbed...')
    peak_to_bigbed(bfilt_overlap_peak, args.peak_type,
                   args.chrsz, args.mem_gb, args.out_dir)

    log.info('Converting peak to starch...')
    peak_to_starch(bfilt_overlap_peak, args.out_dir)

    log.info('Converting peak to hammock...')
    peak_to_hammock(bfilt_overlap_peak, args.mem_gb, args.out_dir)

    if args.ta:  # if TAG-ALIGN is given
        if args.fraglen:  # chip-seq
            log.info('Shifted FRiP with fragment length...')
            frip_shifted(args.ta, bfilt_overlap_peak,
                         args.chrsz, args.fraglen, args.out_dir)
        else:  # atac-seq
            log.info('FRiP without fragment length...')
            frip(args.ta, bfilt_overlap_peak, args.out_dir)

    log.info('List all files in output directory...')
    ls_l(args.out_dir)

    log.info('All done.')


if __name__ == '__main__':
    main()
