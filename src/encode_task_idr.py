#!/usr/bin/env python

# ENCODE DCC IDR wrapper
# Author: Jin Lee (leepc12@gmail.com)

import sys
import os
import argparse
import math
from encode_lib_common import (
    assert_file_not_empty,
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
    bed_clip,
    peak_to_starch,
)
from encode_lib_blacklist_filter import blacklist_filter
from encode_lib_frip import frip, frip_shifted


def parse_arguments():
    parser = argparse.ArgumentParser(
        prog='ENCODE DCC IDR.',
        description='NarrowPeak or RegionPeak only.')
    parser.add_argument('peak1', type=str,
                        help='Peak file 1.')
    parser.add_argument('peak2', type=str,
                        help='Peak file 2.')
    parser.add_argument('peak_pooled', type=str,
                        help='Pooled peak file.')
    parser.add_argument('--prefix', default='idr', type=str,
                        help='Prefix basename for output IDR peak.')
    parser.add_argument('--peak-type', type=str, required=True,
                        choices=['narrowPeak', 'regionPeak',
                                 'broadPeak', 'gappedPeak'],
                        help='Peak file type.')
    parser.add_argument('--idr-thresh', default=0.1, type=float,
                        help='IDR threshold.')
    parser.add_argument('--idr-rank', default='p.value', type=str,
                        choices=['p.value', 'q.value', 'signal.value'],
                        help='IDR ranking method.')
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


def get_npeak_col_by_rank(rank):
    if rank == 'signal.value':
        return 7
    elif rank == 'p.value':
        return 8
    elif rank == 'q.value':
        return 9
    else:
        raise Exception('Invalid score ranking method')

# only for narrowPeak (or regionPeak) type


def idr(basename_prefix, peak1, peak2, peak_pooled, peak_type, chrsz,
        thresh, rank, mem_gb, out_dir):
    prefix = os.path.join(out_dir, basename_prefix)
    prefix += '.idr{}'.format(thresh)
    idr_peak = '{}.{}.gz'.format(prefix, peak_type)
    idr_plot = '{}.unthresholded-peaks.txt.png'.format(prefix)
    idr_stdout = '{}.log'.format(prefix)
    # temporary
    idr_12col_bed = '{}.12-col.bed.gz'.format(peak_type)
    idr_out = '{}.unthresholded-peaks.txt'.format(prefix)
    idr_tmp = '{}.unthresholded-peaks.txt.tmp'.format(prefix)
    idr_out_gz = '{}.unthresholded-peaks.txt.gz'.format(prefix)

    run_shell_cmd(
        'idr --samples {peak1} {peak2} --peak-list {peak_pooled} --input-file-type narrowPeak '
        '--output-file {idr_out} --rank {rank} --soft-idr-threshold {thresh} '
        '--plot --use-best-multisummit-IDR --log-output-file {idr_stdout}'.format(
            peak1=peak1,
            peak2=peak2,
            peak_pooled=peak_pooled,
            idr_out=idr_out,
            rank=rank,
            thresh=thresh,
            idr_stdout=idr_stdout,
        )
    )

    # clip peaks between 0-chromSize.
    bed_clip(idr_out, chrsz, idr_tmp, no_gz=True)

    col = get_npeak_col_by_rank(rank)
    neg_log10_thresh = -math.log10(thresh)
    # LC_COLLATE=C
    run_shell_cmd(
        'awk \'BEGIN{{OFS="\\t"}} $12>={neg_log10_thresh} '
        '{{if ($2<0) $2=0; '
        'print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}}\' {idr_tmp} '
        '| sort {sort_param} | uniq | sort -grk{col},{col} {sort_param} | gzip -nc > {idr_12col_bed}'.format(
            neg_log10_thresh=neg_log10_thresh,
            idr_tmp=idr_tmp,
            sort_param=get_gnu_sort_param(mem_gb * 1024 ** 3, ratio=0.5),
            col=col,
            idr_12col_bed=idr_12col_bed,
        )
    )

    run_shell_cmd(
        'zcat {idr_12col_bed} | '
        'awk \'BEGIN{{OFS="\\t"}} '
        '{{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10}}\' | '
        'gzip -nc > {idr_peak}'.format(
            idr_12col_bed=idr_12col_bed,
            idr_peak=idr_peak,
        )
    )

    run_shell_cmd(
        'cat {idr_tmp} | gzip -nc > {idr_out_gz}'.format(
            idr_tmp=idr_tmp,
            idr_out_gz=idr_out_gz,
        )
    )

    rm_f([idr_out, idr_tmp, idr_12col_bed])
    rm_f('{prefix}.*.noalternatesummitpeaks.png'.format(prefix=prefix))
    return idr_peak, idr_plot, idr_out_gz, idr_stdout


def main():
    # read params
    args = parse_arguments()

    log.info('Initializing and making output directory...')
    mkdir_p(args.out_dir)

    log.info('Do IDR...')
    idr_peak, idr_plot, idr_out_gz, idr_stdout = idr(
        args.prefix,
        args.peak1, args.peak2, args.peak_pooled, args.peak_type,
        args.chrsz,
        args.idr_thresh, args.idr_rank, args.mem_gb, args.out_dir,
    )

    log.info('Checking if output is empty...')
    assert_file_not_empty(idr_peak, help=
        'No IDR peaks found. IDR threshold might be too stringent '
        'or replicates have very poor concordance.')

    log.info('Blacklist-filtering peaks...')
    bfilt_idr_peak = blacklist_filter(
        idr_peak, args.blacklist, args.regex_bfilt_peak_chr_name, args.out_dir)

    log.info('Converting peak to bigbed...')
    peak_to_bigbed(bfilt_idr_peak, args.peak_type, args.chrsz,
                   args.mem_gb, args.out_dir)

    log.info('Converting peak to starch...')
    peak_to_starch(bfilt_idr_peak, args.out_dir)

    log.info('Converting peak to hammock...')
    peak_to_hammock(bfilt_idr_peak, args.mem_gb, args.out_dir)

    if args.ta:  # if TAG-ALIGN is given
        if args.fraglen:  # chip-seq
            log.info('Shifted FRiP with fragment length...')
            frip_shifted(args.ta, bfilt_idr_peak,
                         args.chrsz, args.fraglen, args.out_dir)
        else:  # atac-seq
            log.info('FRiP without fragment length...')
            frip(args.ta, bfilt_idr_peak, args.out_dir)

    log.info('List all files in output directory...')
    ls_l(args.out_dir)

    log.info('All done.')


if __name__ == '__main__':
    main()
