#!/usr/bin/env python

# ENCODE DCC Naive overlap wrapper
# Author: Jin Lee (leepc12@gmail.com)

import sys
import os
import argparse
from encode_common import *
from encode_blacklist_filter import blacklist_filter
from encode_frip import frip, frip_shifted

def parse_arguments():
    parser = argparse.ArgumentParser(prog='ENCODE DCC Naive overlap.',
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
                        choices=['narrowPeak','regionPeak','broadPeak','gappedPeak'],
                        help='Peak file type.')
    parser.add_argument('--nonamecheck', action='store_true',
                        help='bedtools intersect -nonamecheck. \
                        use this if you get bedtools intersect \
                        naming convenction warnings/errors).')
    parser.add_argument('--blacklist', type=str, required=True,
                        help='Blacklist BED file.')
    parser.add_argument('--ta', type=str,
                        help='TAGALIGN file for FRiP.')
    parser.add_argument('--chrsz', type=str,
                        help='2-col chromosome sizes file.')
    parser.add_argument('--fraglen', type=int, default=0,
                        help='Fragment length for TAGALIGN file. \
                        If given, do shifted FRiP (for ChIP-Seq).')
    parser.add_argument('--out-dir', default='', type=str,
                        help='Output directory.')
    parser.add_argument('--log-level', default='INFO', 
                        choices=['NOTSET','DEBUG','INFO',
                            'WARNING','CRITICAL','ERROR','CRITICAL'],
                        help='Log level')
    args = parser.parse_args()
    if args.blacklist.endswith('/dev/null'):
        args.blacklist = ''

    log.setLevel(args.log_level)
    log.info(sys.argv)
    return args

# only for narrowPeak (or regionPeak) type 
def naive_overlap(basename_prefix, peak1, peak2, peak_pooled, peak_type, 
    nonamecheck, out_dir):
    prefix = os.path.join(out_dir, basename_prefix)
    prefix += '.overlap'
    overlap_peak = '{}.{}.gz'.format(prefix, peak_type)

    nonamecheck_param = '-nonamecheck' if nonamecheck else  ''
    # narrowpeak, regionpeak only
    awk_param = '{s1=$3-$2; s2=$13-$12; '
    awk_param += 'if (($21/s1 >= 0.5) || ($21/s2 >= 0.5)) {print $0}}'
    cut_param = '1-10'

    # due to bedtools bug when .gz is given for -a and -b
    tmp1 = gunzip(peak1, 'tmp1', out_dir)
    tmp2 = gunzip(peak2, 'tmp2', out_dir)
    tmp_pooled = gunzip(peak_pooled, 'tmp_pooled', out_dir)

    # Find pooled peaks that overlap peak1 and peak2 
    # where overlap is defined as the fractional overlap 
    # wrt any one of the overlapping peak pairs >= 0.5
    cmd1 = 'intersectBed {} -wo '
    cmd1 += '-a {} -b {} | '
    cmd1 += 'awk \'BEGIN{{FS="\\t";OFS="\\t"}} {}\' | '
    cmd1 += 'cut -f {} | sort | uniq | '
    cmd1 += 'intersectBed {} -wo '
    cmd1 += '-a stdin -b {} | '
    cmd1 += 'awk \'BEGIN{{FS="\\t";OFS="\\t"}} {}\' | '
    cmd1 += 'cut -f {} | sort | uniq | gzip -nc > {}'
    cmd1 = cmd1.format(
        nonamecheck_param,
        tmp_pooled, # peak_pooled
        tmp1, # peak1
        awk_param,
        cut_param,
        nonamecheck_param,
        tmp2, # peak2
        awk_param,
        cut_param,
        overlap_peak)
    run_shell_cmd(cmd1)
    rm_f([tmp1,tmp2,tmp_pooled])
    return overlap_peak

def main():
    # read params
    args = parse_arguments()

    log.info('Initializing and making output directory...')
    mkdir_p(args.out_dir)

    log.info('Do naive overlap...')
    overlap_peak = naive_overlap(
        args.prefix, args.peak1, args.peak2, args.peak_pooled, 
        args.peak_type, args.nonamecheck, args.out_dir)

    log.info('Checking if output is empty...') # bedtools issue    
    assert_file_not_empty(overlap_peak)

    if args.blacklist:    
        log.info('Blacklist-filtering peaks...')
        bfilt_overlap_peak = blacklist_filter(
                overlap_peak, args.blacklist, False, args.out_dir)
    else:
        bfilt_overlap_peak = overlap_peak        

    if args.ta: # if TAG-ALIGN is given
        if args.fraglen: # chip-seq
            log.info('Shifted FRiP with fragment length...')
            frip_qc = frip_shifted( args.ta, bfilt_overlap_peak,
                args.chrsz, args.fraglen, args.out_dir)
        else: # atac-seq
            log.info('FRiP without fragment length...')
            frip_qc = frip( args.ta, bfilt_overlap_peak, args.out_dir)
    else:
        frip_qc = '/dev/null'

    log.info('List all files in output directory...')
    ls_l(args.out_dir)

    log.info('All done.')

if __name__=='__main__':
    main()
