#!/usr/bin/env python

# ENCODE DCC FRiP wrapper
# Author: Jin Lee (leepc12@gmail.com)

import sys
import os
import argparse
from encode_common import *

def parse_arguments():
    parser = argparse.ArgumentParser(prog='ENCODE DCC FRiP.',
                                        description='')
    parser.add_argument('peak', type=str,
                        help='Peak file.')
    parser.add_argument('ta', type=str,
                        help='TAGALIGN file.')
    parser.add_argument('--chrsz', type=str,
                        help='2-col chromosome sizes file. \
                        If given, do shifted FRiP (for ChIP-Seq).')
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

    log.setLevel(args.log_level)
    log.info(sys.argv)
    return args

def frip(ta, peak, out_dir):
    prefix = os.path.join(out_dir, 
        os.path.basename(strip_ext(peak)))
    frip_qc = '{}.frip.qc'.format(prefix)

    # due to bedtools bug when .gz is given for -a and -b
    tmp1 = gunzip(ta, 'tmp1', out_dir)
    tmp2 = gunzip(peak, 'tmp2', out_dir)    

    cmd = 'bedtools intersect -a {} -b {} -wa -u | wc -l'
    cmd = cmd.format(
        tmp1, # ta
        tmp2) # peak
    val1 = run_shell_cmd(cmd)
    val2 = get_num_lines(ta)
    write_txt(frip_qc, str(float(val1)/float(val2)))
    rm_f([tmp1, tmp2])
    return frip_qc

def frip_shifted(ta, peak, chrsz, fraglen, out_dir):
    prefix = os.path.join(out_dir, 
        os.path.basename(strip_ext(peak)))
    frip_qc = '{}.frip.qc'.format(prefix)
    half_fraglen = (fraglen+1)/2

    # due to bedtools bug when .gz is given for -a and -b
    tmp2 = gunzip(peak, 'tmp2', out_dir)    

    cmd = 'bedtools slop -i {} -g {} '
    cmd += '-s -l {} -r {} | '
    cmd += 'awk \'{{if ($2>=0 && $3>=0 && $2<=$3) print $0}}\' | '
    cmd += 'bedtools intersect -a stdin -b {} '
    cmd += '-wa -u | wc -l'
    cmd = cmd.format(
        ta,
        chrsz,
        -half_fraglen,
        half_fraglen,
        tmp2) # peak
    val1 = run_shell_cmd(cmd)
    val2 = get_num_lines(ta)
    write_txt(frip_qc, str(float(val1)/float(val2)))
    rm_f(tmp2)
    return frip_qc

def main():
    # read params
    args = parse_arguments()
    log.info('Initializing and making output directory...')
    mkdir_p(args.out_dir)

    if args.fraglen:
        frip_qc = frip_shifted(args.ta, args.peak, 
            args.chrsz, args.fraglen, args.out_dir)
    else:
        frip_qc = frip(args.ta, args.peak, args.out_dir)

    log.info('List all files in output directory...')
    ls_l(args.out_dir)

    log.info('All done.')

if __name__=='__main__':
    main()