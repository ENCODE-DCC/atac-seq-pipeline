#!/usr/bin/env python

# ENCODE DCC FRiP wrapper
# Author: Jin Lee (leepc12@gmail.com)

import sys
import os
import argparse
import multiprocessing
from encode_dcc_common import *

def parse_arguments():
    parser = argparse.ArgumentParser(prog='ENCODE DCC FRiP.',
                                        description='')
    parser.add_argument('peak', type=str,
                        help='Peak file.')
    parser.add_argument('ta', type=str,
                        help='TAGALIGN file.')
    parser.add_argument('--chrsz', type=str,
                        help='2-col chromosome sizes file.')
    parser.add_argument('--fraglen', type=int,
                        help='Fragment length for TAGALIGN file. \
                        If given, do shifted FRiP (for ChIP-Seq).')
    parser.add_argument('--out-dir', default='.', type=str,
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
        os.path.basename(strip_ext_all_genomic(bed)))
    frip_qc = '{}.frip.qc'.format(prefix)

    cmd = 'bedtools intersect -a <(zcat -f {}) '
    cmd += '-b <(zcat -f {}) -wa -u | wc -l'
    cmd = cmd.format(
        ta,
        peak)
    val1 = run_shell_cmd(cmd)
    val2 = get_num_lines(ta)
    write_txt(frip_qc, float(val1)/float(val2))
    return frip_qc

def frip_shifted(ta, peak, chrsz, fraglen, out_dir):
    prefix = os.path.join(out_dir, 
        os.path.basename(strip_ext_all_genomic(bed)))
    frip_qc = '{}.frip.qc'.format(prefix)
    half_fraglen = (fraglen+1)/2

    cmd = 'bedtools slop -i {} -g {} '
    cmd += '-s -l {} -r {} | '
    cmd += 'awk \'{{if ($2>=0 && $3>=0 && $2<=$3) print $0}}\' | '
    cmd += 'bedtools intersect -a stdin -b <(zcat -f {}) '
    cmd += '-wa -u | wc -l'
    cmd = cmd.format(
        ta,
        chrsz,
        -half_fraglen,
        helf_fraglen,
        peak)
    val1 = run_shell_cmd(cmd)
    val2 = get_num_lines(ta)
    write_txt(frip_qc, float(val1)/float(val2))
    return frip_qc

def main():
    # read params
    args = parse_arguments()
    log.info('Initializing and making output directory...')

    # make out_dir (root of all outputs)
    mkdir_p(args.out_dir)

    if args.fraglen:
        frip_qc = frip_shifted(args.ta, args.peak, 
            args.chrsz, args.fraglen, args.out_dir)
    else:
        frip_qc = frip(args.ta, args.peak, 
            args.out_dir)

    log.info('All done.')

if __name__=='__main__':
    main()