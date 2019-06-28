#!/usr/bin/env python

# ENCODE DCC Count signal track generation
# Author: Jin Lee (leepc12@gmail.com)

import sys
import os
import argparse
from encode_common import *

def parse_arguments():
    parser = argparse.ArgumentParser(prog='ENCODE DCC Count signal track generation',
                                        description='')
    parser.add_argument('ta', type=str,
                        help='Path for TAGALIGN file.')
    parser.add_argument('--chrsz', type=str,
                        help='2-col chromosome sizes file.')
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

def count_signal_track(ta, chrsz, out_dir):
    prefix = os.path.join(out_dir, os.path.basename(strip_ext_ta(ta)))
    pos_bw = '{}.positive.bigwig'.format(prefix)
    neg_bw = '{}.negative.bigwig'.format(prefix)
    # temporary files
    pos_bedgraph = '{}.positive.bedgraph'.format(prefix)
    neg_bedgraph = '{}.negative.bedgraph'.format(prefix)

    temp_files = []

    cmd1 = 'zcat -f {} | sort -k1,1 -k2,2n | '
    cmd1 += 'bedtools genomecov -5 -bg -strand + -g {} -i stdin > {}'
    cmd1 = cmd1.format(ta, chrsz, pos_bedgraph)
    run_shell_cmd(cmd1)

    cmd2 = 'zcat -f {} | sort -k1,1 -k2,2n | '
    cmd2 += 'bedtools genomecov -5 -bg -strand - -g {} -i stdin > {}'
    cmd2 = cmd2.format(ta, chrsz, neg_bedgraph)
    run_shell_cmd(cmd2)

    cmd3 = 'bedGraphToBigWig {} {} {}'
    cmd3 = cmd3.format(pos_bedgraph, chrsz, pos_bw)
    run_shell_cmd(cmd3)    

    cmd4 = 'bedGraphToBigWig {} {} {}'
    cmd4 = cmd4.format(neg_bedgraph, chrsz, neg_bw)
    run_shell_cmd(cmd4)

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
    pos_bw, neg_bw = count_signal_track( args.ta, args.chrsz, args.out_dir)

    log.info('List all files in output directory...')
    ls_l(args.out_dir)

    log.info('All done.')

if __name__=='__main__':
    main()

