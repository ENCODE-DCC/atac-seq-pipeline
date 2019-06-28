#!/usr/bin/env python

# ENCODE DCC blacklist filter wrapper
# Author: Jin Lee (leepc12@gmail.com)

import sys
import os
import argparse
from encode_lib_common import *

def parse_arguments():
    parser = argparse.ArgumentParser(prog='ENCODE DCC Blacklist filter.',
                                        description='')
    parser.add_argument('peak', type=str,
                        help='Peak file.')
    parser.add_argument('--blacklist', type=str, required=True,
                        help='Blacklist BED file.')
    parser.add_argument('--keep-irregular-chr', action="store_true",
                        help='Keep reads with non-canonical chromosome names.')    
    parser.add_argument('--out-dir', default='', type=str,
                        help='Output directory.')    
    parser.add_argument('--log-level', default='INFO', 
                        choices=['NOTSET','DEBUG','INFO',
                            'WARNING','CRITICAL','ERROR','CRITICAL'],
                        help='Log level')
    args = parser.parse_args()
    if args.blacklist.endswith('null'):
        args.blacklist = ''

    log.setLevel(args.log_level)
    log.info(sys.argv)
    return args

def blacklist_filter(peak, blacklist, keep_irregular_chr, out_dir):
    prefix = os.path.join(out_dir, 
        os.path.basename(strip_ext(peak)))
    peak_ext = get_ext(peak)
    filtered = '{}.bfilt.{}.gz'.format(prefix, peak_ext)

    if get_num_lines(peak)==0 or blacklist=='' or get_num_lines(blacklist)==0:
        cmd = 'zcat -f {} | gzip -nc > {}'.format(peak, filtered)
        run_shell_cmd(cmd)
    else:
        # due to bedtools bug when .gz is given for -a and -b
        tmp1 = gunzip(peak, 'tmp1', out_dir)
        tmp2 = gunzip(blacklist, 'tmp2', out_dir)

        cmd = 'bedtools intersect -nonamecheck -v -a {} -b {} | '
        cmd += 'awk \'BEGIN{{OFS="\\t"}} '
        cmd += '{{if ($5>1000) $5=1000; print $0}}\' | '
        if not keep_irregular_chr:
            cmd += 'grep -P \'chr[\\dXY]+\\b\' | '
        cmd += 'gzip -nc > {}'
        cmd = cmd.format(
            tmp1, # peak
            tmp2, # blacklist
            filtered)
        run_shell_cmd(cmd)
        rm_f([tmp1, tmp2])
    return filtered

def main():
    # read params
    args = parse_arguments()
    log.info('Initializing and making output directory...')

    # make out_dir (root of all outputs)
    mkdir_p(args.out_dir)

    # reproducibility QC
    log.info('Filtering peak with blacklist...')
    filtered = blacklist_filter(
                args.peak, args.blacklist, 
                args.keep_irregular_chr, args.out_dir)
    
    log.info('All done.')

if __name__=='__main__':
    main()