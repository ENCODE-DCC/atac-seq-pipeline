#!/usr/bin/env python

# ENCODE DCC pseudo replicator wrapper
# Author: Jin Lee (leepc12@gmail.com)

import sys
import os
import argparse
import multiprocessing
from encode_common import *

def parse_arguments():
    parser = argparse.ArgumentParser(prog='ENCODE DCC pseudo replicator.',
                                        description='')
    parser.add_argument('ta', type=str,
                        help='Path for TAGALIGN file.')
    parser.add_argument('--paired-end', action="store_true",
                        help='Paired-end TAGALIGN.')
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

def spr_se(ta, out_dir):
    prefix = os.path.join(out_dir,
        os.path.basename(strip_ext_ta(ta)))
    tmp_pr1 = '{}.00'.format(prefix)
    tmp_pr2 = '{}.01'.format(prefix)
    ta_pr1 = '{}.pr1.tagAlign.gz'.format(prefix)
    ta_pr2 = '{}.pr2.tagAlign.gz'.format(prefix)
    nlines = int((get_num_lines(ta)+1)/2)
    
    # bash-only
    cmd1 = 'zcat {} | shuf --random-source=<(openssl enc -aes-256-ctr -pass pass:$(zcat -f {} | wc -c) -nosalt </dev/zero 2>/dev/null) | '
    cmd1 += 'split -d -l {} - {}.'
    cmd1 = cmd1.format(
        ta,
        ta,
        nlines,
        prefix)
    run_shell_cmd(cmd1)

    cmd2 = 'gzip -nc {} > {}'
    cmd2 = cmd2.format(
        tmp_pr1,
        ta_pr1)
    run_shell_cmd(cmd2)

    cmd3 = 'gzip -nc {} > {}'
    cmd3 = cmd3.format(
        tmp_pr2,
        ta_pr2)
    run_shell_cmd(cmd3)

    rm_f([tmp_pr1, tmp_pr2])
    return ta_pr1, ta_pr2

def spr_pe(ta, out_dir):
    prefix = os.path.join(out_dir,
        os.path.basename(strip_ext_ta(ta)))
    tmp_pr1 = '{}.00'.format(prefix)
    tmp_pr2 = '{}.01'.format(prefix)
    ta_pr1 = '{}.pr1.tagAlign.gz'.format(prefix)
    ta_pr2 = '{}.pr2.tagAlign.gz'.format(prefix)
    nlines = int((get_num_lines(ta)/2+1)/2)

    # bash-only
    cmd1 = 'zcat -f {} | sed \'N;s/\\n/\\t/\' | '
    cmd1 += 'shuf --random-source=<(openssl enc -aes-256-ctr -pass pass:$(zcat -f {} | wc -c) -nosalt </dev/zero 2>/dev/null) | '
    cmd1 += 'split -d -l {} - {}.'
    cmd1 = cmd1.format(
        ta,
        ta,
        nlines,
        prefix)
    run_shell_cmd(cmd1)

    cmd2 = 'zcat -f {} | '
    cmd2 += 'awk \'BEGIN{{OFS="\\t"}} '
    cmd2 += '{{printf "%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n'
    cmd2 += '%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n",'
    cmd2 += '$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}}\' | '
    cmd2 += 'gzip -nc > {}'
    cmd2 = cmd2.format(
        tmp_pr1,
        ta_pr1)
    run_shell_cmd(cmd2)

    cmd3 = 'zcat -f {} | '
    cmd3 += 'awk \'BEGIN{{OFS="\\t"}} '
    cmd3 += '{{printf "%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n'
    cmd3 += '%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n",'
    cmd3 += '$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}}\' | '
    cmd3 += 'gzip -nc > {}'
    cmd3 = cmd3.format(
        tmp_pr2,
        ta_pr2)
    run_shell_cmd(cmd3)

    rm_f([tmp_pr1, tmp_pr2])
    return ta_pr1, ta_pr2

def main():
    # read params
    args = parse_arguments()
    log.info('Initializing and making output directory...')
    mkdir_p(args.out_dir)

    log.info('Making self-pseudo replicates...')
    if args.paired_end:
        ta_pr1, ta_pr2 = spr_pe(args.ta, args.out_dir)
    else:
        ta_pr1, ta_pr2 = spr_se(args.ta, args.out_dir)

    log.info('List all files in output directory...')
    ls_l(args.out_dir)

    log.info('Checking if output is empty...')
    assert_file_not_empty(ta_pr1)
    assert_file_not_empty(ta_pr2)

    log.info('All done.')

if __name__=='__main__':
    main()