#!/usr/bin/env python

# ENCODE DCC cross-correlation analysis wrapper
# Author: Jin Lee (leepc12@gmail.com)

import sys
import os
import argparse
from encode_common_genomic import *
from encode_common_log_parser import parse_xcor_score

def parse_arguments():
    parser = argparse.ArgumentParser(prog='ENCODE DCC cross-correlation analysis.',
                                        description='')
    parser.add_argument('ta', type=str,
                        help='Path for TAGALIGN file.')
    parser.add_argument('--subsample', type=int, default=0,
                        help='Subsample TAGALIGN.')
    parser.add_argument('--speak', type=int, default=-1,
                        help='User-defined cross-corr. peak strandshift \
                        (-speak= in run_spp.R). Disabled if -1.')    
    parser.add_argument('--paired-end', action="store_true",
                        help='Paired-end TAGALIGN.')
    parser.add_argument('--nth', type=int, default=1,
                        help='Number of threads to parallelize.')
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

def xcor(ta, speak, nth, out_dir):
    prefix = os.path.join(out_dir,
        os.path.basename(strip_ext_ta(ta)))
    xcor_plot_pdf = '{}.cc.plot.pdf'.format(prefix)
    xcor_score = '{}.cc.qc'.format(prefix)
    fraglen_txt = '{}.cc.fraglen.txt'.format(prefix)

    cmd1 = 'Rscript $(which run_spp.R) -rf -c={} -p={} '
    cmd1 += '-filtchr=chrM -savp={} -out={} {}'
    cmd1 = cmd1.format(
        ta,
        nth,
        xcor_plot_pdf,
        xcor_score,
        '-speak={}'.format(speak) if speak>=0 else '')
    run_shell_cmd(cmd1)

    cmd2 = 'sed -r \'s/,[^\\t]+//g\' -i {}'
    cmd2 = cmd2.format(xcor_score)
    run_shell_cmd(cmd2)

    # parse xcor_score and write fraglen (3rd column) to file
    cmd3 = 'echo {} > {}'.format(
        parse_xcor_score(xcor_score)['est_frag_len'],
        fraglen_txt)
    run_shell_cmd(cmd3)

    xcor_plot_png = pdf2png(xcor_plot_pdf, out_dir)
    return xcor_plot_pdf, xcor_plot_png, xcor_score, fraglen_txt

def main():
    # read params
    args = parse_arguments()
    log.info('Initializing and making output directory...')
    mkdir_p(args.out_dir)

    # declare temp arrays
    temp_files = [] # files to deleted later at the end

    if args.subsample:
        log.info('Subsampling TAGALIGN for xcor...')
        if args.paired_end:
            ta_subsampled = subsample_ta_pe(
                args.ta, args.subsample, True, True, args.out_dir)
        else:
            ta_subsampled = subsample_ta_se(
                args.ta, args.subsample, True, args.out_dir)
        temp_files.append(ta_subsampled)
    else:
        ta_subsampled = args.ta

    log.info('Cross-correlation analysis...')
    xcor_plot_pdf, xcor_plot_png, xcor_score, fraglen_txt = xcor(
        ta_subsampled, args.speak, args.nth, args.out_dir)

    log.info('Removing temporary files...')
    rm_f(temp_files)

    log.info('List all files in output directory...')
    ls_l(args.out_dir)

    log.info('All done.')

if __name__=='__main__':
    main()