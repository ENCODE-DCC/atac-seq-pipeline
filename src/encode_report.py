#!/usr/bin/env python

# ENCODE DCC reporting module wrapper
# Author: Jin Lee (leepc12@gmail.com)

import sys
import os
import argparse
from encode_common import *

def parse_arguments():
    parser = argparse.ArgumentParser(prog='ENCODE DCC generate HTML report and QC JSON.',
                                        description='')
    parser.add_argument('--fastqs-json', type=str,
                        help='JSON string for 3-dimensional array of FASTQ files.')
    parser.add_argument('--bams', type=str, nargs='*',
                        help='List of path for BAM files.')
    parser.add_argument('--nodup-bams', type=str, nargs='*',
                        help='List of path for filtered/dedup BAM files.')
    parser.add_argument('--tas', type=str, nargs='*',
                        help='List of path for TAG-ALIGN files.')
    parser.add_argument('--peaks', type=str, nargs='*',
                        help='List of path for PEAK files.')
    parser.add_argument('--peaks-pr1', type=str, nargs='*',
                        help='List of path for 1st pseudo replicated PEAK files.')
    parser.add_argument('--peaks-pr2', type=str, nargs='*',
                        help='List of path for 2nd pseudo replicated PEAK files.')
    parser.add_argument('--peak-ppr1', type=str,
                        help='List of path for files.')
    parser.add_argument('--peak-ppr2', type=str,
                        help='List of path for files.')
    parser.add_argument('--peak-pooled', type=str,
                        help='List of path for files.')



    parser.add_argument('--', type=str, nargs='*',
                        help='List of path for files.')
    parser.add_argument('--paired-end', action="store_true",
                        help='Paired-end TAGALIGN.')

    parser.add_argument('--out-report-html', default='report.html', type=str,
                        help='Output report HTML file.')
    parser.add_argument('--out-qc-json', default='qc.json', type=str,
                        help='Output QC JSON file.')
    parser.add_argument('--out-files-json', default='files.json', type=str,
                        help='Output files path/URI JSON file.')
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
    tmp_peak = '{}.tmp.gz'.format(prefix)
    xcor_plot = '{}.cc.plot.pdf'.format(prefix)
    xcor_score = '{}.cc.qc'.format(prefix)

    cmd1 = 'Rscript $(which run_spp.R) -rf -c={} -p={} '
    cmd1 += '-filtchr=chrM -savp={} -out={} {}'
    cmd1 = cmd1.format(
        ta,
        nth,
        xcor_plot,
        xcor_score,
        '-speak={}'.format(speak) if speak>=0 else '')
    run_shell_cmd(cmd1)

    cmd2 = 'sed -r \'s/,[^\\t]+//g\' -i {}'
    cmd2 = cmd2.format(xcor_score)
    run_shell_cmd(cmd2)
    return xcor_plot, xcor_score

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
    xcor_plot, xcor_score = xcor(
        ta_subsampled, args.speak, args.nth, args.out_dir)

    # remove temporary/intermediate files
    log.info('Removing temporary files...')
    rm_f(temp_files)

    log.info('List all files in output directory...')
    ls_l(args.out_dir)

    log.info('All done.')

if __name__=='__main__':
    main()