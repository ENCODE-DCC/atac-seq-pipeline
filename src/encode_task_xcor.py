#!/usr/bin/env python

# ENCODE DCC cross-correlation analysis wrapper
# Author: Jin Lee (leepc12@gmail.com)

import sys
import os
import argparse
from encode_lib_common import (
    log, ls_l, mkdir_p, pdf2png, rm_f, run_shell_cmd,
    strip_ext_ta)
from encode_lib_genomic import (
    subsample_ta_pe, subsample_ta_se)
from encode_lib_log_parser import parse_xcor_score


def parse_arguments():
    parser = argparse.ArgumentParser(
        prog='ENCODE DCC cross-correlation analysis.')
    parser.add_argument('ta', type=str,
                        help='Path for TAGALIGN file.')
    parser.add_argument('--mito-chr-name', default='chrM',
                        help='Mito chromosome name.')
    parser.add_argument('--subsample', type=int, default=0,
                        help='Subsample TAGALIGN.')
    parser.add_argument('--speak', type=int, default=-1,
                        help='User-defined cross-corr. peak strandshift \
                        (-speak= in run_spp.R). Disabled if -1.')
    parser.add_argument('--exclusion-range-min', type=int,
                        help='User-defined exclusion range minimum used for '
                             '-x=${xcor_exclusion_range_min}:'
                             '${xcor_exclusion_range_max}')
    parser.add_argument('--exclusion-range-max', type=int,
                        help='User-defined exclusion range maximum used for '
                             '-x=${xcor_exclusion_range_min}:'
                             '${xcor_exclusion_range_max}')
    parser.add_argument('--chip-seq-type', choices=['tf', 'histone'],
                        help='Type of ChIP-seq pipeline (histone of tf)')
    parser.add_argument('--paired-end', action="store_true",
                        help='Paired-end TAGALIGN.')
    parser.add_argument('--nth', type=int, default=1,
                        help='Number of threads to parallelize.')
    parser.add_argument('--out-dir', default='', type=str,
                        help='Output directory.')
    parser.add_argument('--log-level', default='INFO',
                        choices=['NOTSET', 'DEBUG', 'INFO',
                                 'WARNING', 'CRITICAL', 'ERROR',
                                 'CRITICAL'],
                        help='Log level')
    args = parser.parse_args()

    log.setLevel(args.log_level)
    log.info(sys.argv)
    return args


def get_exclusion_range_max(ta, chip_seq_type):
    # estimate read length from TA
    cmd0 = "zcat -f {} > {}.tmp".format(ta, ta)
    run_shell_cmd(cmd0)

    cmd = "head -n 100 {}.tmp | awk 'function abs(v) "
    cmd += "{{return v < 0 ? -v : v}} BEGIN{{sum=0}} "
    cmd += "{{sum+=abs($3-$2)}} END{{print int(sum/NR)}}'"
    cmd = cmd.format(ta)
    read_len = int(run_shell_cmd(cmd))
    rm_f(ta+'.tmp')

    if chip_seq_type == 'tf':
        return max(read_len + 10, 50)
    elif chip_seq_type == 'histone':
        return max(read_len + 10, 100)
    else:
        raise NotImplementedError('chip_seq_type not supported')


def xcor(ta, speak, mito_chr_name,
         nth, out_dir, chip_seq_type=None,
         exclusion_range_min=None, exclusion_range_max=None):
    prefix = os.path.join(out_dir,
                          os.path.basename(strip_ext_ta(ta)))
    xcor_plot_pdf = '{}.cc.plot.pdf'.format(prefix)
    xcor_score = '{}.cc.qc'.format(prefix)
    fraglen_txt = '{}.cc.fraglen.txt'.format(prefix)

    if chip_seq_type is not None and exclusion_range_min is not None:
        if exclusion_range_max is None:
            exclusion_range_max = get_exclusion_range_max(ta, chip_seq_type)

        exclusion_range_param = ' -x={}:{}'.format(
            exclusion_range_min, exclusion_range_max)
    else:
        exclusion_range_param = ''

    cmd1 = 'Rscript --max-ppsize=500000 $(which run_spp.R) -rf -c={} -p={} '
    cmd1 += '-filtchr="{}" -savp={} -out={} {}'
    cmd1 += exclusion_range_param
    cmd1 = cmd1.format(
        ta,
        nth,
        mito_chr_name,
        xcor_plot_pdf,
        xcor_score,
        '-speak={}'.format(speak) if speak >= 0 else '')
    run_shell_cmd(cmd1)

    cmd2 = 'sed -r \'s/,[^\\t]+//g\' -i {}'
    cmd2 = cmd2.format(xcor_score)
    run_shell_cmd(cmd2)

    # parse xcor_score and write fraglen (3rd column) to file
    cmd3 = 'echo {} > {}'.format(
        parse_xcor_score(xcor_score)['estimated_fragment_len'],
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
    temp_files = []  # files to deleted later at the end

    log.info('Subsampling TAGALIGN for xcor...')
    if args.paired_end:
        ta_subsampled = subsample_ta_pe(
            args.ta, args.subsample, True,
            args.mito_chr_name, True, args.out_dir)
    else:
        ta_subsampled = subsample_ta_se(
            args.ta, args.subsample, True,
            args.mito_chr_name, args.out_dir)
    temp_files.append(ta_subsampled)

    log.info('Cross-correlation analysis...')
    xcor_plot_pdf, xcor_plot_png, xcor_score, fraglen_txt = xcor(
        ta_subsampled, args.speak, args.mito_chr_name, args.nth, args.out_dir,
        args.chip_seq_type, args.exclusion_range_min, args.exclusion_range_max)

    log.info('Removing temporary files...')
    rm_f(temp_files)

    log.info('List all files in output directory...')
    ls_l(args.out_dir)

    log.info('All done.')


if __name__ == '__main__':
    main()
