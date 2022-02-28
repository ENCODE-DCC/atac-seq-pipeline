#!/usr/bin/env python

# ENCODE DCC spp call peak wrapper
# Author: Jin Lee (leepc12@gmail.com)

import sys
import os
import argparse
from encode_lib_common import (
    assert_file_not_empty, human_readable_number, log,
    ls_l, mkdir_p, rm_f, run_shell_cmd, strip_ext_ta)
from encode_lib_genomic import (
    subsample_ta_se, subsample_ta_pe, bed_clip)


def parse_arguments():
    parser = argparse.ArgumentParser(
        prog='ENCODE spp call_peak')
    parser.add_argument(
        'tas', type=str, nargs=2,
        help='Path for TAGALIGN file and control TAGALIGN file.')
    parser.add_argument('--chrsz', type=str,
                        help='2-col chromosome sizes file.')
    parser.add_argument('--fraglen', type=int, required=True,
                        help='Fragment length.')
    parser.add_argument('--fdr-thresh', default=0.01, type=float,
                        help='FDR threshold for run_spp.R -fdr parameter.')
    parser.add_argument('--cap-num-peak', default=300000, type=int,
                        help='Capping number of peaks by taking top N peaks.')
    parser.add_argument('--ctl-subsample', default=0, type=int,
                        help='Subsample control to this read depth '
                             '(0: no subsampling).')
    parser.add_argument('--ctl-paired-end', action="store_true",
                        help='Paired-end control TA.')
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


def spp(ta, ctl_ta, chrsz, fraglen, cap_num_peak, fdr_thresh,
        ctl_subsample, ctl_paired_end, nth, out_dir):
    basename_ta = os.path.basename(strip_ext_ta(ta))

    if ctl_subsample:
        if ctl_paired_end:
            ctl_ta = subsample_ta_pe(
                ctl_ta, ctl_subsample,
                non_mito=False, mito_chr_name=None, r1_only=False,
                out_dir=out_dir)
        else:
            ctl_ta = subsample_ta_se(
                ctl_ta, ctl_subsample,
                non_mito=False, mito_chr_name=None,
                out_dir=out_dir)
    basename_ctl_ta = os.path.basename(strip_ext_ta(ctl_ta))
    basename_prefix = '{}_x_{}'.format(basename_ta, basename_ctl_ta)
    if len(basename_prefix) > 200:  # UNIX cannot have filename > 255
        basename_prefix = '{}_x_control'.format(basename_ta)
    nth_param = '-p={}'.format(nth) if nth >= 2 else ''
    prefix = os.path.join(out_dir, basename_prefix)
    rpeak = '{}.{}.regionPeak.gz'.format(
        prefix,
        human_readable_number(cap_num_peak))
    rpeak_tmp_prefix = '{}.tmp'.format(rpeak)
    rpeak_tmp_gz = '{}.tmp.gz'.format(rpeak)
    rpeak_tmp2 = '{}.tmp2'.format(rpeak)

    cmd0 = 'Rscript --max-ppsize=500000 $(which run_spp.R) -c={} -i={} '
    cmd0 += '-npeak={} -odir={} -speak={} -savr={} -fdr={} -rf {}'
    cmd0 = cmd0.format(
        ta,
        ctl_ta,
        cap_num_peak,
        os.path.abspath(out_dir),
        fraglen,
        rpeak_tmp_prefix,
        fdr_thresh,
        nth_param)
    run_shell_cmd(cmd0)

    # if we have scientific representation of chr coord. then convert it to int
    cmd1 = 'zcat -f {} | awk \'BEGIN{{OFS="\\t"}}'
    cmd1 += '{{if ($2<0) $2=0; '
    cmd1 += 'print $1,int($2),int($3),$4,$5,$6,$7,$8,$9,$10;}}\' > {}'
    cmd1 = cmd1.format(
        rpeak_tmp_gz,
        rpeak_tmp2)
    run_shell_cmd(cmd1)
    rm_f(rpeak_tmp_gz)

    # clip peaks between 0-chromSize.
    bed_clip(rpeak_tmp2, chrsz, rpeak)
    rm_f(rpeak_tmp2)

    return rpeak


def main():
    # read params
    args = parse_arguments()

    log.info('Initializing and making output directory...')
    mkdir_p(args.out_dir)

    log.info('Calling peaks with spp...')
    rpeak = spp(args.tas[0], args.tas[1], args.chrsz,
                args.fraglen, args.cap_num_peak, args.fdr_thresh,
                args.ctl_subsample, args.ctl_paired_end, args.nth, args.out_dir)

    log.info('Checking if output is empty...')
    assert_file_not_empty(rpeak, help=
        'No peaks found. FDR threshold (fdr_thresh in your input JSON) '
        'might be too stringent or poor quality sample?')

    log.info('List all files in output directory...')
    ls_l(args.out_dir)

    log.info('All done.')


if __name__ == '__main__':
    main()
