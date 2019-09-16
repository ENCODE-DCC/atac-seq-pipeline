#!/usr/bin/env python

# ENCODE DCC spp call peak wrapper
# Author: Jin Lee (leepc12@gmail.com)

import sys
import os
import argparse
from encode_lib_common import *

def parse_arguments():
    parser = argparse.ArgumentParser(prog='ENCODE spp call_peak',
                                        description='')
    parser.add_argument('tas', type=str, nargs=2,
                        help='Path for TAGALIGN file and control TAGALIGN file.')
    parser.add_argument('--chrsz', type=str,
                        help='2-col chromosome sizes file.')
    parser.add_argument('--fraglen', type=int, required=True,
                        help='Fragment length.')
    parser.add_argument('--cap-num-peak', default=300000, type=int,
                        help='Capping number of peaks by taking top N peaks.')
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

def spp(ta, ctl_ta, fraglen, cap_num_peak, nth, out_dir):
    basename_ta = os.path.basename(strip_ext_ta(ta))
    basename_ctl_ta = os.path.basename(strip_ext_ta(ctl_ta))
    basename_prefix = '{}_x_{}'.format(basename_ta, basename_ctl_ta)
    if len(basename_prefix) > 200: # UNIX cannot have filename > 255
        basename_prefix = '{}_x_control'.format(basename_ta)
    nth_param = '-p={}'.format(nth) if nth<2 else ''
    prefix = os.path.join(out_dir, basename_prefix)
    rpeak = '{}.{}.regionPeak.gz'.format(
        prefix,
        human_readable_number(cap_num_peak))
    rpeak_tmp = '{}.tmp'.format(rpeak)
    rpeak_tmp_gz = '{}.tmp.gz'.format(rpeak)

    cmd0 = 'Rscript --max-ppsize=500000 $(which run_spp.R) -c={} -i={} '
    cmd0 += '-npeak={} -odir={} -speak={} -savr={} -rf {}'
    cmd0 = cmd0.format(
        ta,
        ctl_ta,
        cap_num_peak,
        os.path.abspath(out_dir),
        fraglen,
        rpeak_tmp,
        nth_param)
    run_shell_cmd(cmd0)

    # if we have scientific representation of chr coord. then convert it to int
    cmd1 = 'zcat -f {} | awk \'BEGIN{{OFS="\\t"}}'
    cmd1 += '{{if ($2<0) $2=0; print $1,int($2),int($3),$4,$5,$6,$7,$8,$9,$10;}}\' | '
    cmd1 += 'gzip -f -nc > {}'
    cmd1 = cmd1.format(
        rpeak_tmp,
        rpeak)
    run_shell_cmd(cmd1)
    rm_f([rpeak_tmp, rpeak_tmp_gz])

    return rpeak

def main():
    # read params
    args = parse_arguments()

    log.info('Initializing and making output directory...')
    mkdir_p(args.out_dir)

    log.info('Calling peaks with spp...')
    rpeak = spp(args.tas[0], args.tas[1], 
        args.fraglen, args.cap_num_peak, args.nth, args.out_dir)

    log.info('Checking if output is empty...')
    assert_file_not_empty(rpeak)

    log.info('List all files in output directory...')
    ls_l(args.out_dir)

    log.info('All done.')

if __name__=='__main__':
    main()
