#!/usr/bin/env python

# ENCODE DCC MACS2 call peak wrapper
# Author: Jin Lee (leepc12@gmail.com)

import sys
import os
import argparse
from encode_lib_common import *

def parse_arguments():
    parser = argparse.ArgumentParser(prog='ENCODE DCC MACS2 callpeak',
                                        description='')
    parser.add_argument('tas', type=str, nargs='+',
                        help='Path for TAGALIGN file (first) and control TAGALIGN file (second; optional).')
    parser.add_argument('--fraglen', type=int, required=True,
                        help='Fragment length.')
    parser.add_argument('--shift', type=int, default=0,
                        help='macs2 callpeak --shift.')
    parser.add_argument('--chrsz', type=str,
                        help='2-col chromosome sizes file.')
    parser.add_argument('--gensz', type=str,
                        help='Genome size (sum of entries in 2nd column of \
                            chr. sizes file, or hs for human, ms for mouse).')
    parser.add_argument('--pval-thresh', default=0.01, type=float,
                        help='P-Value threshold.')
    parser.add_argument('--cap-num-peak', default=500000, type=int,
                        help='Capping number of peaks by taking top N peaks.')
    parser.add_argument('--out-dir', default='', type=str,
                        help='Output directory.')
    parser.add_argument('--log-level', default='INFO', 
                        choices=['NOTSET','DEBUG','INFO',
                            'WARNING','CRITICAL','ERROR','CRITICAL'],
                        help='Log level')
    args = parser.parse_args()
    if len(args.tas)==1:
        args.tas.append('')
    log.setLevel(args.log_level)
    log.info(sys.argv)
    return args

def macs2(ta, ctl_ta, chrsz, gensz, pval_thresh, shift, fraglen, cap_num_peak,
        out_dir):
    basename_ta = os.path.basename(strip_ext_ta(ta))
    if ctl_ta:
        basename_ctl_ta = os.path.basename(strip_ext_ta(ctl_ta))
        basename_prefix = '{}_x_{}'.format(basename_ta, basename_ctl_ta)
        if len(basename_prefix) > 200: # UNIX cannot have len(filename) > 255
            basename_prefix = '{}_x_control'.format(basename_ta)
    else:
        basename_prefix = basename_ta
    prefix = os.path.join(out_dir, basename_prefix)
    npeak = '{}.{}.{}.narrowPeak.gz'.format(
        prefix,
        'pval{}'.format(pval_thresh),
        human_readable_number(cap_num_peak))
    npeak_tmp = '{}.tmp'.format(npeak)
    fc_bigwig = '{}.fc.signal.bigwig'.format(prefix)
    pval_bigwig = '{}.pval.signal.bigwig'.format(prefix)
    # temporary files
    fc_bedgraph = '{}.fc.signal.bedgraph'.format(prefix)
    fc_bedgraph_srt = '{}.fc.signal.srt.bedgraph'.format(prefix)
    pval_bedgraph = '{}.pval.signal.bedgraph'.format(prefix)
    pval_bedgraph_srt = '{}.pval.signal.srt.bedgraph'.format(prefix)

    temp_files = []

    cmd0 = ' macs2 callpeak '
    cmd0 += '-t {} {} -f BED -n {} -g {} -p {} '
    cmd0 += '--nomodel --shift {} --extsize {} --keep-dup all -B --SPMR'
    cmd0 = cmd0.format(
        ta,
        '-c {}'.format(ctl_ta) if ctl_ta else '',
        prefix,
        gensz,
        pval_thresh,
        0,
        fraglen)
    run_shell_cmd(cmd0)

    cmd1 = 'LC_COLLATE=C sort -k 8gr,8gr "{}"_peaks.narrowPeak | '
    cmd1 += 'awk \'BEGIN{{OFS="\\t"}}'
    cmd1 += '{{$4="Peak_"NR; if ($2<0) $2=0; if ($3<0) $3=0; if ($10==-1) '
    cmd1 += '$10=$2+int(($3-$2+1)/2.0); print $0}}\' > {}'
    cmd1 = cmd1.format(
        prefix,
        npeak_tmp)
    run_shell_cmd(cmd1)

    cmd2 = 'head -n {} {} | gzip -nc > {}'.format(
        cap_num_peak,
        npeak_tmp,
        npeak)
    run_shell_cmd(cmd2)
    rm_f(npeak_tmp)
    
    # remove temporary files
    temp_files.append("{}_*".format(prefix))
    rm_f(temp_files)

    return npeak

def main():
    # read params
    args = parse_arguments()

    log.info('Initializing and making output directory...')
    mkdir_p(args.out_dir)

    log.info('Calling peaks with macs2...')
    npeak = macs2(
        args.tas[0], args.tas[1], args.chrsz, args.gensz, args.pval_thresh,
        args.shift, args.fraglen, args.cap_num_peak, args.out_dir)

    log.info('Checking if output is empty...')
    assert_file_not_empty(npeak)

    log.info('List all files in output directory...')
    ls_l(args.out_dir)

    log.info('All done.')

if __name__=='__main__':
    main()