#!/usr/bin/env python

# ENCODE DCC MACS2 signal track wrapper
# Author: Jin Lee (leepc12@gmail.com)

import sys
import os
import argparse
from encode_lib_common import *

def parse_arguments():
    parser = argparse.ArgumentParser(prog='ENCODE DCC MACS2 signal track',
                                        description='')
    parser.add_argument('ta', type=str,
                        help='Path for TAGALIGN file.')
    parser.add_argument('--chrsz', type=str,
                        help='2-col chromosome sizes file.')
    parser.add_argument('--gensz', type=str,
                        help='Genome size (sum of entries in 2nd column of \
                            chr. sizes file, or hs for human, ms for mouse).')
    parser.add_argument('--pval-thresh', default=0.01, type=float,
                        help='P-Value threshold.')
    parser.add_argument('--smooth-win', default=73, type=int,
                        help='Smoothing window size.')
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

def macs2_signal_track(ta, chrsz, gensz, pval_thresh, smooth_win, out_dir):
    prefix = os.path.join(out_dir,
        os.path.basename(strip_ext_ta(ta)))
    fc_bigwig = '{}.fc.signal.bigwig'.format(prefix)
    pval_bigwig = '{}.pval.signal.bigwig'.format(prefix)
    # temporary files
    fc_bedgraph = '{}.fc.signal.bedgraph'.format(prefix)
    fc_bedgraph_srt = '{}.fc.signal.srt.bedgraph'.format(prefix)
    pval_bedgraph = '{}.pval.signal.bedgraph'.format(prefix)
    pval_bedgraph_srt = '{}.pval.signal.srt.bedgraph'.format(prefix)

    shiftsize = -int(round(float(smooth_win)/2.0))
    temp_files = []

    cmd0 = 'macs2 callpeak '
    cmd0 += '-t {} -f BED -n {} -g {} -p {} '
    cmd0 += '--shift {} --extsize {} '
    cmd0 += '--nomodel -B --SPMR '
    cmd0 += '--keep-dup all --call-summits '
    cmd0 = cmd0.format(
        ta,
        prefix,
        gensz,
        pval_thresh,
        shiftsize,
        smooth_win)
    run_shell_cmd(cmd0)

    cmd3 = 'macs2 bdgcmp -t "{}"_treat_pileup.bdg '
    cmd3 += '-c "{}"_control_lambda.bdg '
    cmd3 += '--o-prefix "{}" -m FE '
    cmd3 = cmd3.format(
        prefix, 
        prefix, 
        prefix)
    run_shell_cmd(cmd3)

    cmd4 = 'bedtools slop -i "{}"_FE.bdg -g {} -b 0 | '
    cmd4 += 'bedClip stdin {} {}'
    cmd4 = cmd4.format(
        prefix, 
        chrsz, 
        chrsz, 
        fc_bedgraph)
    run_shell_cmd(cmd4)
  
    # sort and remove any overlapping regions in bedgraph by comparing two lines in a row
    cmd5 = 'LC_COLLATE=C sort -k1,1 -k2,2n {} | ' \
        'awk \'BEGIN{{OFS="\\t"}}{{if (NR==1 || NR>1 && (prev_chr!=$1 '\
        '|| prev_chr==$1 && prev_chr_e<=$2)) ' \
        '{{print $0}}; prev_chr=$1; prev_chr_e=$3;}}\' > {}'.format(
        fc_bedgraph,
        fc_bedgraph_srt)
    run_shell_cmd(cmd5)
    rm_f(fc_bedgraph)

    cmd6 = 'bedGraphToBigWig {} {} {}'
    cmd6 = cmd6.format(
        fc_bedgraph_srt,
        chrsz,
        fc_bigwig)
    run_shell_cmd(cmd6)
    rm_f(fc_bedgraph_srt)

    # sval counts the number of tags per million in the (compressed) BED file
    sval = float(get_num_lines(ta))/1000000.0
    
    cmd7 = 'macs2 bdgcmp -t "{}"_treat_pileup.bdg '
    cmd7 += '-c "{}"_control_lambda.bdg '
    cmd7 += '--o-prefix {} -m ppois -S {}'
    cmd7 = cmd7.format(
        prefix,
        prefix,
        prefix,
        sval)
    run_shell_cmd(cmd7)

    cmd8 = 'bedtools slop -i "{}"_ppois.bdg -g {} -b 0 | '
    cmd8 += 'bedClip stdin {} {}'
    cmd8 = cmd8.format(
        prefix,
        chrsz,
        chrsz,
        pval_bedgraph)
    run_shell_cmd(cmd8)

    # sort and remove any overlapping regions in bedgraph by comparing two lines in a row
    cmd9 = 'LC_COLLATE=C sort -k1,1 -k2,2n {} | ' \
        'awk \'BEGIN{{OFS="\\t"}}{{if (NR==1 || NR>1 && (prev_chr!=$1 '\
        '|| prev_chr==$1 && prev_chr_e<=$2)) ' \
        '{{print $0}}; prev_chr=$1; prev_chr_e=$3;}}\' > {}'.format(
        pval_bedgraph,
        pval_bedgraph_srt)
    run_shell_cmd(cmd9)
    rm_f(pval_bedgraph)

    cmd10 = 'bedGraphToBigWig {} {} {}'
    cmd10 = cmd10.format(
        pval_bedgraph_srt,
        chrsz,
        pval_bigwig)
    run_shell_cmd(cmd10)
    rm_f(pval_bedgraph_srt)
    
    # remove temporary files
    temp_files.append("{}_*".format(prefix))
    rm_f(temp_files)

    return fc_bigwig, pval_bigwig

def main():
    # read params
    args = parse_arguments()

    log.info('Initializing and making output directory...')
    mkdir_p(args.out_dir)

    log.info('Calling peaks and generating signal tracks with MACS2...')
    fc_bigwig, pval_bigwig = macs2_signal_track(
        args.ta, args.chrsz, args.gensz, args.pval_thresh,
        args.smooth_win, args.out_dir)

    log.info('List all files in output directory...')
    ls_l(args.out_dir)

    log.info('All done.')

if __name__=='__main__':
    main()
