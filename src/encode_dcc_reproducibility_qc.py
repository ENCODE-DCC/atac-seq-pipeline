#!/usr/bin/env python

# ENCODE DCC FRiP and reproducibility QC (IDR and naive-overlap) wrapper
# Author: Jin Lee (leepc12@gmail.com)

import sys
import os
import argparse
import multiprocessing
from encode_dcc_common import *

def parse_arguments():
    parser = argparse.ArgumentParser(prog='ENCODE DCC FRiP and reproducibility QC.',
                                        description='')
    parser.add_argument('--peaks', type=str, nargs='+', required=True,
                        help='List of peak files from true replicates.')
    parser.add_argument('--peaks-pr1', type=str, nargs='+', required=True,
                        help='List of peak files from 1st pseudo replicates.')
    parser.add_argument('--peaks-pr2', type=str, nargs='+', required=True,
                        help='List of peak files from 2nd pseudo replicates.')
    parser.add_argument('--peak-pooled', type=str, required=True,
                        help='Peak file from pooled replicate.')
    parser.add_argument('--peak-ppr1', type=str, required=True,
                        help='Peak file from 1st pooled pseudo replicate.')
    parser.add_argument('--peak-ppr2', type=str, required=True,
                        help='Peak file from 2nd pooled pseudo replicate.')
    parser.add_argument('--method', type=str, required=True,
                        choices=['idr','overlap']
                        help='Post-processing method for raw peaks.')
    parser.add_argument('--idr-thresh', default=0.1, type=float,
                        help='IDR threshold.')
    parser.add_argument('--idr-rank', default='p.value', type=str,
                        choices=['p.value','q.value','signal.value'],
                        help='IDR ranking method.')
    parser.add_argument('--blacklist', type=str,
                        help='Blacklist BED file.')
    parser.add_argument('--out-dir', default='.', type=str,
                        help='Output directory.')
    parser.add_argument('--log-level', default='INFO', 
                        choices=['NOTSET','DEBUG','INFO',
                            'WARNING','CRITICAL','ERROR','CRITICAL'],
                        help='Log level')
    args = parser.parse_args()

    if len(args.peaks)!=len(args.peaks_pr1) or \
        len(args.peaks)!=len(args.peaks_pr2):
        raise ValueError('--peaks, --peaks-pr1 and --peaks-pr2 '+
            'must have the same dimension.')

    log.setLevel(args.log_level)
    log.info(sys.argv)
    return args

def filter_out_blacklist(bed, blacklist):
    prefix = os.path.join(out_dir, 
        os.path.basename(strip_ext_all_genomic(bed)))
    return

def filter_out_chr(bed, pattern):
    return

def idr(basename_prefix, peak1, peak2, peak_pooled, 
    thresh, rank, out_dir):
    prefix = os.path.join(out_dir, basename_prefix)

    npeak = '{}.narrowPeak.gz'.format(prefix)
    fc_bigwig = '{}.fc.signal.bigwig'.format(prefix)
    pval_bigwig = '{}.pval.signal.bigwig'.format(prefix)
    # temporary files
    fc_bedgraph = '{}.fc.signal.bedgraph'.format(prefix)
    fc_bedgraph_srt = '{}.fc.signal.srt.bedgraph'.format(prefix)
    pval_bedgraph = '{}.pval.signal.bedgraph'.format(prefix)
    pval_bedgraph_srt = '{}.pval.signal.srt.bedgraph'.format(prefix)

    shiftsize = -round(float(smooth_win)/2.0)
    temp_files = []

    cmd1 = 'export LC_COLLATE=C && macs2 callpeak '
    cmd1 += '-t {} -f BED -g {} -p {} '
    cmd1 += '--shift {} --extsize {} '
    cmd1 += '-n  '
    cmd1 += '--nomodel -B --SPMR '
    cmd1 += '--keep-dup all --call-summits '
    cmd1 = cmd1.format(
        ta,
        gensz,
        pval_thresh,
        shiftsize,
        smooth_window,
        prefix)
    run_shell_cmd(cmd1)
    return

def overlap(basename_prefix, peak1, peak2, peak_pooled, out_dir):
    filter_out_blacklist()
    return

def frip(ta, peak, out_dir):
    prefix = os.path.join(out_dir, 
        os.path.basename(strip_ext_all_genomic(bed)))
    frip_qc = '{}.frip.qc'.format(prefix)

    cmd1 = 'bedtools intersect -a <(zcat -f {}) '
    cmd1 += '-b <(zcat -f {}) -wa -u | wc -l'
    cmd1 = cmd1.format(
        ta,
        peak)
    val1 = run_shell_cmd(cmd1)    
    val2 = get_num_lines(ta)
    write_txt(frip_qc, float(val1)/float(val2))
    return frip_qc

def frip_shifted(ta, peak, chrsz, fraglen, out_dir):
    prefix = os.path.join(out_dir, 
        os.path.basename(strip_ext_all_genomic(bed)))
    frip_qc = '{}.frip.qc'.format(prefix)
    half_fraglen = (fraglen+1)/2

    cmd1 = 'bedtools slop -i {} -g {} '
    cmd1 += '-s -l {} -r {} | '
    cmd1 += 'awk \'{{if ($2>=0 && $3>=0 && $2<=$3) print $0}}\' | '
    cmd1 += 'bedtools intersect -a stdin -b <(zcat -f {}) '
    cmd1 += '-wa -u | wc -l'
    cmd1 = cmd1.format(
        ta,
        chrsz,
        -half_fraglen,
        helf_fraglen,
        peak)
    val1 = run_shell_cmd(cmd1)
    val2 = get_num_lines(ta)
    write_txt(frip_qc, float(val1)/float(val2))
    return frip_qc

def main():
    # read params
    args = parse_arguments()
    log.info('Initializing and making output directory...')

    # make out_dir (root of all outputs)
    mkdir_p(args.out_dir)

    # initialize multithreading
    log.info('Initializing multithreading...')
    num_process = min(3,args.nth)
    log.info('Number of threads={}.'.format(num_process))
    pool = multiprocessing.Pool(num_process)

    log.info('Post-processing raw peaks...')
    ret_val_true={}
    ret_val_pr={}
    for i in range(len(args.peaks)):
        # for every pair of true replicates
        for j in range(i+1, len(args.peaks)):
            basename_prefix = 
                '{}.rep{}_vs_rep{}'.format(t,i+1,j+1)
            if t=='idr':
                ret_val_true[(i,j)] = pool.apply_sync(idr,
                    (basename_prefix, 
                    args.peaks[i], args.peaks[j], args.peak_pooled,
                    args.idr_thresh, args.idr_rank, args.out_dir))
            elif t=='overlap':
                ret_val_true[(i,j)] = pool.apply_sync(naive_overlap,
                    (basename_prefix, 
                    args.peaks[i], args.peaks[j], args.peak_pooled,
                    args.out_dir))
    
        # for pseudo replicates
        basename_prefix = 
            'IDR.rep{}_pr'.format(i+1)
        if t=='idr':
            ret_val_pr[i] = pool.apply_sync(idr, 
                (basename_prefix, 
                args.peaks_pr1[i], args.peaks_pr2[j], args.peaks[i],
                args.idr_thresh, args.idr_rank, args.out_dir))
        elif t=='overlap'
            ret_val_pr[i] = pool.apply_sync(naive_overlap, 
                (basename_prefix, 
                args.peaks_pr1[i], args.peaks_pr2[j], args.peaks[i],
                args.out_dir))

    # for pooled replicates
    basename_prefix = 'IDR.ppr'
    if t=='idr':
        ret_val_ppr = idr(
            basename_prefix, 
            args.peak_ppr1, args.peak_ppr2, args.peak_pooled,
            args.idr_thresh, args.idr_rank, args.out_dir)
    elif t=='overlap'
        ret_val_ppr = naive_overlap(
            basename_prefix, 
            args.peak_ppr1, args.peak_ppr2, args.peak_pooled,
            args.out_dir)

    # gather processed peaks
    out_true={}
    out_pr={}
    for i in range(len(args.peaks)):
        # for every pair of true replicates
        for j in range(i+1, len(args.peaks)):
            out_true[(i,j)] = ret_val_true[(i,j)].get(BIG_INT)
        # for pseudo replicates
        out_pr[i] = ret_val_pr[i].get(BIG_INT)
    # for pooled replicates
    out_ppr = ret_val_ppr.get(BIG_INT)

    # reproducibility QC
    log.info('Reproducibility QC...')

    # calculate FRiP
    ret_val_true = {}
    for t in args.method: # for each type (IDR and NAIVE-OVERLAP)
        log.info('Cacluating FRiP...')
        for i in range(len(args.peaks)):
            ret_val_truefrip(
            # for every pair of true replicates
            for j in range(i+1, len(args.peaks)):

        # calculate FRiP
        for i in range(len(args.peaks)):
            for j in range(i, len(args.peaks)):
            for j, peak in enumerate(args.peaks):
            args.peaks[i]


    # close multithreading
    pool.close()
    pool.join()

    log.info('All done.')

if __name__=='__main__':
    main()