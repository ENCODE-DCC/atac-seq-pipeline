#!/usr/bin/env python

# ENCODE DCC reproducibility QC wrapper
# Author: Jin Lee (leepc12@gmail.com)

import sys
import os
import argparse
import multiprocessing
from encode_common import *
from encode_common_genomic import peak_to_bigbed

def parse_arguments():
    parser = argparse.ArgumentParser(prog='ENCODE DCC reproducibility QC.',
                        description='IDR peak or overlap peak only.')
    parser.add_argument('peaks', type=str, nargs='*',
                        help='List of peak files \
                        from true replicates in a sorted order. \
                        For example of 4 true replicates, \
                         0,1 0,2 0,3 1,2 1,3 2,3. \
                         x,y means peak file from rep-x vs rep-y.')
    parser.add_argument('--peaks-pr', type=str, nargs='+', required=True,
                        help='List of peak files from pseudo replicates.')
    parser.add_argument('--peak-ppr', type=str,
                        help='Peak file from pooled pseudo replicate.')
    parser.add_argument('--peak-type', type=str, default='narrowPeak',
                        choices=['narrowPeak','regionPeak','broadPeak','gappedPeak'],
                        help='Peak file type.')
    parser.add_argument('--chrsz', type=str,
                        help='2-col chromosome sizes file.')
    parser.add_argument('--prefix', type=str,
                        help='Basename prefix for reproducibility QC file.')
    parser.add_argument('--out-dir', default='', type=str,
                        help='Output directory.')
    parser.add_argument('--log-level', default='INFO', 
                        choices=['NOTSET','DEBUG','INFO',
                            'WARNING','CRITICAL','ERROR','CRITICAL'],
                        help='Log level')
    args = parser.parse_args()
    if len(args.peaks_pr)!=infer_n_from_nC2(len(args.peaks)):
        raise argparse.ArgumentTypeError(
            'Invalid number of peak files or --peak-pr.')

    log.setLevel(args.log_level)
    log.info(sys.argv)
    return args

def main():
    # read params
    args = parse_arguments()
    log.info('Initializing and making output directory...')
    mkdir_p(args.out_dir)

    log.info('Reproducibility QC...')
    # description for variables
    # N: list of number of peaks in peak files from pseudo replicates
    # Nt: top number of peaks in peak files from true replicates (rep-x_vs_rep-y)
    # Np: number of peaks in peak files from pooled pseudo replicate
    N = [get_num_lines(peak) for peak in args.peaks_pr]
    if len(args.peaks):
        # multiple replicate case
        num_rep = infer_n_from_nC2(len(args.peaks))
        num_peaks_tr = [get_num_lines(peak) for peak in args.peaks]

        Nt = max(num_peaks_tr)
        Np = get_num_lines(args.peak_ppr)
        rescue_ratio = float(max(Np,Nt))/float(min(Np,Nt))
        self_consistency_ratio = float(max(N))/float(min(N))

        Nt_idx = num_peaks_tr.index(Nt)
        label_tr = infer_pair_label_from_idx(num_rep, Nt_idx)

        conservative_set = label_tr
        conservative_peak = make_hard_link(args.peaks[Nt_idx], args.out_dir)
        N_conservative = Nt
        if Nt>Np:
            optimal_set = conservative_set
            optimal_peak = conservative_peak
            N_optimal = N_conservative
        else:
            optimal_set = "ppr"
            optimal_peak = make_hard_link(args.peak_ppr, args.out_dir)
            N_optimal = Np
    else:
        # single replicate case
        num_rep = 1
        
        Nt = 0
        Np = 0
        rescue_ratio = 0.0
        self_consistency_ratio = 1.0

        conservative_set = 'rep1-pr'
        conservative_peak = make_hard_link(args.peaks_pr[0], args.out_dir)
        N_conservative = N[0]
        optimal_set = conservative_set
        optimal_peak = conservative_peak
        N_optimal = N_conservative

    reproducibility = 'pass'
    if rescue_ratio>2.0 or self_consistency_ratio>2.0:
        reproducibility = 'borderline'
    if rescue_ratio>2.0 and self_consistency_ratio>2.0:
        reproducibility = 'fail'

    log.info('Writing optimal/conservative peak files...')
    optimal_peak_file = os.path.join(args.out_dir, 'optimal_peak.gz')
    conservative_peak_file = os.path.join(args.out_dir, 'conservative_peak.gz')
    copy_f_to_f(optimal_peak, optimal_peak_file)
    copy_f_to_f(conservative_peak, conservative_peak_file)

    log.info('Converting peak to bigbed...')
    if args.chrsz:
        peak_to_bigbed(optimal_peak_file, args.peak_type, args.chrsz, args.out_dir)
        peak_to_bigbed(conservative_peak_file, args.peak_type, args.chrsz, args.out_dir)

    log.info('Writing reproducibility QC log...')
    if args.prefix:
        reproducibility_qc = '{}.reproducibility.qc'.format(args.prefix)
    else:
        reproducibility_qc = 'reproducibility.qc'
    reproducibility_qc = os.path.join(args.out_dir, reproducibility_qc)
    
    with open(reproducibility_qc,'w') as fp:
        header = '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
            'Nt',
            '\t'.join(['N{}'.format(i+1) for i in range(num_rep)]),
            'Np',
            'N_opt',
            'N_consv',
            'opt_set',
            'consv_set',
            'rescue_ratio',
            'self_consistency_ratio',
            'reproducibility',
            )
        line = '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
            Nt,
            '\t'.join([str(i) for i in N]),
            Np,
            N_optimal,
            N_conservative,
            optimal_set,
            conservative_set,
            rescue_ratio,
            self_consistency_ratio,
            reproducibility)
        fp.write(header)
        fp.write(line)

    log.info('All done.')

if __name__=='__main__':
    main()