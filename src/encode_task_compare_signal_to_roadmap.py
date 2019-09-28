#!/usr/bin/env python

# ENCODE DCC compare signal to roadmap wrapper
# Author: Daniel Kim, Jin Lee (leepc12@gmail.com)

import warnings
from matplotlib import pyplot as plt
import sys
import os
import argparse
from encode_lib_common import (
    strip_ext_bigwig, ls_l, log, mkdir_p)
import numpy as np
import pandas as pd
import scipy.stats
import matplotlib as mpl
mpl.use('Agg')

warnings.filterwarnings("ignore")


def parse_arguments():
    parser = argparse.ArgumentParser(prog='ENCODE compare signal to roadmap')
    parser.add_argument('--bigwig', type=str,
                        help='BIGWIG file (from task macs2).')
    parser.add_argument('--dnase', type=str, help='DNase file.')
    parser.add_argument('--reg2map', type=str, help='Reg2map file.')
    parser.add_argument('--reg2map-bed', type=str, help='Reg2map bed file.')
    parser.add_argument('--roadmap-meta', type=str,
                        help='Roadmap metadata file.')
    parser.add_argument('--out-dir', default='', type=str,
                        help='Output directory.')
    parser.add_argument('--log-level', default='INFO', help='Log level',
                        choices=['NOTSET', 'DEBUG', 'INFO', 'WARNING',
                                 'CRITICAL', 'ERROR', 'CRITICAL'])
    args = parser.parse_args()
    log.setLevel(args.log_level)
    log.info(sys.argv)
    return args


def compare_to_roadmap(bw_file, regions_file, reg2map_file,
                       metadata, output_prefix):
    '''
    Takes a bigwig file and signal file, gets the bwAverageOverBed,
    then compares that signal with the signal in the Roadmap
    regions
    '''

    out_file = '{0}.signal'.format(output_prefix)
    log_file = '{0}.roadmap_compare.log'.format(output_prefix)

    # First get the signal vals for the peak regions
    # remember to use a UCSC formatted bed file for regions
    bw_average_over_bed = 'bigWigAverageOverBed {0} {1} {2}'.format(
        bw_file, regions_file, out_file)
    log.info(bw_average_over_bed)
    os.system(bw_average_over_bed)

    # Read the file back in
    sample_data = pd.read_table(out_file, header=None)
    sample_mean0_col = np.array(sample_data.iloc[:, 5])

    # Then, calculate correlations with all other Roadmap samples and rank
    # the correlations
    roadmap_signals = pd.read_table(reg2map_file, compression='gzip')
    (nrow, ncol) = roadmap_signals.shape

    results = pd.DataFrame(columns=('eid', 'corr'))
    with open(log_file, 'w') as fp:
        for i in range(ncol):
            # Slice, run correlation
            roadmap_i = roadmap_signals.iloc[:, i]
            spearman_corr = scipy.stats.spearmanr(np.array(roadmap_i),
                                                  sample_mean0_col)
            results.loc[i] = [roadmap_i.name, spearman_corr[0]]
            s = '{0}\t{1}'.format(roadmap_i.name, spearman_corr)
            log.info(s)
            fp.write(s + '\n')

    # Read in metadata to make the chart more understandable
    metadata = pd.read_table(metadata)
    metadata.columns = ['eid', 'mnemonic']

    merged = pd.merge(metadata, results, on='eid')

    sorted_results = merged.sort_values('corr', ascending=True)

    # Plot results
    pos = np.array(range(ncol)) + 0.5
    fig = plt.figure(figsize=(5, int(ncol/4)))
    plt.barh(pos, sorted_results['corr'], align='center', height=1.0)
    plt.yticks(pos, sorted_results['mnemonic'].tolist(), fontsize=7)
    plt.xlabel('Spearmans correlation')
    plt.title('Signal correlation to Roadmap DNase')
    plt.axis('tight')
    ax = plt.axes()
    ax.yaxis.set_ticks_position('none')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    plot_img = output_prefix + '.roadmap_compare_plot.png'
    fig.savefig(plot_img, format='png', bbox_inches='tight')

    return plot_img


def main():
    # read params
    args = parse_arguments()
    BIGWIG = args.bigwig
    DNASE = args.dnase
    OUTPUT_PREFIX = os.path.join(
        args.out_dir,
        os.path.basename(strip_ext_bigwig(BIGWIG)))

    REG2MAP_BED = args.reg2map_bed if args.reg2map_bed and os.path.basename(
        args.reg2map_bed) != 'null' else DNASE
    REG2MAP = args.reg2map if args.reg2map and os.path.basename(
        args.reg2map) != 'null' else ''
    ROADMAP_META = args.roadmap_meta if args.roadmap_meta and os.path.basename(
        args.roadmap_meta) != 'null' else ''

    log.info('Initializing and making output directory...')
    mkdir_p(args.out_dir)

    compare_to_roadmap(BIGWIG, REG2MAP_BED, REG2MAP,
                       ROADMAP_META, OUTPUT_PREFIX)

    log.info('List all files in output directory...')
    ls_l(args.out_dir)

    log.info('All done.')


if __name__ == '__main__':
    main()
