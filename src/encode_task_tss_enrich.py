#!/usr/bin/env python

# ENCODE DCC TSS enrich wrapper
# Author: Daniel Kim, Jin Lee (leepc12@gmail.com)

import pybedtools
import numpy as np
from matplotlib import mlab
from matplotlib import pyplot as plt
import sys
import os
import argparse
from encode_lib_common import (
    strip_ext_bam, ls_l, log, logging, mkdir_p, rm_f)
from encode_lib_genomic import (
    remove_read_group, samtools_index)
import matplotlib as mpl
mpl.use('Agg')
import metaseq
import warnings
warnings.filterwarnings("ignore")


def parse_arguments():
    parser = argparse.ArgumentParser(prog='ENCODE TSS enrichment.')
    parser.add_argument('--read-len-log', type=str,
                        help='Read length log file (from aligner task).')
    parser.add_argument('--nodup-bam', type=str,
                        help='Raw BAM file (from task filter).')
    parser.add_argument('--chrsz', type=str,
                        help='2-col chromosome sizes file.')
    parser.add_argument('--tss', type=str, help='TSS definition bed file.')
    parser.add_argument('--out-dir', default='', type=str,
                        help='Output directory.')
    parser.add_argument('--log-level', default='INFO', help='Log level',
                        choices=['NOTSET', 'DEBUG', 'INFO', 'WARNING',
                                 'CRITICAL', 'ERROR', 'CRITICAL'])
    args = parser.parse_args()
    log.setLevel(args.log_level)
    log.info(sys.argv)
    return args


def make_tss_plot(bam_file, tss, prefix, chromsizes,
                  read_len, bins=400, bp_edge=2000,
                  processes=8, greenleaf_norm=True):
    '''
    Take bootstraps, generate tss plots, and get a mean and
    standard deviation on the plot. Produces 2 plots. One is the
    aggregation plot alone, while the other also shows the signal
    at each TSS ordered by strength.
    '''
    logging.info('Generating tss plot...')
    tss_plot_file = '{0}.tss_enrich.png'.format(prefix)
    tss_plot_large_file = '{0}.large_tss_enrich.png'.format(prefix)
    tss_log_file = '{0}.tss_enrich.qc'.format(prefix)

    # Load the TSS file
    tss = pybedtools.BedTool(tss)
    tss_ext = tss.slop(b=bp_edge, g=chromsizes)

    # Load the bam file
    # Need to shift reads and just get ends, just load bed file?
    bam = metaseq.genomic_signal(bam_file, 'bam')
    # Shift to center the read on the cut site
    bam_array = bam.array(tss_ext, bins=bins, shift_width=-read_len/2,
                          processes=processes, stranded=True)
    # Normalization (Greenleaf style): Find the avg height
    # at the end bins and take fold change over that
    if greenleaf_norm:
        # Use enough bins to cover 100 bp on either end
        num_edge_bins = int(100/(2*bp_edge/bins))
        bin_means = bam_array.mean(axis=0)
        avg_noise = (sum(bin_means[:num_edge_bins]) +
                     sum(bin_means[-num_edge_bins:]))/(2*num_edge_bins)
        bam_array /= avg_noise
    else:
        bam_array /= bam.mapped_read_count() / 1e6

    # Generate a line plot
    fig = plt.figure()
    ax = fig.add_subplot(111)
    x = np.linspace(-bp_edge, bp_edge, bins)

    ax.plot(x, bam_array.mean(axis=0), color='r', label='Mean')
    ax.axvline(0, linestyle=':', color='k')

    # Note the middle high point (TSS)
    tss_point_val = max(bam_array.mean(axis=0))

    # write tss_point_val to file
    with open(tss_log_file, 'w') as fp:
        fp.write(str(tss_point_val))

    ax.set_xlabel('Distance from TSS (bp)')
    if greenleaf_norm:
        ax.set_ylabel('TSS Enrichment')
    else:
        ax.set_ylabel('Average read coverage (per million mapped reads)')
    ax.legend(loc='best')

    fig.savefig(tss_plot_file)

    # Print a more complicated plot with lots of info

    # Find a safe upper percentile - we can't use X if the Xth percentile is 0
    upper_prct = 99
    if mlab.prctile(bam_array.ravel(), upper_prct) == 0.0:
        upper_prct = 100.0

    plt.rcParams['font.size'] = 8
    fig = metaseq.plotutils.imshow(bam_array,
                                   x=x,
                                   figsize=(5, 10),
                                   vmin=5, vmax=upper_prct, percentile=True,
                                   line_kwargs=dict(color='k', label='All'),
                                   fill_kwargs=dict(color='k', alpha=0.3),
                                   sort_by=bam_array.mean(axis=1))

    # And save the file
    fig.savefig(tss_plot_large_file)

    return tss_plot_file, tss_plot_large_file, tss_log_file


def main():
    # read params
    args = parse_arguments()

    CHROMSIZES = args.chrsz
    TSS = args.tss if args.tss and os.path.basename(args.tss) != 'null' else ''
    FINAL_BAM = args.nodup_bam
    OUTPUT_PREFIX = os.path.join(
        args.out_dir,
        os.path.basename(strip_ext_bam(FINAL_BAM)))
    samtools_index(FINAL_BAM)  # make an index first
    RG_FREE_FINAL_BAM = remove_read_group(FINAL_BAM)

    log.info('Initializing and making output directory...')
    mkdir_p(args.out_dir)

    # Also get read length
    # read_len = get_read_length(FASTQ)
    if args.read_len_log:
        with open(args.read_len_log, 'r') as fp:
            read_len = int(fp.read().strip())
    else:
        read_len = None

    # Enrichments: V plot for enrichment
    # Use final to avoid duplicates
    tss_plot, tss_large_plot, tss_enrich_qc = \
        make_tss_plot(FINAL_BAM,
                      TSS,
                      OUTPUT_PREFIX,
                      CHROMSIZES,
                      read_len)

    # remove temporary files
    rm_f(RG_FREE_FINAL_BAM)

    log.info('List all files in output directory...')
    ls_l(args.out_dir)

    log.info('All done.')


if __name__ == '__main__':
    main()
