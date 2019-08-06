#!/usr/bin/env python

# ENCODE DCC GC bias wrapper
# Author: Daniel Kim, Jin Lee (leepc12@gmail.com)

import sys
import os
import argparse
from encode_lib_common import strip_ext_bam, ls_l, log, logging, mkdir_p, rm_f
from encode_lib_genomic import remove_read_group, locate_picard
import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use('Agg')

from matplotlib import pyplot as plt
import warnings
warnings.filterwarnings("ignore")

def parse_arguments():
    parser = argparse.ArgumentParser(prog='ENCODE GC bias')
    parser.add_argument('--read-len-log', type=str, help='Read length log file (from aligner task).')
    parser.add_argument('--nodup-bam', type=str, help='Raw BAM file (from task filter).')
    parser.add_argument('--ref-fa', type=str, help='Reference fasta file.')
    parser.add_argument('--out-dir', default='', type=str, help='Output directory.')
    parser.add_argument('--log-level', default='INFO', help='Log level',
                        choices=['NOTSET','DEBUG','INFO','WARNING','CRITICAL','ERROR','CRITICAL'])
    args = parser.parse_args()
    log.setLevel(args.log_level)
    log.info(sys.argv)
    return args

def get_gc(qsorted_bam_file, reference_fasta, prefix):
    '''
    Uses picard tools (CollectGcBiasMetrics). Note that the reference
    MUST be the same fasta file that generated the bowtie indices.
    Assumes picard was already loaded into space (module add picard-tools)
    '''
    # remove redundant (or malformed) info (read group) from bam
    logging.info('Getting GC bias...')
    output_file = '{0}.gc.txt'.format(prefix)
    plot_file = '{0}.gcPlot.pdf'.format(prefix)
    summary_file = '{0}.gcSummary.txt'.format(prefix)
    get_gc_metrics = ('java -Xmx6G -XX:ParallelGCThreads=1 -jar '
                      '{5} '
                      'CollectGcBiasMetrics R={0} I={1} O={2} '
                      'USE_JDK_DEFLATER=TRUE USE_JDK_INFLATER=TRUE '
                      'VERBOSITY=ERROR QUIET=TRUE '
                      'ASSUME_SORTED=FALSE '
                      'CHART={3} S={4}').format(reference_fasta,
                                                qsorted_bam_file,
                                                output_file,
                                                plot_file,
                                                summary_file,
                                                locate_picard())
    logging.info(get_gc_metrics)
    os.system(get_gc_metrics)
    return output_file, plot_file, summary_file

def plot_gc(data_file, prefix):
    '''
    Replot the Picard output as png file to put into the html
    '''
    # Load data
    data = pd.read_table(data_file, comment="#")

    # Plot the data
    fig = plt.figure()
    ax = fig.add_subplot(111)

    plt.xlim((0, 100))

    lin1 = ax.plot(data['GC'], data['NORMALIZED_COVERAGE'],
                   label='Normalized coverage', color='r')
    ax.set_ylabel('Normalized coverage')

    ax2 = ax.twinx()
    lin2 = ax2.plot(data['GC'], data['MEAN_BASE_QUALITY'],
                    label='Mean base quality at GC%', color='b')
    ax2.set_ylabel('Mean base quality at GC%')

    ax3 = ax.twinx()
    lin3 = ax3.plot(data['GC'], data['WINDOWS']/np.sum(data['WINDOWS']),
                    label='Windows at GC%', color='g')
    ax3.get_yaxis().set_visible(False)

    lns = lin1 + lin2 + lin3
    labs = [l.get_label() for l in lns]
    ax.legend(lns, labs, loc='best')

    # plot_img = BytesIO()
    # fig.savefig(plot_img, format='png')
    prefix = data_file.rstrip('.gc.txt')
    plot_png = prefix + '.gc_plot.png'
    fig.savefig(plot_png, format='png')

    return plot_png

def main():
    # read params
    args = parse_arguments()

    REF = args.ref_fa
    FINAL_BAM = args.nodup_bam
    OUTPUT_PREFIX = os.path.join(
        args.out_dir,
        os.path.basename(strip_ext_bam(FINAL_BAM)))
    RG_FREE_FINAL_BAM = remove_read_group(FINAL_BAM)

    gc_out, gc_plot_pdf, gc_summary = get_gc(RG_FREE_FINAL_BAM,
                                         REF,
                                         OUTPUT_PREFIX)
    # will generate PNG format from gc_out
    gc_plot_png = plot_gc(gc_out, OUTPUT_PREFIX)

    rm_f(RG_FREE_FINAL_BAM)

    log.info('List all files in output directory...')
    ls_l(args.out_dir)

    log.info('All done.')

if __name__=='__main__':
    main()
