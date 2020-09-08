#!/usr/bin/env python

# ENCODE preseq wrapper
# Author: Daniel Kim, Jin Lee (leepc12@gmail.com)

import warnings
import numpy as np
from matplotlib import pyplot as plt
import sys
import os
import argparse
from encode_lib_common import (
    strip_ext_bam, ls_l, log, logging, rm_f)
from encode_lib_genomic import (
    remove_read_group, locate_picard, samtools_sort)
import matplotlib as mpl
mpl.use('Agg')

warnings.filterwarnings("ignore")


def parse_arguments():
    parser = argparse.ArgumentParser(prog='ENCODE preseq')
    parser.add_argument('--paired-end', action="store_true",
                        help='Paired-end BAM.')
    parser.add_argument('--bam', type=str, help='Raw BAM file.')
    parser.add_argument('--picard-java-heap',
                        help='Picard\'s Java max. heap: java -jar picard.jar '
                             '-Xmx[MAX_HEAP]')
    parser.add_argument('--nth', type=int, default=1,
                        help='Number of threads to parallelize.')
    parser.add_argument('--mem-gb', type=float,
                        help='Max. memory for samtools sort in GB. '
                        'It should be total memory for this task (not memory per thread).')
    parser.add_argument('--out-dir', default='', type=str,
                        help='Output directory.')
    parser.add_argument('--log-level', default='INFO', help='Log level',
                        choices=['NOTSET', 'DEBUG', 'INFO', 'WARNING',
                                 'CRITICAL', 'ERROR', 'CRITICAL'])
    args = parser.parse_args()
    log.setLevel(args.log_level)
    log.info(sys.argv)
    return args


def get_picard_complexity_metrics(aligned_bam, prefix, java_heap=None):
    '''
    Picard EsimateLibraryComplexity
    '''
    # remove redundant (or malformed) info (read group) from bam
    out_file = '{0}.picardcomplexity.qc'.format(prefix)
    if java_heap is None:
        java_heap_param = '-Xmx6G'
    else:
        java_heap_param = '-Xmx{}'.format(java_heap)
    get_gc_metrics = (
        'mkdir -p tmp_java && java -Djava.io.tmpdir=$PWD/tmp_java '
        '{3} -XX:ParallelGCThreads=1 -jar '
        '{2} '
        'EstimateLibraryComplexity INPUT={0} OUTPUT={1} '
        'USE_JDK_DEFLATER=TRUE USE_JDK_INFLATER=TRUE '
        'VERBOSITY=ERROR '
        'QUIET=TRUE && rm -rf tmp_java').format(
        aligned_bam, out_file, locate_picard(), java_heap_param)
    os.system(get_gc_metrics)

    # Extract the actual estimated library size
    header_seen = False
    est_library_size = 0
    with open(out_file, 'r') as fp:
        for line in fp:
            if header_seen:
                est_library_size = int(float(line.strip().split()[-1]))
                break
            if 'ESTIMATED_LIBRARY_SIZE' in line:
                header_seen = True

    return est_library_size


def run_preseq(bam_w_dups, prefix, nth=1, mem_gb=None):
    '''
    Runs preseq. Look at preseq data output to get PBC/NRF.
    '''
    # First sort because this file no longer exists...

    sort_bam = samtools_sort(bam_w_dups, nth, mem_gb)

    logging.info('Running preseq...')
    preseq_data = '{0}.preseq.dat'.format(prefix)
    preseq_log = '{0}.preseq.log'.format(prefix)

    run_shell_cmd(
        'preseq lc_extrap -P -B -o {preseq_data} {sort_bam} '
        '-seed 1 -v 2> {preseq_log}'.format(
            preseq_data=preseq_data,
            sort_bam=sort_bam,
            preseq_log=preseq_log,
        )
    )
    rm_f(sort_bam)

    return preseq_data, preseq_log


def get_preseq_plot(data_file, prefix):
    '''
    Generate a preseq plot
    '''
    try:
        data = np.loadtxt(data_file, skiprows=1)
    except IOError:
        return ''
    data /= 1e6  # scale to millions of reads

    fig = plt.figure()

    # Plot the average expected yield
    plt.plot(data[:, 0], data[:, 1], 'r-')

    # Plot confidence intervals
    ci_lower, = plt.plot(data[:, 0], data[:, 2], 'b--')
    ci_upper, = plt.plot(data[:, 0], data[:, 3], 'b--')
    plt.legend([ci_lower], ['95% confidence interval'], loc=4)

    plt.title('Preseq estimated yield')
    plt.xlabel('Sequenced fragments [ millions ]')
    plt.ylabel('Expected distinct fragments [ millions ]')

    # plot_img = BytesIO()
    # fig.savefig(plot_img, format='png')
    plot_img = prefix + '.preseq.png'
    fig.savefig(plot_img, format='png')

    return plot_img


def main():
    # read params
    args = parse_arguments()

    ALIGNED_BAM = args.bam
    OUTPUT_PREFIX = os.path.join(
        args.out_dir,
        os.path.basename(strip_ext_bam(ALIGNED_BAM)))
    RG_FREE_ALIGNED_BAM = remove_read_group(ALIGNED_BAM)
    JAVA_HEAP = args.picard_java_heap
    # Library complexity: Preseq results, NRF, PBC1, PBC2
    if args.paired_end:
        picard_est_lib_size = get_picard_complexity_metrics(
            RG_FREE_ALIGNED_BAM, OUTPUT_PREFIX, JAVA_HEAP)
    else:
        picard_est_lib_size = None
    preseq_data, preseq_log = run_preseq(
        ALIGNED_BAM, OUTPUT_PREFIX, args.nth, args.mem_gb)  # SORTED BAM

    get_preseq_plot(preseq_data, OUTPUT_PREFIX)

    # write picard_est_lib_size to file
    if picard_est_lib_size is not None:
        picard_est_lib_size_file = OUTPUT_PREFIX + '.picard_est_lib_size.qc'
        with open(picard_est_lib_size_file, 'w') as fp:
            fp.write(str(picard_est_lib_size)+'\n')

    rm_f(RG_FREE_ALIGNED_BAM)

    log.info('List all files in output directory...')
    ls_l(args.out_dir)

    log.info('All done.')


if __name__ == '__main__':
    main()
