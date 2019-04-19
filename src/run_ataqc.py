#!/usr/bin/env python2

# Daniel Kim, CS Foo
# 2016-03-28
# Script to run ataqc, all parts

import matplotlib
matplotlib.use('Agg')

import os
import sys
import pysam
import pybedtools
import metaseq
import subprocess
import multiprocessing
import timeit
import datetime
import gzip
import numpy as np
import pandas as pd
import scipy.stats
import argparse
import logging
import re
import signal

from base64 import b64encode
from collections import namedtuple
from collections import OrderedDict
from io import BytesIO
from scipy.signal import find_peaks_cwt
from jinja2 import Template
from matplotlib import pyplot as plt
from matplotlib import mlab
from encode_common_genomic import *


# utils

def run_shell_cmd(cmd):
    """Taken from ENCODE DCC ATAC pipeline
    """
    p = subprocess.Popen(['/bin/bash','-o','pipefail'], # to catch error in pipe
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        universal_newlines=True,
        preexec_fn=os.setsid) # to make a new process with a new PGID
    pid = p.pid
    pgid = os.getpgid(pid)
    log.info('run_shell_cmd: PID={}, PGID={}, CMD={}'.format(pid, pgid, cmd))
    stdout, stderr = p.communicate(cmd)
    rc = p.returncode
    err_str = 'PID={}, PGID={}, RC={}\nSTDERR={}\nSTDOUT={}'.format(pid, pgid, rc,
        stderr.strip(), stdout.strip())
    if rc:
        # kill all child processes
        try:
            os.killpg(pgid, signal.SIGKILL)
        except:
            pass
        finally:
            raise Exception(err_str)
    else:
        log.info(err_str)
    return stdout.strip('\n')

def get_num_lines(f):
    """Taken from ENCODE DCC ATAC pipeline
    """
    cmd = 'zcat -f {} | wc -l'.format(f)
    return int(run_shell_cmd(cmd))


# QC STUFF

QCResult = namedtuple('QCResult', ['metric', 'qc_pass', 'message'])
INF = float("inf")


class QCCheck(object):
    def __init__(self, metric):
        self.metric = metric

    def check(self, value):
        return True

    def message(self, value, qc_pass):
        return ('{} - OK'.format(value) if qc_pass
                else '{} - Failed'.format(value))

    def __call__(self, value):
        qc_pass = self.check(value)
        return QCResult(self.metric, qc_pass, self.message(value, qc_pass))


class QCIntervalCheck(QCCheck):
    def __init__(self, metric, lower, upper):
        super(QCIntervalCheck, self).__init__(metric)
        self.lower = lower
        self.upper = upper

    def check(self, value):
        return self.lower <= value <= self.upper

    def message(self, value, qc_pass):
        return ('{} - OK'.format(value) if qc_pass else
                '{} out of range [{}, {}]'.format(value, self.lower,
                                                  self.upper))


class QCLessThanEqualCheck(QCIntervalCheck):
    def __init__(self, metric, upper):
        super(QCLessThanEqualCheck, self).__init__(metric, -INF, upper)


class QCGreaterThanEqualCheck(QCIntervalCheck):
    def __init__(self, metric, lower):
        super(QCGreaterThanEqualCheck, self).__init__(metric, lower, INF)


class QCHasElementInRange(QCCheck):
    def __init__(self, metric, lower, upper):
        super(QCHasElementInRange, self).__init__(metric)
        self.lower = lower
        self.upper = upper

    def check(self, elems):
        return (len([elem for elem in elems
                    if self.lower <= elem <= self.upper]) > 0)

    def message(self, elems, qc_pass):
        return ('OK' if qc_pass else
                'Cannot find element in range [{}, {}]'.format(
                    self.lower, self.upper))

# HELPER FUNCTIONS

def getFileHandle(filename, mode="r"):
    if (re.search('.gz$',filename) or re.search('.gzip',filename)):
        if (mode=="r"):
            mode="rb";
        return gzip.open(filename,mode)
    else:
        return open(filename,mode)


# QC FUNCTIONS

def determine_paired(bam_file):
    '''
    Quick function to determine if the BAM file is paired end or single end
    '''
    num_paired_reads = int(subprocess.check_output(['samtools',
                                                    'view', '-f', '0x1',
                                                    '-c', bam_file]).strip())
    if num_paired_reads > 1:
        return "Paired-ended"
    else:
        return "Single-ended"


def get_read_length(fastq_file):
    '''
    Get read length out of fastq file
    '''
    total_reads_to_consider = 1000000
    line_num = 0
    total_reads_considered = 0
    max_length = 0
    with getFileHandle(fastq_file, 'rb') as fp:
        for line in fp:
            if line_num % 4 == 1:
                if len(line.strip()) > max_length:
                    max_length = len(line.strip())
                total_reads_considered += 1
            if total_reads_considered >= total_reads_to_consider:
                break
            line_num += 1

    return int(max_length)


def get_bowtie_stats(bowtie_alignment_log):
    '''
    From the Bowtie alignment log, get relevant stats and return
    the file in a list format where each line is an element in
    the list. Can be parsed further if desired.
    '''
    logging.info('Reading bowtie alignment log...')
    bowtie_text = ''
    with open(bowtie_alignment_log, 'rb') as fp:
        for line in fp:
            logging.info(line.strip())
            bowtie_text += line
    return bowtie_text


def get_chr_m(sorted_bam_file, mito_chr_name):
    '''
    Get fraction of reads that are mitochondrial (chr M).
    '''
    logging.info('Getting mitochondrial chromosome fraction...')
    chrom_list = pysam.idxstats(sorted_bam_file, split_lines=True)
    tot_reads = 0
    chr_m_reads = 0
    for chrom in chrom_list:
        chrom_stats = chrom.split('\t')
        if chrom_stats[0] == mito_chr_name:
            chr_m_reads = int(chrom_stats[2])
        tot_reads += int(chrom_stats[2])
    if tot_reads==0:
        fract_chr_m = 0
    else:
        fract_chr_m = float(chr_m_reads) / tot_reads

    return chr_m_reads, fract_chr_m

def get_gc(qsorted_bam_file, reference_fasta, prefix):
    '''
    Uses picard tools (CollectGcBiasMetrics). Note that the reference
    MUST be the same fasta file that generated the bowtie indices.
    Assumes picard was already loaded into space (module add picard-tools)
    '''
    # remove redundant (or malformed) info (read group) from bam
    logging.info('Getting GC bias...')
    output_file = '{0}_gc.txt'.format(prefix)
    plot_file = '{0}_gcPlot.pdf'.format(prefix)
    summary_file = '{0}_gcSummary.txt'.format(prefix)
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


def plot_gc(data_file):
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

    plot_img = BytesIO()
    fig.savefig(plot_img, format='png')

    return b64encode(plot_img.getvalue())


def run_preseq(bam_w_dups, prefix):
    '''
    Runs preseq. Look at preseq data output to get PBC/NRF.
    '''
    # First sort because this file no longer exists...
    sort_bam = 'samtools sort -o {1}.sorted.bam -T {1} -@ 2 {0}'.format(
        bam_w_dups, prefix)
    os.system(sort_bam)

    logging.info('Running preseq...')
    preseq_data = '{0}.preseq.dat'.format(prefix)
    preseq_log = '{0}.preseq.log'.format(prefix)
    preseq = ('preseq lc_extrap '
              '-P -B -o {0} {1}.sorted.bam -seed 1 -v 2> {2}').format(preseq_data,
                                                              prefix,
                                                              preseq_log)
    logging.info(preseq)
    os.system(preseq)
    os.system('rm {0}.sorted.bam'.format(prefix))
    return preseq_data, preseq_log


def get_encode_complexity_measures(pbc_output):
    '''
    Gets the unique count statistics from the filtered bam file,
    which is consistent with ENCODE metrics for ChIP-seq
    '''
    with open(pbc_output, 'rb') as fp:
        for line in fp:
            l_list = line.strip().split('\t')
            NRF = float(l_list[4])
            PBC1 = float(l_list[5])
            PBC2 = float(l_list[6])
            break

    # QC check
    results = []
    results.append(QCGreaterThanEqualCheck('NRF', 0.8)(NRF))
    results.append(QCGreaterThanEqualCheck('PBC1', 0.8)(PBC1))
    results.append(QCGreaterThanEqualCheck('PBC2', 1.0)(PBC2))

    return results


def get_encode_complexity_measures_OLD(preseq_log):
    '''
    Use info from the preseq log to calculate NRF, PBC1, and PBC2
    '''
    with open(preseq_log, 'rb') as fp:
        for line in fp:
            if line.startswith('TOTAL READS'):
                tot_reads = float(line.strip().split("= ")[1])
            elif line.startswith('DISTINCT READS'):
                distinct_reads = float(line.strip().split('= ')[1])
            elif line.startswith('1\t'):
                one_pair = float(line.strip().split()[1])
            elif line.startswith('2\t'):
                two_pair = float(line.strip().split()[1])

    NRF = distinct_reads/tot_reads
    PBC1 = one_pair/distinct_reads
    PBC2 = one_pair/two_pair

    # QC check
    results = []
    results.append(QCGreaterThanEqualCheck('NRF', 0.8)(NRF))
    results.append(QCGreaterThanEqualCheck('PBC1', 0.8)(PBC1))
    results.append(QCGreaterThanEqualCheck('PBC2', 1.0)(PBC2))

    return results


def get_picard_complexity_metrics(aligned_bam, prefix):
    '''
    Picard EsimateLibraryComplexity
    '''
    # remove redundant (or malformed) info (read group) from bam
    out_file = '{0}.picardcomplexity.qc'.format(prefix)
    get_gc_metrics = ('mkdir -p tmp_java && java -Djava.io.tmpdir=$PWD/tmp_java -Xmx6G -XX:ParallelGCThreads=1 -jar '
                      '{2} '
                      'EstimateLibraryComplexity INPUT={0} OUTPUT={1} '
                      'USE_JDK_DEFLATER=TRUE USE_JDK_INFLATER=TRUE '
                      'VERBOSITY=ERROR '
                      'QUIET=TRUE && rm -rf tmp_java').format(aligned_bam,
                                           out_file,
                                           locate_picard())
    os.system(get_gc_metrics)

    # Extract the actual estimated library size
    header_seen = False
    est_library_size = 0
    with open(out_file, 'rb') as fp:
        for line in fp:
            if header_seen:
                est_library_size = int(float(line.strip().split()[-1]))
                break
            if 'ESTIMATED_LIBRARY_SIZE' in line:
                header_seen = True

    return est_library_size


def preseq_plot(data_file):
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

    plot_img = BytesIO()
    fig.savefig(plot_img, format='png')

    return b64encode(plot_img.getvalue())


def make_tss_plot(bam_file, tss, prefix, chromsizes, read_len, bins=400, bp_edge=2000,
                  processes=8, greenleaf_norm=True):
    '''
    Take bootstraps, generate tss plots, and get a mean and
    standard deviation on the plot. Produces 2 plots. One is the
    aggregation plot alone, while the other also shows the signal
    at each TSS ordered by strength.
    '''
    logging.info('Generating tss plot...')
    tss_plot_file = '{0}_tss-enrich.png'.format(prefix)
    tss_plot_large_file = '{0}_large_tss-enrich.png'.format(prefix)

    # Load the TSS file
    tss = pybedtools.BedTool(tss)
    tss_ext = tss.slop(b=bp_edge, g=chromsizes)

    # Load the bam file
    bam = metaseq.genomic_signal(bam_file, 'bam') # Need to shift reads and just get ends, just load bed file?
    bam_array = bam.array(tss_ext, bins=bins, shift_width = -read_len/2, # Shift to center the read on the cut site
                          processes=processes, stranded=True)

    # Actually first build an "ends" file
    #get_ends = '''zcat {0} | awk -F '\t' 'BEGIN {{OFS="\t"}} {{if ($6 == "-") {{$2=$3-1; print}} else {{$3=$2+1; print}} }}' | gzip -c > {1}_ends.bed.gz'''.format(bed_file, prefix)
    #print(get_ends)
    #os.system(get_ends)

    #bed_reads = metaseq.genomic_signal('{0}_ends.bed.gz'.format(prefix), 'bed')
    #bam_array = bed_reads.array(tss_ext, bins=bins,
    #                      processes=processes, stranded=True)

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

    return tss_plot_file, tss_plot_large_file, tss_point_val


def get_picard_dup_stats(picard_dup_file, paired_status):
    '''
    Parse Picard's MarkDuplicates metrics file
    '''
    logging.info('Running Picard MarkDuplicates...')
    mark = 0
    dup_stats = {}
    with open(picard_dup_file) as fp:
        for line in fp:
            if '##' in line:
                if 'METRICS CLASS' in line:
                    mark = 1
                continue

            if mark == 2:
                line_elems = line.strip().split('\t')
                dup_stats['PERCENT_DUPLICATION'] = line_elems[7]
                dup_stats['READ_PAIR_DUPLICATES'] = line_elems[5]
                dup_stats['READ_PAIRS_EXAMINED'] = line_elems[2]
                if paired_status == 'Paired-ended':
                    return float(line_elems[5]), float(line_elems[7])
                else:
                    return float(line_elems[4]), float(line_elems[7])

            if mark > 0:
                mark += 1
    return None


def get_sambamba_dup_stats(sambamba_dup_file, paired_status):
    '''
    Parse sambamba markdup's metrics file
    '''
    logging.info('Running sambamba markdup...')
    with open(sambamba_dup_file, 'r') as fp:
        lines = fp.readlines()

    end_pairs = int(lines[1].strip().split()[1])
    single_ends = int(lines[2].strip().split()[1])
    ends_marked_dup = int(lines[4].strip().split()[1])
    if paired_status == 'Paired-ended':
        pairs_marked_dup = 0.5 * float(ends_marked_dup)
        prct_dup = pairs_marked_dup / float(end_pairs)
        return pairs_marked_dup, prct_dup
    else:
        prct_dup = float(ends_marked_dup) / float(single_ends)
        return ends_marked_dup, prct_dup


def get_mito_dups(sorted_bam, prefix, mito_chr_name, endedness='Paired-ended', use_sambamba=False):
    '''
    Marks duplicates in the original aligned bam file and then determines
    how many reads are duplicates AND from chrM

    To use sambamba markdup instead of picard MarkDuplicates, set
    use_sambamba to True (default False).
    '''

    out_file = '{0}.dupmark.ataqc.bam'.format(prefix)
    metrics_file = '{0}.dup.ataqc'.format(prefix)

    # Filter bam on the flag 0x002
    tmp_filtered_bam = '{0}.filt.bam'.format(prefix)
    tmp_filtered_bam_prefix = tmp_filtered_bam.replace('.bam', '')
    if endedness == 'Paired-ended':
        filter_bam = ('samtools view -F 1804 -f 2 -u {0} | '
                      'samtools sort - {1}'.format(sorted_bam, tmp_filtered_bam_prefix))
    else:
        filter_bam = ('samtools view -F 1804 -u {0} | '
                      'samtools sort - {1}'.format(sorted_bam, tmp_filtered_bam_prefix))
    os.system(filter_bam)

    # Run Picard MarkDuplicates
    mark_duplicates = ('java -Xmx6G -XX:ParallelGCThreads=1 -jar '
                       '{0} '
                       'MarkDuplicates INPUT={1} OUTPUT={2} '
                       'METRICS_FILE={3} '
                       'VALIDATION_STRINGENCY=LENIENT '
                       'ASSUME_SORTED=TRUE '
                       'REMOVE_DUPLICATES=FALSE '
                       'VERBOSITY=ERROR '
                       'USE_JDK_DEFLATER=TRUE USE_JDK_INFLATER=TRUE '
                       'QUIET=TRUE').format(locate_picard(),
                                            tmp_filtered_bam,
                                            out_file,
                                            metrics_file)
    if use_sambamba:
        mark_duplicates = ('sambamba markdup -t 8 '
                           '--hash-table-size=17592186044416 '
                           '--overflow-list-size=20000000 '
                           '--io-buffer-size=256 '
                           '{0} '
                           '{1} '
                           '2> {2}').format(tmp_filtered_bam,
                                            out_file,
                                            metrics_file)
    os.system(mark_duplicates)

    # Index the file
    index_file = 'samtools index {0}'.format(out_file)
    os.system(index_file)

    # Get the mitochondrial reads that are marked duplicates
    mito_dups = int(subprocess.check_output(['samtools',
                                             'view', '-f', '1024',
                                             '-c', out_file, mito_chr_name]).strip())

    total_dups = int(subprocess.check_output(['samtools',
                                              'view', '-f', '1024',
                                              '-c', out_file]).strip())

    # Clean up
    remove_bam = 'rm {0}'.format(out_file)
    os.system(remove_bam)
    remove_metrics_file = 'rm {0}'.format(metrics_file)
    os.system(remove_metrics_file)
    remove_tmp_filtered_bam = 'rm {0}'.format(tmp_filtered_bam)
    os.system(remove_tmp_filtered_bam)

    return mito_dups, float(mito_dups) / total_dups


def get_samtools_flagstat(bam_file):
    '''
    Runs samtools flagstat to get read metrics
    '''
    logging.info('samtools flagstat...')
    results = pysam.flagstat(bam_file, split_lines=True)
    flagstat = ''
    for line in results:
        logging.info(line.strip())
        flagstat += line
        if "mapped" in line and "mate" not in line:
            mapped_reads = int(line.split('+')[0].strip())
    return flagstat, mapped_reads


def get_fract_mapq(bam_file, q=30):
    '''
    Runs samtools view to get the fraction of reads of a certain
    map quality.
    '''
    # Current bug in pysam.view module...
    logging.info('samtools mapq 30...')

    # There is a bug in pysam.view('-c'), so just use subprocess
    num_qreads = int(subprocess.check_output(['samtools',
                                              'view', '-c',
                                              '-q', str(q), bam_file]).strip())
    tot_reads = int(subprocess.check_output(['samtools',
                                             'view', '-c',
                                             bam_file]).strip())
    fract_good_mapq = float(num_qreads)/tot_reads
    return num_qreads, fract_good_mapq


def get_final_read_count(first_bam, last_bam):
    '''
    Get final mapped reads compared to initial reads
    '''
    logging.info('final read counts...')
    # Bug in pysam.view
    num_reads_last_bam = int(subprocess.check_output(['samtools',
                                                      'view', '-c',
                                                      last_bam]).strip())
    num_reads_first_bam = int(subprocess.check_output(['samtools',
                                                       'view', '-c',
                                                       first_bam]).strip())
    fract_reads_left = float(num_reads_last_bam)/num_reads_first_bam

    return num_reads_first_bam, num_reads_last_bam, fract_reads_left


def get_insert_distribution(final_bam, prefix):
    '''
    Calls Picard CollectInsertSizeMetrics
    '''
    logging.info('insert size distribution...')
    insert_data = '{0}.inserts.hist_data.log'.format(prefix)
    insert_plot = '{0}.inserts.hist_graph.pdf'.format(prefix)
    graph_insert_dist = ('java -Xmx6G -XX:ParallelGCThreads=1 -jar '
                         '{3} '
                         'CollectInsertSizeMetrics '
                         'INPUT={0} OUTPUT={1} H={2} '
                         'VERBOSITY=ERROR QUIET=TRUE '
                         'USE_JDK_DEFLATER=TRUE USE_JDK_INFLATER=TRUE '
                         'W=1000 STOP_AFTER=5000000').format(final_bam,
                                                             insert_data,
                                                             insert_plot,
                                                             locate_picard())
    logging.info(graph_insert_dist)
    os.system(graph_insert_dist)
    return insert_data, insert_plot


def get_fract_reads_in_regions_old(reads_bed, regions_bed):
    '''
    Function that takes in bed file of reads and bed file of regions and
    gets fraction of reads sitting in said regions
    '''
    reads_bedtool = pybedtools.BedTool(reads_bed)
    regions_bedtool = pybedtools.BedTool(regions_bed)

    reads = regions_bedtool.sort().merge().intersect(reads_bedtool, c=True, nonamecheck=True)

    read_count = 0
    for interval in reads:
        read_count += int(interval[-1])
    fract_reads = float(read_count)/reads_bedtool.count()

    return read_count, fract_reads


def get_fract_reads_in_regions(reads_bed, regions_bed):
    """Function that takes in bed file of reads and bed file of regions and
    gets fraction of reads sitting in said regions
    """
    # uses new run_shell_cmd
    cmd = "bedtools sort -i {}  | "
    cmd += "bedtools merge -i stdin | "
    cmd += "bedtools intersect -u -nonamecheck -a {} -b stdin | "
    cmd += "wc -l"
    #cmd += "bedtools intersect -c -nonamecheck -a stdin -b {} | "
    #cmd += "awk '{{ sum+=$4 }} END {{ print sum }}'"
    cmd = cmd.format(regions_bed, reads_bed)
    intersect_read_count = int(run_shell_cmd(cmd))
    total_read_count = get_num_lines(reads_bed)
    fract_reads = float(intersect_read_count) / total_read_count

    return intersect_read_count, fract_reads


def get_signal_to_noise(final_bed, dnase_regions, blacklist_regions,
                        prom_regions, enh_regions, peaks):
    '''
    Given region sets, determine whether reads are
    falling in or outside these regions
    '''
    logging.info('signal to noise...')

    # Dnase regions
    reads_dnase, fract_dnase = get_fract_reads_in_regions(final_bed,
                                                          dnase_regions)

    # Blacklist regions
    reads_blacklist, \
        fract_blacklist = get_fract_reads_in_regions(final_bed,
                                                     blacklist_regions)

    # Prom regions
    reads_prom, fract_prom = get_fract_reads_in_regions(final_bed,
                                                        prom_regions)

    # Enh regions
    reads_enh, fract_enh = get_fract_reads_in_regions(final_bed, enh_regions)

    # Peak regions
    reads_peaks, fract_peaks = get_fract_reads_in_regions(final_bed, peaks)

    return reads_dnase, fract_dnase, reads_blacklist, fract_blacklist, \
        reads_prom, fract_prom, reads_enh, fract_enh, reads_peaks, \
        fract_peaks


def get_region_size_metrics(peak_file):
    '''
    From the peak file, return a plot of the region size distribution and
    the quartile metrics (summary from R)
    '''

    peak_size_summ = OrderedDict([
        ('Min size', 0),
        ('25 percentile', 0),
        ('50 percentile (median)', 0),
        ('75 percentile', 0),
        ('Max size', 0),
        ('Mean', 0),
    ])

    # If peak file is none, return nothing
    if peak_file == None:
        return peak_size_summ, ''

    # Load peak file. If it fails, return nothing as above
    try:
        peak_df = pd.read_table(peak_file, compression='gzip', header=None)
    except:
        return peak_size_summ, ''

    # Subtract third column from second to get summary
    region_sizes = peak_df.ix[:,2] - peak_df.ix[:,1]

    # Summarize and store in ordered dict
    peak_summary_stats = region_sizes.describe()

    peak_size_summ = OrderedDict([
        ('Min size', peak_summary_stats['min']),
        ('25 percentile', peak_summary_stats['25%']),
        ('50 percentile (median)', peak_summary_stats['50%']),
        ('75 percentile', peak_summary_stats['75%']),
        ('Max size', peak_summary_stats['max']),
        ('Mean', peak_summary_stats['mean']),
    ])

    # Plot density diagram using matplotlib
    fig = plt.figure()
    ax = fig.add_subplot(111)

    y, binEdges = np.histogram(region_sizes, bins=100)
    bincenters = 0.5 * (binEdges[1:] + binEdges[:-1])

    # density = gaussian_kde(y) # from scipy.stats import gaussian_kde
    # density.covariance_factor = lambda : .25
    # density._compute_covariance()

    plt.plot(bincenters, y, '-')
    filename = peak_file.split('/')[-1]
    ax.set_title('Peak width distribution for {0}'.format(filename))
    #ax.set_yscale('log')

    plot_img = BytesIO()
    fig.savefig(plot_img, format='png')

    return peak_size_summ, b64encode(plot_img.getvalue())


def get_peak_counts(raw_peaks, naive_overlap_peaks=None, idr_peaks=None):
    '''
    Return a table with counts for raw peaks, IDR peaks, and naive
    overlap peaks
    '''

    # Count peaks
    raw_count = sum(1 for line in getFileHandle(raw_peaks))
    if naive_overlap_peaks != None:
        naive_count = sum(1 for line in getFileHandle(naive_overlap_peaks))
    else:
        naive_count = 0

    if idr_peaks != None:
        idr_count = sum(1 for line in getFileHandle(idr_peaks))
    else:
        idr_count = 0

    # Literally just throw these into a QC table
    results = []
    # results.append(QCGreaterThanEqualCheck('Raw peaks', 10000)(raw_count))
    results.append(QCGreaterThanEqualCheck('Naive overlap peaks',
                                           10000)(naive_count))
    results.append(QCGreaterThanEqualCheck('IDR peaks', 10000)(idr_count))

    return results


def track_reads(reads_list, labels):
    '''
    This function takes in read counts for different stages
    and generates a bar chart to show where reads are lost
    '''
    # Initial bam, filters (q30), dups, chrM
    ind = np.arange(len(reads_list))
    width = 0.35

    fig, ax = plt.subplots()
    ax.bar(ind, reads_list, width, color='b')
    ax.set_ylabel('Read count')
    ax.set_title('Reads at each processing step')
    ax.set_xticks(ind+width)
    ax.set_xticklabels(labels)

    plot_img = BytesIO()
    fig.savefig(plot_img, format='png')

    return b64encode(plot_img.getvalue())


def read_picard_histogram(data_file):
    with open(data_file) as fp:
        for line in fp:
            if line.startswith('## HISTOGRAM'):
                break
        data = np.loadtxt(fp, skiprows=1)

    return data


def fragment_length_qc(data):
    results = []

    NFR_UPPER_LIMIT = 150
    MONO_NUC_LOWER_LIMIT = 150
    MONO_NUC_UPPER_LIMIT = 300

    # % of NFR vs res
    nfr_reads = data[data[:,0] < NFR_UPPER_LIMIT][:,1]
    percent_nfr = nfr_reads.sum() / data[:,1].sum()
    results.append(
        QCGreaterThanEqualCheck('Fraction of reads in NFR', 0.4)(percent_nfr))

    # % of NFR vs mononucleosome
    mono_nuc_reads = data[
        (data[:,0] > MONO_NUC_LOWER_LIMIT) &
        (data[:,0] <= MONO_NUC_UPPER_LIMIT)][:,1]
    
    percent_nfr_vs_mono_nuc = (
        nfr_reads.sum() /
        mono_nuc_reads.sum())
    results.append(
        QCGreaterThanEqualCheck('NFR / mono-nuc reads', 2.5)(
            percent_nfr_vs_mono_nuc))

    # peak locations
    pos_start_val = data[0,0] # this may be greater than 0
    peaks = find_peaks_cwt(data[:, 1], np.array([25]))
    nuc_range_metrics = [('Presence of NFR peak', 20 - pos_start_val, 90 - pos_start_val),
                         ('Presence of Mono-Nuc peak', 120 - pos_start_val, 250 - pos_start_val),
                         ('Presence of Di-Nuc peak', 300 - pos_start_val, 500 - pos_start_val)]
    for range_metric in nuc_range_metrics:
        results.append(QCHasElementInRange(*range_metric)(peaks))

    return results


def fragment_length_plot(data_file, peaks=None):
    try:
        data = read_picard_histogram(data_file)
    except IOError:
        return ''
    except TypeError:
        return ''

    fig = plt.figure()
    plt.bar(data[:, 0], data[:, 1])
    plt.xlim((0, 1000))

    if peaks:
        peak_vals = [data[peak_x, 1] for peak_x in peaks]
        plt.plot(peaks, peak_vals, 'ro')

    plot_img = BytesIO()
    fig.savefig(plot_img, format='png')

    return b64encode(plot_img.getvalue())


def compare_to_roadmap(bw_file, regions_file, reg2map_file,
                       metadata, output_prefix):
    '''
    Takes a bigwig file and signal file, gets the bwAverageOverBed,
    then compares that signal with the signal in the Roadmap
    regions
    '''

    out_file = '{0}.signal'.format(output_prefix)

    # First get the signal vals for the peak regions
    # remember to use a UCSC formatted bed file for regions
    bw_average_over_bed = 'bigWigAverageOverBed {0} {1} {2}'.format(
                            bw_file, regions_file, out_file)
    logging.info(bw_average_over_bed)
    os.system(bw_average_over_bed)

    # Read the file back in
    sample_data = pd.read_table(out_file, header=None)
    sample_mean0_col = np.array(sample_data.iloc[:, 5])

    # Then, calculate correlations with all other Roadmap samples and rank
    # the correlations
    roadmap_signals = pd.read_table(reg2map_file, compression='gzip')
    (nrow, ncol) = roadmap_signals.shape

    results = pd.DataFrame(columns=('eid', 'corr'))
    for i in range(ncol):
        # Slice, run correlation
        roadmap_i = roadmap_signals.iloc[:, i]
        spearman_corr = scipy.stats.spearmanr(np.array(roadmap_i),
                                              sample_mean0_col)
        results.loc[i] = [roadmap_i.name, spearman_corr[0]]
        logging.info('{0}\t{1}'.format(roadmap_i.name, spearman_corr))

    # Read in metadata to make the chart more understandable
    metadata = pd.read_table(metadata)
    metadata.columns = ['eid', 'mnemonic']

    merged = pd.merge(metadata, results, on='eid')

    sorted_results = merged.sort_values('corr', ascending=True)

    # Plot results
    pos = np.array(range(ncol)) + 0.5
    fig = plt.figure(figsize=(5,int(ncol/4)))
    plt.barh(pos, sorted_results['corr'], align='center', height=1.0)
    plt.yticks(pos, sorted_results['mnemonic'].tolist(), fontsize=7)
    plt.xlabel('Spearmans correlation')
    plt.title('Signal correlation to Roadmap DNase')
    plt.axis('tight')
    ax = plt.axes()
    ax.yaxis.set_ticks_position('none')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    plot_img = BytesIO()
    fig.savefig(plot_img, format='png', bbox_inches='tight')
    fig.savefig('test.png', format='png', bbox_inches='tight')

    return b64encode(plot_img.getvalue())


html_template = Template("""
{% macro inline_img(base64_img, img_type='png') -%}
    {% if base64_img == '' %}
        <pre>Metric failed.</pre>
    {% else %}
        <img src="data:image/{{ img_type }};base64,{{ base64_img }}">
    {% endif %}
{%- endmacro %}

{% macro qc_table(qc_results) -%}
  <table class='qc_table'>
    <thead>
        <tr>
            <th scope='col'>Metric</th>
            <th scope='col'>Result</th>
        </tr>
    </thead>
    <tbody>
    {% for result in qc_results %}
    <tr>
        <td>{{ result.metric }}</td>
        <td {% if not result.qc_pass %} class='fail' {% endif %}>
          {{ result.message }}
        </td>
    </tr>
    {% endfor %}
    </tbody>
  </table>
{%- endmacro %}

<html>

<head>
  <title>{{ sample['name'] }} - ATAqC report</title>
  <style>
  .qc_table{
      font-family:"Lucida Sans Unicode", "Lucida Grande", Sans-Serif;
      font-size:12px;
      width:480px;
      text-align:left;
      border-collapse:collapse;
      margin:20px;
  }

  .qc_table th{
      font-size:14px;
      font-weight:normal;
      background:#8c1515;
      border-top:4px solid #700000;
      border-bottom:1px solid #fff;
      color:white;
      padding:8px;
  }

  .qc_table td{
      background:#f2f1eb;
      border-bottom:1px solid #fff;
      color:black;
      border-top:1px solid transparent;
      padding:8px;
  }

  .qc_table .fail{
      color:#ff0000;
      font-weight:bold;
  }
  </style>
</head>

<body>
  <h2>ATAqC</h2>

{% if 'basic_info' in sample %}
  <h2>Sample Information</h2>
  <table class='qc_table'>
    <tbody>
      {% for field, value in sample['basic_info'].iteritems() %}
      <tr>
        <td>{{ field }}</td>
        <td>{{ value }}</td>
      </tr>
      {% endfor %}
    </tbody>
  </table>
{% endif %}

{% if 'summary_stats' in sample %}
  <h2>Summary</h2>
  <table class='qc_table'>
    <tbody>
      {% for field, value in sample['summary_stats'].iteritems() %}
      <tr>
        <td>{{ field }}</td>
        <td>{{ '{0:,}'.format(value) }}</td>
      </tr>
      {% endfor %}
    </tbody>
  </table>
<pre>
Note that all these read counts are determined using 'samtools view' - as such,
these are all reads found in the file, whether one end of a pair or a single
end read. In other words, if your file is paired end, then you should divide
these counts by two. Each step follows the previous step; for example, the
duplicate reads were removed after reads were removed for low mapping quality.
</pre>
{% endif %}

{% if 'read_tracker' in sample %}
  {{ inline_img(sample['read_tracker']) }}
<pre>
This bar chart also shows the filtering process and where the reads were lost
over the process. Note that each step is sequential - as such, there may
have been more mitochondrial reads which were already filtered because of
high duplication or low mapping quality. Note that all these read counts are
determined using 'samtools view' - as such, these are all reads found in
the file, whether one end of a pair or a single end read. In other words,
if your file is paired end, then you should divide these counts by two.
</pre>
{% endif %}

{% if 'bowtie_stats' in sample %}
  <h2>Alignment statistics</h2>
  <h3>Bowtie alignment log</h3>
  <pre>
{{ sample['bowtie_stats'] }}
  </pre>
{% endif %}

{% if 'samtools_flagstat' in sample %}
  <h3>Samtools flagstat</h3>
  <pre>
{{ sample['samtools_flagstat'] }}
  </pre>
{% endif %}

{% if 'filtering_stats' in sample %}
  <h2>Filtering statistics</h2>
  <table class='qc_table'>
    <tbody>
      {% for field, value in sample['filtering_stats'].iteritems() %}
      <tr>
        <td>{{ field }}</td>
        <td>{{ '{0:,}'.format(value[0]) }}</td>
        <td>{{ '{0:.3f}'.format(value[1]) }}</td>
      </tr>
      {% endfor %}
    </tbody>
  </table>
  <pre>
Mapping quality refers to the quality of the read being aligned to that
particular location in the genome. A standard quality score is > 30.
Duplications are often due to PCR duplication rather than two unique reads
mapping to the same location. High duplication is an indication of poor
libraries. Mitochondrial reads are often high in chromatin accessibility
assays because the mitochondrial genome is very open. A high mitochondrial
fraction is an indication of poor libraries. Based on prior experience, a
final read fraction above 0.70 is a good library.
  </pre>
{% endif %}

{% if 'encode_lib_complexity' in sample %}
  <h2>Library complexity statistics</h2>
  <h3>ENCODE library complexity metrics</h3>
  {{ qc_table(sample['encode_lib_complexity']) }}
<pre>
The non-redundant fraction (NRF) is the fraction of non-redundant mapped reads
in a dataset; it is the ratio between the number of positions in the genome
that uniquely mapped reads map to and the total number of uniquely mappable
reads. The NRF should be > 0.8. The PBC1 is the ratio of genomic locations
with EXACTLY one read pair over the genomic locations with AT LEAST one read
pair. PBC1 is the primary measure, and the PBC1 should be close to 1.
Provisionally 0-0.5 is severe bottlenecking, 0.5-0.8 is moderate bottlenecking,
0.8-0.9 is mild bottlenecking, and 0.9-1.0 is no bottlenecking. The PBC2 is
the ratio of genomic locations with EXACTLY one read pair over the genomic
locations with EXACTLY two read pairs. The PBC2 should be significantly
greater than 1.
</pre>
{% endif %}

{% if 'picard_est_library_size' in sample %}
  <h3>Picard EstimateLibraryComplexity</h3>
  {{ '{0:,}'.format(sample['picard_est_library_size']) }}
{% endif %}

{% if 'yield_prediction' in sample %}
  <h3>Yield prediction</h3>
  {% if sample['yield_prediction'] == 'Tm9uZQ==' %}
    {{ 'Preseq did not converge (or failed in some other way)'}}
  {% else %}
    {{ inline_img(sample['yield_prediction']) }}
  {% endif %}
<pre>
Preseq performs a yield prediction by subsampling the reads, calculating the
number of distinct reads, and then extrapolating out to see where the
expected number of distinct reads no longer increases. The confidence interval
gives a gauge as to the validity of the yield predictions.
</pre>
{% endif %}

{% if 'fraglen_dist' in sample and 'nucleosomal' in sample %}
  <h2>Fragment length statistics</h2>
  {{ inline_img(sample['fraglen_dist']) }}
  {{ qc_table(sample['nucleosomal']) }}
<pre>
Open chromatin assays show distinct fragment length enrichments, as the cut
sites are only in open chromatin and not in nucleosomes. As such, peaks
representing different n-nucleosomal (ex mono-nucleosomal, di-nucleosomal)
fragment lengths will arise. Good libraries will show these peaks in a
fragment length distribution and will show specific peak ratios.
</pre>

{% endif %}

  <h2>Peak statistics</h2>

{% if 'peak_counts' in sample %}
  {{ qc_table(sample['peak_counts']) }}
{% endif %}

{% if 'raw_peak_summ' in sample %}
  <h3>Raw peak file statistics</h3>
  <table class='qc_table'>
    <tbody>
      {% for field, value in sample['raw_peak_summ'].iteritems() %}
      <tr>
        <td>{{ field }}</td>
        <td>{{ value }}</td>
      </tr>
      {% endfor %}
    </tbody>
  </table>
{% endif %}

{% if 'raw_peak_dist' in sample %}
  {{ inline_img(sample['raw_peak_dist']) }}
{% endif %}

{% if 'naive_peak_summ' in sample %}
  <h3>Naive overlap peak file statistics</h3>

  <table class='qc_table'>
    <tbody>
      {% for field, value in sample['naive_peak_summ'].iteritems() %}
      <tr>
        <td>{{ field }}</td>
        <td>{{ value }}</td>
      </tr>
      {% endfor %}
    </tbody>
  </table>
{% endif %}

{% if 'naive_peak_dist' in sample %}
  {{ inline_img(sample['naive_peak_dist']) }}
{% endif %}

{% if 'idr_peak_summ' in sample %}
  <h3>IDR peak file statistics</h3>

  <table class='qc_table'>
    <tbody>
      {% for field, value in sample['idr_peak_summ'].iteritems() %}
      <tr>
        <td>{{ field }}</td>
        <td>{{ value }}</td>
      </tr>
      {% endfor %}
    </tbody>
  </table>
{% endif %}

{% if 'idr_peak_dist' in sample %}
  {{ inline_img(sample['idr_peak_dist']) }}

<pre>
For a good ATAC-seq experiment in human, you expect to get 100k-200k peaks
for a specific cell type.
</pre>
{% endif %}


{% if 'gc_bias' in sample %}
  <h2>Sequence quality metrics</h2>
  <h3>GC bias</h3>
  {{ inline_img(sample['gc_bias']) }}
<pre>
Open chromatin assays are known to have significant GC bias. Please take this
into consideration as necessary.
</pre>
{% endif %}

{% if 'enrichment_plots' in sample or 'annot_enrichments' in sample %}
  <h2>Annotation-based quality metrics</h2>
{% endif %}

{% if 'enrichment_plots' in sample %}
  <h3>Enrichment plots (TSS)</h3>
  {{ inline_img(sample['enrichment_plots']['tss']) }}
  <pre>
Open chromatin assays should show enrichment in open chromatin sites, such as
TSS's. An average TSS enrichment is above 6-7. A strong TSS enrichment is
above 10.
  </pre>
{% endif %}

{% if 'annot_enrichments' in sample %}
  <h3>Annotated genomic region enrichments</h3>
  <table class='qc_table'>
    <tbody>
      {% for field, value in sample['annot_enrichments'].iteritems() %}
      <tr>
        <td>{{ field }}</td>
        <td>{{ '{0:,}'.format(value[0]) }}</td>
        <td>{{ '{0:.3f}'.format(value[1]) }}</td>
      </tr>
      {% endfor %}
    </tbody>
  </table>
<pre>
Signal to noise can be assessed by considering whether reads are falling into
known open regions (such as DHS regions) or not. A high fraction of reads
should fall into the universal (across cell type) DHS set. A small fraction
should fall into the blacklist regions. A high set (though not all) should
fall into the promoter regions. A high set (though not all) should fall into
the enhancer regions. The promoter regions should not take up all reads, as
it is known that there is a bias for promoters in open chromatin assays.
</pre>
{% endif %}

{% if 'roadmap_plot' in sample %}
  <h2>Comparison to Roadmap DNase</h2>
  {{ inline_img(sample['roadmap_plot']) }}
<pre>
This bar chart shows the correlation between the Roadmap DNase samples to
your sample, when the signal in the universal DNase peak region sets are
compared. The closer the sample is in signal distribution in the regions
to your sample, the higher the correlation.
</pre>
{% endif %}


</body>

</html>
""")


# ===========================================================

def parse_args():
    '''
    Set up the package to be run from the command line
    '''

    parser = argparse.ArgumentParser(description='ATAC-seq QC package')

    # Directories and prefixes
    parser.add_argument('--workdir', help='Working directory')
    parser.add_argument('--outdir', help='Output directory')
    parser.add_argument('--outprefix', help='Output prefix')

    # Annotation files
    parser.add_argument('--genome', help='Genome build used')
    parser.add_argument('--chromsizes', help='chromsizes file')
    parser.add_argument('--ref', help='Reference fasta file')
    parser.add_argument('--tss', help='TSS file')
    parser.add_argument('--dnase', help='Open chromatin region file')
    parser.add_argument('--blacklist', help='Blacklisted region file')
    parser.add_argument('--prom', help='Promoter region file')
    parser.add_argument('--enh', help='Enhancer region file')
    parser.add_argument('--reg2map_bed', help='file of regions used to generate reg2map signals')
    parser.add_argument('--reg2map', help='file with cell type signals')
    parser.add_argument('--meta', help='Roadmap metadata')

    # Choose which mode
    parser.add_argument('--pipeline',
                        default=None,
                        help='Specific pipeline was used')

    # Mode 1: Provide an input prefix
    parser.add_argument('--inprefix',
                        default='test',
                        help='Input file prefix')

    # Mode 2: Define every possible QC file
    parser.add_argument('--fastq1',
                        help='First set of reads if paired end, \
                              or the single end reads')
    parser.add_argument('--fastq2',
                        help='Second set of reads if paired end')
    parser.add_argument('--alignedbam', help='BAM file from the aligner')
    parser.add_argument('--alignmentlog', help='Alignment log')
    parser.add_argument('--coordsortbam', help='BAM file sorted by coordinate')
    parser.add_argument('--duplog', help='Picard duplicate metrics file')
    parser.add_argument('--pbc', help='ENCODE library complexity metrics file')
    parser.add_argument('--finalbam', help='Final filtered BAM file')
    parser.add_argument('--finalbed',
                        help='Final filtered alignments in BED format')
    parser.add_argument('--bigwig',
                        help='Final bigwig')
    parser.add_argument('--peaks',
                        help='Peak file')
    parser.add_argument('--naive_overlap_peaks',
                        default=None, help='Naive overlap peak file')
    parser.add_argument('--idr_peaks',
                        default=None, help='IDR peak file')
    parser.add_argument('--use_sambamba_markdup', action='store_true',
                        help='Use sambamba markdup instead of Picard')

    args = parser.parse_args()

    # Set up all variables
    INPUT_PREFIX = os.path.join(args.workdir, args.inprefix)
    OUTPUT_PREFIX = os.path.join(args.outdir, args.outprefix)
    os.system('mkdir -p {0}'.format(args.outdir))
    NAME = args.outprefix

    # Set up annotations
    GENOME = args.genome
    CHROMSIZES = args.chromsizes
    REF = args.ref
    TSS = args.tss
    DNASE = args.dnase
    BLACKLIST = args.blacklist
    PROM = args.prom
    ENH = args.enh
    REG2MAP_BED = args.reg2map_bed
    REG2MAP = args.reg2map
    ROADMAP_META = args.meta

    # If mode 1 - TO BE DEPRECATED. In this case, the module is run with
    # Jin's pipeline
    if args.pipeline == 'kundajelab':
        FASTQ = args.fastq1
        ALIGNED_BAM = '{0}.bam'.format(INPUT_PREFIX)
        ALIGNMENT_LOG = '{0}.align.log'.format(INPUT_PREFIX)
        COORDSORT_BAM = '{0}.nodup.bam'.format(INPUT_PREFIX)
        DUP_LOG = '{0}.dup.qc'.format(INPUT_PREFIX)
        PBC_LOG = '{0}.nodup.pbc.qc'.format(INPUT_PREFIX)
        FINAL_BAM = '{0}.nodup.nonchrM.bam'.format(INPUT_PREFIX)
        FINAL_BED = '{0}.nodup.nonchrM.tn5.bed.gz'.format(INPUT_PREFIX)
        BIGWIG = '{0}.nodup.nonchrM.tn5.pf.pval.signal.bigwig'.format(
            INPUT_PREFIX)
        PEAKS = '{0}.nodup.nonchrM.tn5.pf_peaks.narrowPeak'.format(
            INPUT_PREFIX)
    else:  # mode 2
        FASTQ = args.fastq1
        ALIGNED_BAM = args.alignedbam
        ALIGNMENT_LOG = args.alignmentlog
        COORDSORT_BAM = args.coordsortbam
        DUP_LOG = args.duplog
        PBC_LOG = args.pbc
        FINAL_BAM = args.finalbam
        FINAL_BED = args.finalbed
        BIGWIG = args.bigwig
        PEAKS = args.peaks
        NAIVE_OVERLAP_PEAKS = args.naive_overlap_peaks
        IDR_PEAKS = args.idr_peaks
        USE_SAMBAMBA_MARKDUP = args.use_sambamba_markdup

    return NAME, OUTPUT_PREFIX, REF, TSS, DNASE, BLACKLIST, PROM, ENH, \
        REG2MAP_BED, REG2MAP, ROADMAP_META, GENOME, CHROMSIZES, FASTQ, ALIGNED_BAM, \
        ALIGNMENT_LOG, COORDSORT_BAM, DUP_LOG, PBC_LOG, FINAL_BAM, \
        FINAL_BED, BIGWIG, PEAKS, NAIVE_OVERLAP_PEAKS, IDR_PEAKS, \
        USE_SAMBAMBA_MARKDUP


def main():

    # Parse args
    [NAME, OUTPUT_PREFIX, REF, TSS, DNASE, BLACKLIST, PROM, ENH, REG2MAP_BED, REG2MAP,
     ROADMAP_META, GENOME, CHROMSIZES, FASTQ, ALIGNED_BAM, ALIGNMENT_LOG, COORDSORT_BAM,
     DUP_LOG, PBC_LOG, FINAL_BAM, FINAL_BED, BIGWIG, PEAKS,
     NAIVE_OVERLAP_PEAKS, IDR_PEAKS, USE_SAMBAMBA_MARKDUP] = parse_args()
    MITO_CHR_NAME = 'chrM'

    # Set up the log file and timing
    logging.basicConfig(filename='test.log', level=logging.DEBUG)
    start = timeit.default_timer()

    # First check if paired/single
    paired_status = determine_paired(FINAL_BAM)

    # Also get read length
    read_len = get_read_length(FASTQ)

    # Sequencing metrics: Bowtie1/2 alignment log, chrM, GC bias
    BOWTIE_STATS = get_bowtie_stats(ALIGNMENT_LOG)
    chr_m_reads, fraction_chr_m = get_chr_m(COORDSORT_BAM, MITO_CHR_NAME)
    gc_out, gc_plot, gc_summary = get_gc(FINAL_BAM,
                                         REF,
                                         OUTPUT_PREFIX)

    # Library complexity: Preseq results, NRF, PBC1, PBC2
    picard_est_library_size = get_picard_complexity_metrics(ALIGNED_BAM,
                                                            OUTPUT_PREFIX)
    preseq_data, preseq_log = run_preseq(ALIGNED_BAM, OUTPUT_PREFIX) # SORTED BAM
    #preseq_data = '/srv/scratch/dskim89/ataqc/results/2016-03-27.ENCODE_Hardison_Wold/wold/forebrain_e14-5.b1/forebrain_e14-5.b1.preseq.dat'
    #preseq_log = '/srv/scratch/dskim89/ataqc/results/2016-03-27.ENCODE_Hardison_Wold/wold/forebrain_e14-5.b1/forebrain_e14-5.b1.preseq.log'

    encode_lib_metrics = get_encode_complexity_measures(PBC_LOG)

    # Filtering metrics: duplicates, map quality
    num_mapq, fract_mapq = get_fract_mapq(ALIGNED_BAM)

    if USE_SAMBAMBA_MARKDUP:
        read_dups, percent_dup = get_sambamba_dup_stats(DUP_LOG, paired_status)
    else:
        read_dups, percent_dup = get_picard_dup_stats(DUP_LOG, paired_status)

    mito_dups, fract_dups_from_mito = get_mito_dups(ALIGNED_BAM,
                                                    OUTPUT_PREFIX,
                                                    MITO_CHR_NAME,
                                                    paired_status,
                                                    use_sambamba=USE_SAMBAMBA_MARKDUP)
    [flagstat, mapped_count] = get_samtools_flagstat(ALIGNED_BAM)

    # Final read statistics
    first_read_count, final_read_count, \
        fract_reads_left = get_final_read_count(ALIGNED_BAM,
                                                FINAL_BAM)

    # Insert size distribution - CAN'T GET THIS FOR SE FILES
    if paired_status == "Paired-ended":
        insert_data, insert_plot = get_insert_distribution(FINAL_BAM,
                                                           OUTPUT_PREFIX)
    else:
        insert_data = ''
        insert_plot = ''

    # Enrichments: V plot for enrichment
    tss_plot_file, tss_plot_large_file, tss_point_val = make_tss_plot(FINAL_BAM, # Use final to avoid duplicates
                                                                      TSS,
                                                                      OUTPUT_PREFIX,
                                                                      CHROMSIZES,
                                                                      read_len)

    # Signal to noise: reads in DHS regions vs not, reads falling
    # into blacklist regions
    reads_dnase, fract_dnase, reads_blacklist, fract_blacklist, \
        reads_prom, fract_prom, reads_enh, fract_enh, \
        reads_peaks, fract_peaks = get_signal_to_noise(FINAL_BED,
                                                       DNASE,
                                                       BLACKLIST,
                                                       PROM,
                                                       ENH,
                                                       PEAKS)

    # Also need to run n-nucleosome estimation
    if paired_status == 'Paired-ended':
        nucleosomal_qc = fragment_length_qc(read_picard_histogram(insert_data))
    else:
        nucleosomal_qc = ''

    # Peak metrics
    peak_counts = get_peak_counts(PEAKS, NAIVE_OVERLAP_PEAKS, IDR_PEAKS)
    raw_peak_summ, raw_peak_dist = get_region_size_metrics(PEAKS)
    naive_peak_summ, naive_peak_dist = get_region_size_metrics(NAIVE_OVERLAP_PEAKS)
    idr_peak_summ, idr_peak_dist = get_region_size_metrics(IDR_PEAKS)

    # Compare to roadmap
    roadmap_compare_plot = compare_to_roadmap(BIGWIG, REG2MAP_BED, REG2MAP,
                                              ROADMAP_META, OUTPUT_PREFIX)

    # Finally output the bar chart of reads
    read_count_data = [first_read_count, first_read_count*fract_mapq,
                       first_read_count*fract_mapq*(1-float(percent_dup)),
                       final_read_count]
    read_count_labels = ['Start', 'q>30', 'dups removed',
                         'chrM removed (final)']
    read_tracker_plot = track_reads(read_count_data, read_count_labels)

    # Take all this info and render the html file
    SAMPLE_INFO = OrderedDict([
        ('Sample', NAME),
        ('Genome', GENOME),
        ('Paired/Single-ended', paired_status),
        ('Read length', read_len),
    ])

    SUMMARY_STATS = OrderedDict([
        ('Read count from sequencer', first_read_count),
        ('Read count successfully aligned', mapped_count),
        ('Read count after filtering for mapping quality', num_mapq),
        ('Read count after removing duplicate reads',
            int(num_mapq - read_dups)),
        ('Read count after removing mitochondrial reads (final read count)',
            final_read_count),
    ])

    FILTERING_STATS = OrderedDict([
        ('Mapping quality > q30 (out of total)', (num_mapq, fract_mapq)),
        ('Duplicates (after filtering)', (read_dups, percent_dup)),
        ('Mitochondrial reads (out of total)', (chr_m_reads, fraction_chr_m)),
        ('Duplicates that are mitochondrial (out of all dups)',
            (mito_dups, fract_dups_from_mito)),
        ('Final reads (after all filters)', (final_read_count,
                                             fract_reads_left)),
    ])

    ENRICHMENT_PLOTS = {
        'tss': b64encode(open(tss_plot_large_file, 'rb').read())
    }

    ANNOT_ENRICHMENTS = OrderedDict([
        ('Fraction of reads in universal DHS regions', (reads_dnase,
                                                        fract_dnase)),
        ('Fraction of reads in blacklist regions', (reads_blacklist,
                                                    fract_blacklist)),
        ('Fraction of reads in promoter regions', (reads_prom, fract_prom)),
        ('Fraction of reads in enhancer regions', (reads_enh, fract_enh)),
        ('Fraction of reads in called peak regions', (reads_peaks,
                                                      fract_peaks)),
    ])

    SAMPLE = OrderedDict([
        ('Name', NAME),
        ('basic_info', SAMPLE_INFO),

        # Summary
        ('summary_stats', SUMMARY_STATS),
        ('read_tracker', read_tracker_plot),

        # Alignment statistics
        ('bowtie_stats', BOWTIE_STATS),
        ('samtools_flagstat', flagstat),

        # Filtering statistics
        ('filtering_stats', FILTERING_STATS),

        # Library complexity statistics
        ('encode_lib_complexity', encode_lib_metrics),
        ('picard_est_library_size', picard_est_library_size),
        ('yield_prediction', preseq_plot(preseq_data)),

        # Fragment length statistics
        ('fraglen_dist', fragment_length_plot(insert_data)),
        ('nucleosomal', nucleosomal_qc),

        # Peak metrics
        ('peak_counts', peak_counts),
        ('raw_peak_summ', raw_peak_summ),
        ('naive_peak_summ', naive_peak_summ),
        ('idr_peak_summ', idr_peak_summ),
        ('raw_peak_dist', raw_peak_dist),
        ('naive_peak_dist', naive_peak_dist),
        ('idr_peak_dist', idr_peak_dist),

        # GC
        ('gc_bias', plot_gc(gc_out)),

        # Annotation based statistics
        ('enrichment_plots', ENRICHMENT_PLOTS),
        ('TSS_enrichment', tss_point_val),
        ('annot_enrichments', ANNOT_ENRICHMENTS),

        # Roadmap plot
        ('roadmap_plot', roadmap_compare_plot),
    ])

    results = open('{0}_qc.html'.format(OUTPUT_PREFIX), 'w')
    results.write(html_template.render(sample=SAMPLE))
    results.close()

    # Also produce a text file of relevant stats (so that another module
    # can combine the stats) and put in using ordered dictionary
    textfile = open('{0}_qc.txt'.format(OUTPUT_PREFIX), 'w')
    for key, value in SAMPLE.iteritems():
        # Make sure to not get b64encode
        if isinstance(value, str) and (len(value) < 300) and (len(value) > 0):
            textfile.write('{0}\t{1}\n'.format(key, value))
        elif isinstance(value, int):
            textfile.write('{0}\t{1}\n'.format(key, value))
        elif isinstance(value, float):
            textfile.write('{0}\t{1}\n'.format(key, value))
        elif isinstance(value, OrderedDict):
            for dict_key, dict_value in value.iteritems():
                if isinstance(dict_value, tuple):
                    textfile.write('{0}'.format(dict_key))
                    for tuple_val in dict_value:
                        textfile.write('\t{0}'.format(tuple_val))
                    textfile.write('\n')
                else:
                    textfile.write('{0}\t{1}\n'.format(dict_key, dict_value))
        # QC tables go here
        elif isinstance(value, list):
            if 'bowtie' in value[0]: # Hack, fix this
                continue
            for result in value:
                textfile.write('{0}\t{1}\n'.format(result.metric,
                                                   result.message))
        else:
            pass
    textfile.close()

    stop = timeit.default_timer()
    print("Run time:", str(datetime.timedelta(seconds=int(stop - start))))

    return None

if __name__=='__main__':
    main()
