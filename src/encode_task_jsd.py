#!/usr/bin/env python

# ENCODE DCC fingerprint/JSD plot wrapper
# Author: Jin Lee (leepc12@gmail.com)

import sys
import os
import argparse
from encode_lib_genomic import *

def parse_arguments():
    parser = argparse.ArgumentParser(prog='ENCODE DCC Fingerprint/JSD plot.',
                                        description='')
    parser.add_argument('bams', nargs='+', type=str,
                        help='List of paths for filtered experiment BAM files.')
    parser.add_argument('--ctl-bam', type=str, default='',
                        help='Path for filtered control BAM file.')
    parser.add_argument('--blacklist', type=str, default='',
                        help='Blacklist BED file.')
    parser.add_argument('--mapq-thresh', default=30, type=int,
                        help='Threshold for low MAPQ reads removal.')
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

def fingerprint(bams, ctl_bam, blacklist, mapq_thresh, nth, out_dir):
    # make bam index (.bai) first
    for bam in bams:
        samtools_index(bam, nth)
    if ctl_bam:
        samtools_index(ctl_bam, nth)

    prefix = os.path.join(out_dir,
        os.path.basename(strip_ext_bam(bams[0])))
    plot_png =  '{}.jsd_plot.png'.format(prefix)
    tmp_log = '{}.jsd.tmp'.format(prefix)

    labels = []
    bam_paths = []
    jsd_qcs = []
    for i, bam in enumerate(bams):
        prefix_ = os.path.join(out_dir,
            os.path.basename(strip_ext_bam(bam)))
        jsd_qcs.append('rep{}.{}.jsd.qc'.format(i+1,prefix_))
        labels.append('rep{}'.format(i+1)) # repN
        bam_paths.append(bam)
    # add control
    if ctl_bam:
        labels.append('ctl1')
        bam_paths.append(ctl_bam)

    cmd = 'LC_ALL=en_US.UTF-8 LANG=en_US.UTF-8 plotFingerprint -b {} '
    if ctl_bam:
        cmd += '--JSDsample {} '.format(ctl_bam)
    cmd += '--labels {} '
    cmd += '--outQualityMetrics {} '
    cmd += '--minMappingQuality {} '
    cmd += '-T "Fingerprints of different samples" '
    if blacklist:
        cmd += '--blackListFileName {} '.format(blacklist)
    cmd += '--numberOfProcessors {} '
    cmd += '--plotFile {}'
    cmd = cmd.format(
        ' '.join(bam_paths),
        ' '.join(labels),
        tmp_log,
        mapq_thresh,
        nth,
        plot_png)
    run_shell_cmd(cmd)

    # parse tmp_log to get jsd_qc for each exp replicate
    with open(tmp_log,'r') as fp:
        for i, line in enumerate(fp.readlines()): # i is rep_id-1
            if i==0: continue
            if i>len(jsd_qcs): break
            with open(jsd_qcs[i-1],'w') as fp2:
                # removing repN from lines
                fp2.write('\t'.join(line.strip().split('\t')[1:]))
    rm_f(tmp_log)
    return plot_png, jsd_qcs

def main():
    # read params
    args = parse_arguments()

    log.info('Initializing and making output directory...')
    mkdir_p(args.out_dir)

    # declare temp arrays
    temp_files = [] # files to deleted later at the end

    log.info('Plotting Fingerprint on BAMs and calculating JSD...')
    plot_png, jsd_qcs = fingerprint(
        args.bams, args.ctl_bam, args.blacklist, args.mapq_thresh,
        args.nth, args.out_dir)

    log.info('List all files in output directory...')
    ls_l(args.out_dir)

    log.info('All done.')

if __name__=='__main__':
    main()