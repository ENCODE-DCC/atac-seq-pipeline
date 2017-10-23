#!/usr/bin/env python

# ENCODE DCC bowtie2 aligner python script
# Author: Jin Lee (leepc12@gmail.com)

import sys
import os
import re
import json
import logging
import gzip
import collections
import argparse
import encode_dcc_common

logger = logging.getLogger(__name__)

def parse_arguments():
    parser = argparse.ArgumentParser(prog='ENCODE DCC bowtie2 aligner python script',
                                        description='')
    parser.add_argument('index_prefix_or_tar', type=str,
                        help='Path for prefix (or a tarball) for reference bowtie2 indices.')
    parser.add_argument('fastqs', nargs='+', type=str,
                        help='TSV file path or \
                            list of FASTQ files delimited by space. Use \
                            TSV for multiple techincal replicates).')
    parser.add_argument('--nth', type=int, default=1,
                        help='No. threads to parallelize bowtie2.')
    parser.add_argument('--paired-end', action="store_true",
                        help='FASTQs have paired-end')
    group_out = parser.add_mutually_exclusive_group()
    group_out.add_argument('--out-dir', default='.', type=str,
                            help='Output directory. Prefix will be taken from FASTQs.')
    group_out.add_argument('--out-prefix', type=str,
                            help='Output prefix path.')
    group_bowtie2 = parser.add_argument_group(title='bowtie2',
                        description='Bowtie2 settings.')
    group_bowtie2.add_argument('--multimapping', default=0, type=int,
                        help='Multimapping for bowtie2 -k.')
    group_bowtie2.add_argument('--score-min', default='', type=str,
                        help='--score-min for bowtie2 -k.')
    args = parser.parse_args()
    # parse fastqs command line
    if args.fastqs[0].endswith('.gz'): # it's fastq
        args.fastqs = [args.fastqs]
    else: # it's TSV
        args.fastqs = read_tsv(args.fastqs[0])
    # parse --adapters command line
    if args.adapters:
        if os.path.exists(args.adapters[0]): # it's TSV
            args.adapters = read_tsv(args.adapters[0])
        else:
            args.adapters = [args.adapters]
        if args.adapters and len(args.adapters) != len(args.fastqs):
            raise ValueError(
                'fastqs and adapters must have the same dimension.')
    return args

## functions returning filename

def merge_fastqs(fastqs, out_dir):
    prefix1 = os.path.join(out_dir,
        os.path.basename(strip_ext_fastq(fastq1)))
    merged = '{}.merged.fastq.gz'.format(prefix)
    cmd = 'zcat {} | gzip -nc > {}'.format(
        ' '.join(fastqs),
        merged)
    print(cmd)
    print(subprocess.check_output(cmd, shell=True))
    return merged

def bowtie2_se(fastq, ref_prefix, 
        multimapping, score_min, nth, out_dir):
    prefix = os.path.join(out_dir,
        os.path.basename(strip_ext_fastq(fastq)))
    bam = '{}.bam'.format(prefix)
    bai = '{}.bam.bai'.format(prefix)
    align_log = '{}.align.log'.format(prefix)
    cmd1 = 'bowtie2 {} {} --local --threads {} -x {} -U {} 2> {} '
    cmd1 += '| samtools view -Su /dev/stdin | samtools sort - {}'
    cmd1 = cmd1.format(
        '-k {}'.format(multimapping) if multimapping else '',
        '--score-min {}'.format(score_min) if score_min else '',
        nth,
        ref_prefix,
        fastq,
        align_log,
        prefix)
    print(cmd1)
    print(subprocess.check_output(shlex.split(cmd1)))
    cmd2 = 'samtools index {}'.format(bam)
    print(cmd2)
    print(subprocess.check_output(shlex.split(cmd2)))
    return bam, bai, align_log

def bowtie2_pe(fastq1, fastq2, ref_prefix, 
        multimapping, score_min, nth, out_dir):
    prefix = os.path.join(out_dir,
        os.path.basename(strip_ext_fastq(fastq1)))
    print(fastq1,fastq2,prefix)
    bam = '{}.bam'.format(prefix)
    bai = '{}.bam.bai'.format(prefix)
    align_log = '{}.align.log'.format(prefix)
    cmd1 = 'bowtie2 {} {} -X2000 --mm --local --threads {} -x {} '
    cmd1 += '-1 {} -2 {} 2>{} | '
    cmd1 += 'samtools view -Su /dev/stdin | samtools sort - {}'
    cmd1 = cmd1.format(
        '-k {}'.format(multimapping) if multimapping else '',
        '--score-min {}'.format(score_min) if score_min else '',
        nth,
        ref_prefix,
        fastq1,
        fastq2,
        align_log,
        prefix)
    print(cmd1)
    print(subprocess.check_output(shlex.split(cmd1)))
    cmd2 = 'samtools index {}'.format(bam)
    print(cmd2)
    print(subprocess.check_output(shlex.split(cmd2)))
    return bam, bai, align_log

def samtools_flagstat(bam, out_dir):
    prefix = os.path.join(out_dir,
        os.path.basename(strip_ext_fastq(bam)))
    flagstat_qc = '{}.flagstat.qc'.format(prefix)
    cmd = 'samtools flagstat {} > {}'.format(
        bam,
        flagstat_qc)
    print(cmd)
    print(subprocess.check_output(cmd, shell=True))
    return flagstat_qc

def main():
    args = parse_arguments()
    temp_files = []
    R1_to_merge = []
    R2_to_merge = []

    # detect/trim adapters
    for i in range(len(args.fastqs)):
        fastqs = args.fastqs[i] # R1 and R2
        adapters = args.adapters[i] if args.adapters else [None,None]
        if args.paired_end:
            if not args.auto_detect_adapter and \
                not (adapters[0] and adapters[1]):
                trimmed_fastqs = fastqs
            else:
                if args.auto_detect_adapter and \
                    not (adapters[0] and adapters[1]): # detect adapters
                    adapters[0] = detect_adapter(fastqs[0])
                    adapters[1] = detect_adapter(fastqs[1])
                if adapters[0] and adapters[1]:
                    trimmed_fastqs = trim_adapter_pe(
                        fastqs[0], fastqs[1], 
                        adapters[1], adapter[2],
                        args.min_trim_len,
                        args.err_rate,
                        args.out_dir)
                    temp_files.append(trimmed_fastqs)
                else:
                    trimmed_fastqs = fastqs
            R1_to_merge.append(trimmed_fastqs[0])
            R2_to_merge.append(trimmed_fastqs[1])
        else:
            if not args.auto_detect_adapter and \
                not adapters[0]:
                trimmed_fastq = fastq[0]
            else:
                if args.auto_detect_adapter and not adapters[0]:
                    adapters[0] = detect_adapter(fastqs[0])
                if adapters[0]:
                    trimmed_fastqs = trim_adapter_se(
                        fastqs[0],
                        adapters[0],
                        args.min_trim_len,
                        args.err_rate,
                        args.out_dir)
                    temp_files.append(trimmed_fastqs)
                else:
                    trimmed_fastq = fastqs[0]
            R1_to_merge.append(trimmed_fastq)
    # merge fastqs    
    if len(R1_to_merge)>1:
        R1_merged = merge_fastqs(R1_to_merge)
        temp_files.append(R1_merged)
    else:
        R1_merged = R1_to_merge[0]
    if len(R2_to_merge)>1:
        R2_merged = merge_fastqs(R2_to_merge)
        temp_files.append(R2_merged)
    else:
        R2_merged = R2_to_merge[0]
    # generate read length file
    if R1_merged:
        R1_read_length_file = get_read_length_file(R1_merged, args.out_dir)
    if R2_merged:
        R2_read_length_file = get_read_length_file(R2_merged, args.out_dir)
    # bowtie2
    if args.paired_end:
        bam, bai, align_log, flagstat_qc = bowtie2_pe(
            R1_merged, R2_merged, args.ref_prefix,
            args.multimapping, args.score_min, args.nth,
            args.out_dir)
    else:
        bam, bai, align_log, flagstat_qc = bowtie2_se(
            R1_merged, args.ref_prefix,
            args.multimapping, args.score_min, args.nth,
            args.out_dir)
    # remove temporary/intermediate files
    subprocess.check_output(shlex.split(
        'rm -rf {}'.format(' '.join(temp_files))))

if __name__=='__main__':
    main()