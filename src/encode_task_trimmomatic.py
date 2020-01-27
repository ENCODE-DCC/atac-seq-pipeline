#!/usr/bin/env python

# ENCODE DCC Trimmomatic wrapper
# Author: Jin Lee (leepc12@gmail.com)

import sys
import os
import argparse
from encode_lib_common import (
    log, ls_l, mkdir_p, rm_f,
    run_shell_cmd, strip_ext_fastq)
from encode_lib_genomic import (
    locate_trimmomatic)


def parse_arguments(debug=False):
    parser = argparse.ArgumentParser(
        prog='ENCODE DCC Trimmomatic wrapper.')
    parser.add_argument('--fastq1',
                        help='FASTQ R1 to be trimmed.')
    parser.add_argument('--fastq2',
                        help='FASTQ R2 to be trimmed.')
    parser.add_argument('--paired-end', action="store_true",
                        help='Paired-end FASTQs.')
    parser.add_argument('--crop-length', type=int,
                        help='Number of basepair to crop.')
    parser.add_argument('--out-dir-R1', default='', type=str,
                        help='Output directory for cropped R1 fastq.')
    parser.add_argument('--out-dir-R2', default='', type=str,
                        help='Output directory for cropped R2 fastq.')
    parser.add_argument('--trimmomatic-java-heap',
                        help='Trimmomatic\'s Java max. heap: java -jar Trimmomatic.jar '
                             '-Xmx[MAX_HEAP]')
    parser.add_argument('--nth', type=int, default=1,
                        help='Number of threads to parallelize.')    
    parser.add_argument('--log-level', default='INFO',
                        choices=['NOTSET', 'DEBUG', 'INFO',
                                 'WARNING', 'CRITICAL', 'ERROR',
                                 'CRITICAL'],
                        help='Log level')
    args = parser.parse_args()
    if not args.crop_length:
        raise ValueError('Crop length must be > 0.')

    log.setLevel(args.log_level)
    log.info(sys.argv)
    return args


def trimmomatic_se(fastq1, crop_length, out_dir,
                   nth=1, java_heap=None):
    prefix = os.path.join(out_dir,
                          os.path.basename(strip_ext_fastq(fastq1)))
    cropped = '{}.crop_{}bp.fastq.gz'.format(prefix, crop_length)

    if java_heap is None:
        java_heap_param = '-Xmx6G'
    else:
        java_heap_param = '-Xmx{}'.format(java_heap)

    cmd = 'java -XX:ParallelGCThreads=1 {} -jar {} SE -threads {} '
    cmd += '{} {} MINLEN:{} CROP:{}'
    cmd = cmd.format(
        java_heap_param,
        locate_trimmomatic(),
        nth,
        fastq1,
        cropped,
        crop_length,
        crop_length)
    run_shell_cmd(cmd)

    return cropped


def trimmomatic_pe(fastq1, fastq2, crop_length, out_dir_R1, out_dir_R2,
                   nth=1, java_heap=None):
    prefix_R1 = os.path.join(
        out_dir_R1, os.path.basename(strip_ext_fastq(fastq1)))
    prefix_R2 = os.path.join(
        out_dir_R2, os.path.basename(strip_ext_fastq(fastq2)))
    cropped_R1 = '{}.crop_{}bp.fastq.gz'.format(prefix_R1, crop_length)
    cropped_R2 = '{}.crop_{}bp.fastq.gz'.format(prefix_R2, crop_length)
    tmp_cropped_R1 = '{}.tmp'.format(cropped_R1)
    tmp_cropped_R2 = '{}.tmp'.format(cropped_R2)

    if java_heap is None:
        java_heap_param = '-Xmx6G'
    else:
        java_heap_param = '-Xmx{}'.format(java_heap)

    cmd = 'java -XX:ParallelGCThreads=1 {} -jar {} PE -threads {} '
    cmd += '{} {} {} {} {} {} MINLEN:{} CROP:{}'
    cmd = cmd.format(
        java_heap_param,
        locate_trimmomatic(),
        nth,
        fastq1,
        fastq2,
        cropped_R1,
        tmp_cropped_R1,
        cropped_R2,
        tmp_cropped_R2,
        crop_length,
        crop_length)
    run_shell_cmd(cmd)
    rm_f([tmp_cropped_R1, tmp_cropped_R2])

    return cropped_R1, cropped_R2


def main():
    # read params
    args = parse_arguments()

    log.info('Initializing and making output directory...')
    mkdir_p(args.out_dir_R1)
    if args.paired_end:
        mkdir_p(args.out_dir_R2)

    log.info('Cropping fastqs ({} bp) with Trimmomatic...'.format(args.crop_length))
    if args.paired_end:
        trimmomatic_pe(
            args.fastq1, args.fastq2,
            args.crop_length,
            args.out_dir_R1, args.out_dir_R2,
            args.nth,
            args.trimmomatic_java_heap)
    else:
        trimmomatic_se(
            args.fastq1,
            args.crop_length,
            args.out_dir_R1,
            args.nth,
            args.trimmomatic_java_heap)

    log.info('List all files in output directory...')
    ls_l(args.out_dir_R1)
    if args.paired_end:
        ls_l(args.out_dir_R2)

    log.info('All done.')


if __name__ == '__main__':
    main()
