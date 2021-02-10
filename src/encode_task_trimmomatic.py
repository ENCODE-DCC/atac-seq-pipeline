#!/usr/bin/env python

# ENCODE DCC Trimmomatic wrapper
# Author: Jin Lee (leepc12@gmail.com)

import sys
import os
import argparse
from encode_lib_common import (
    assert_file_not_empty, log, ls_l, mkdir_p, rm_f,
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
    parser.add_argument('--crop-length', type=int, required=True,
                        help='Number of basepair to crop.'
                             'Trimmomatic\'s parameter CROP.')
    parser.add_argument('--crop-length-tol', type=int, default=2,
                        help='Crop length tolerance to keep shorter reads '
                             'around the crop length. '
                             'Trimmomatic\'s parameter MINLEN will be --crop-length '
                             '- abs(--crop-length-tol).')
    parser.add_argument('--phred-score-format',
                        default='auto',
                        choices=['auto', 'phred33', 'phred64'],
                        help='Base encoding for Phred scores in FASTQs. '
                             'If it is not auto then -phred33 or -phred64 to '
                             'Trimmomatic\'s command line.')
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


def trimmomatic_se(fastq1, crop_length, crop_length_tol,
                   phred_score_format, out_dir,
                   nth=1, java_heap=None):
    prefix = os.path.join(out_dir,
                          os.path.basename(strip_ext_fastq(fastq1)))
    crop_length_tol = abs(crop_length_tol)
    min_length = crop_length - crop_length_tol
    cropped = '{p}.crop_{cl}-{tol}bp.fastq.gz'.format(
        p=prefix, cl=crop_length, tol=crop_length_tol)

    if java_heap is None:
        java_heap_param = '-Xmx6G'
    else:
        java_heap_param = '-Xmx{}'.format(java_heap)

    phred_score_format = phred_score_format.lower()
    if phred_score_format == 'auto':
        phred_score_param = ''
    elif phred_score_format == 'phred33':
        phred_score_param = '-phred33'
    elif phred_score_format == 'phred64':
        phred_score_param = '-phred64'
    else:
        raise ValueError('Wrong phred_score_format!')

    cmd = 'java -XX:ParallelGCThreads=1 {param} -jar {jar} SE -threads {nth} {phred_score_param} ' \
          '{fq1} {cropped} MINLEN:{ml} CROP:{cl}'.format(
        param=java_heap_param,
        jar=locate_trimmomatic(),
        nth=nth,
        phred_score_param=phred_score_param,
        fq1=fastq1,
        cropped=cropped,
        ml=min_length,
        cl=crop_length,
    )
    run_shell_cmd(cmd)

    return cropped


def trimmomatic_pe(fastq1, fastq2, crop_length, crop_length_tol,
                   phred_score_format, out_dir_R1, out_dir_R2,
                   nth=1, java_heap=None):
    prefix_R1 = os.path.join(
        out_dir_R1, os.path.basename(strip_ext_fastq(fastq1)))
    prefix_R2 = os.path.join(
        out_dir_R2, os.path.basename(strip_ext_fastq(fastq2)))

    crop_length_tol = abs(crop_length_tol)
    min_length = crop_length - crop_length_tol

    cropped_R1 = '{p}.crop_{cl}-{tol}bp.fastq.gz'.format(
        p=prefix_R1, cl=crop_length, tol=crop_length_tol)
    cropped_R2 = '{p}.crop_{cl}-{tol}bp.fastq.gz'.format(
        p=prefix_R2, cl=crop_length, tol=crop_length_tol)
    tmp_cropped_R1 = '{}.tmp'.format(cropped_R1)
    tmp_cropped_R2 = '{}.tmp'.format(cropped_R2)

    if java_heap is None:
        java_heap_param = '-Xmx6G'
    else:
        java_heap_param = '-Xmx{}'.format(java_heap)

    phred_score_format = phred_score_format.lower()
    if phred_score_format == 'auto':
        phred_score_param = ''
    elif phred_score_format == 'phred33':
        phred_score_param = '-phred33'
    elif phred_score_format == 'phred64':
        phred_score_param = '-phred64'
    else:
        raise ValueError('Wrong phred_score_format!')

    cmd = 'java -XX:ParallelGCThreads=1 {param} -jar {jar} PE -threads {nth} {phred_score_param} ' \
          '{fq1} {fq2} {cropped1} {tmp_cropped1} {cropped2} {tmp_cropped2} ' \
          'MINLEN:{ml} CROP:{cl}'.format(
        param=java_heap_param,
        jar=locate_trimmomatic(),
        nth=nth,
        phred_score_param=phred_score_param,
        fq1=fastq1,
        fq2=fastq2,
        cropped1=cropped_R1,
        tmp_cropped1=tmp_cropped_R1,
        cropped2=cropped_R2,
        tmp_cropped2=tmp_cropped_R2,
        ml=min_length,
        cl=crop_length,
    )
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

    log.info(
        'Cropping fastqs with Trimmomatic... '
        'crop_length={cl}, crop_length_tol={clt}'.format(
            cl=args.crop_length,
            clt=args.crop_length_tol))
    if args.paired_end:
        cropped_R1, cropped_R2 = trimmomatic_pe(
            args.fastq1, args.fastq2,
            args.crop_length, args.crop_length_tol,
            args.phred_score_format,
            args.out_dir_R1, args.out_dir_R2,
            args.nth,
            args.trimmomatic_java_heap)
    else:
        cropped_R1 = trimmomatic_se(
            args.fastq1,
            args.crop_length, args.crop_length_tol,
            args.phred_score_format,
            args.out_dir_R1,
            args.nth,
            args.trimmomatic_java_heap)

    log.info('List all files in output directory...')
    ls_l(args.out_dir_R1)
    if args.paired_end:
        ls_l(args.out_dir_R2)

    log.info('Checking if output is empty...')
    assert_file_not_empty(cropped_R1, help=
        'No reads in FASTQ after cropping. crop_length might be too high? '
        'While cropping, Trimmomatic (with MINLEN=crop_length-abs(crop_length_tol)) '
        'excludes all reads SHORTER than crop_length.')

    log.info('All done.')


if __name__ == '__main__':
    main()
