#!/usr/bin/env python

# ENCODE DCC filter python script
# Author: Jin Lee (leepc12@gmail.com)

import sys
import os
import re
import json
import collections
import argparse
import multiprocessing
import subprocess
import copy
import encode_dcc_common

def parse_arguments():
    parser = argparse.ArgumentParser(prog='ENCODE DCC filter python script',
                                        description='')
    parser.add_argument('bam', type=str,
                        help='Path for BAM file.')
    parser.add_argument('--nth', type=int, default=1,
                        help='No. threads to parallelize bowtie2.')
    parser.add_argument('--paired-end', action="store_true",
                        help='FASTQs have paired-end')
    parser.add_argument('--out-dir', default='.', type=str,
                            help='Output directory. Prefix will be taken from FASTQs.')
    parser.add_argument('--multimapping', default=0, type=int,
                        help='Multimapping for bowtie2 -k.')
    parser.add_argument('--mapq-thresh', default=30, type=int,
                        help='')
    parser.add_argument('--no-dup-removal', action="store_true",
                        help='')
    args = parser.parse_args()
    return args



def bowtie2_se(fastq, ref_prefix, 
        multimapping, score_min, nth, out_dir):
    prefix = os.path.join(out_dir,
        os.path.basename(strip_ext_fastq(fastq)))
    bam = '{}.bam'.format(prefix)
    align_log = '{}.align.log'.format(prefix)
    cmd = 'bowtie2 {} {} --local --threads {} -x {} -U {} 2> {} '
    cmd += '| samtools view -Su /dev/stdin | samtools sort - {}'
    cmd = cmd.format(
        '-k {}'.format(multimapping) if multimapping else '',
        '--score-min {}'.format(score_min) if score_min else '',
        nth,
        ref_prefix,
        fastq,
        align_log,
        prefix)
    print(cmd)
    print(subprocess.check_output(shlex.split(cmd)))
    return bam, align_log

def bowtie2_pe(fastq1, fastq2, ref_prefix, 
        multimapping, score_min, nth, out_dir):
    prefix = os.path.join(out_dir,
        os.path.basename(strip_ext_fastq(fastq1)))
    print(fastq1,fastq2,prefix)
    bam = '{}.bam'.format(prefix)
    bai = '{}.bam.bai'.format(prefix)
    align_log = '{}.align.log'.format(prefix)
    cmd = 'bowtie2 {} {} -X2000 --mm --local --threads {} -x {} '
    cmd += '-1 {} -2 {} 2>{} | '
    cmd += 'samtools view -Su /dev/stdin | samtools sort - {}'
    cmd = cmd.format(
        '-k {}'.format(multimapping) if multimapping else '',
        '--score-min {}'.format(score_min) if score_min else '',
        nth,
        ref_prefix,
        fastq1,
        fastq2,
        align_log,
        prefix)
    print(cmd)
    print(subprocess.check_output(shlex.split(cmd)))
    return bam, align_log

def samtools_index(bam, out_dir):
    prefix = os.path.join(out_dir,
        os.path.basename(strip_ext_bam(bam)))
    bai = '{}.bam.bai'.format(prefix)
    cmd = 'samtools index {}'.format(bam)
    print(cmd)
    print(subprocess.check_output(shlex.split(cmd)))
    return bai    

def samtools_flagstat(bam, out_dir):
    prefix = os.path.join(out_dir,
        os.path.basename(strip_ext_bam(bam)))
    flagstat_qc = '{}.flagstat.qc'.format(prefix)
    cmd = 'samtools flagstat {} > {}'.format(
        bam,
        flagstat_qc)
    print(cmd)
    print(subprocess.check_output(cmd, shell=True))
    return flagstat_qc

def main():
    # read params
    args = parse_arguments()
    if args.debug:
        print(sys.argv)
    # for multithreading and keyboard interrupt handling (SIGINT)
    sigint_handler = signal.signal(signal.SIGINT, signal.SIG_IGN)
    pool = multiprocessing.Pool(args.nth)
    signal.signal(signal.SIGINT, sigint_handler)

    try:
        # make out_dir
        mkdir_p(args.out_dir)
        # generate read length file
        ret_val1 = pool.apply_async(get_read_length_file,
                                (args.fastqs[0], args.out_dir))
        ret_val2 = pool.apply_async(get_read_length_file,
                                (args.fastqs[1], args.out_dir))    
        R1_read_length_file = ret_val1.get()
        R2_read_length_file = ret_val2.get()
        # bowtie2
        if args.paired_end:
            ret_val3 = pool.apply_async(bowtie2_pe,(
                args.fastqs[0], args.fastqs[1], args.ref_prefix,
                args.multimapping, args.score_min, args.nth,
                args.out_dir))
        else:
            ret_val3 = pool.apply_async(bowtie2_se,(
                args.fastqs[0], args.ref_prefix,
                args.multimapping, args.score_min, args.nth,
                args.out_dir))
        # get()
        bam, align_log = ret_val3.get()
        # samtools index
        ret_val1 = pool.apply_async(samtools_index,(bam, args.out_dir))
        # samtools flagstat qc
        ret_val2 = pool.apply_async(samtools_flagstat,(bam, args.out_dir))
        # get()
        bai = ret_val1.get()
        flagstat_qc = ret_val2.get()
    except KeyboardInterrupt:
        print("Keyboard interrupted, terminating workers...")
        pool.terminate()
    except:
        print("Unknown exception, terminating workers...")
        pool.terminate()
    else:
        if args.debug:
            print("Normal termination.")
        pool.close()

    pool.join()

if __name__=='__main__':
    main()