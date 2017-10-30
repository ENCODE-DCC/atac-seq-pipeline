#!/usr/bin/env python

# ENCODE DCC bowtie2 aligner python script
# Author: Jin Lee (leepc12@gmail.com), Daniel Kim

import sys
import os
import re
import argparse
import multiprocessing
from encode_dcc_common import *

def parse_arguments():
    parser = argparse.ArgumentParser(prog='ENCODE DCC bowtie2 aligner python script',
                                        description='')
    parser.add_argument('bowtie2_index_prefix_or_tar', type=str,
                        help='Path for prefix (or a tarball .tar) \
                            for reference bowtie2 index. \
                            Prefix must be [PREFIX].1.bt2. \
                            Tar ball must be packed without compression \
                            by tar cvf [TAR] [TAR_PREFIX].*.bt2.')
    parser.add_argument('fastqs', nargs='+', type=str,
                        help='List of FASTQs delimited by space. \
                            FASTQs must be compressed with gzip (with .gz).')
    parser.add_argument('--nth', type=int, default=1,
                        help='Number of threads to parallelize.')
    parser.add_argument('--paired-end', action="store_true",
                        help='Paired-end FASTQs.')
    parser.add_argument('--out-dir', default='.', type=str,
                            help='Output directory.')
    parser.add_argument('--multimapping', default=0, type=int,
                        help='Multimapping for bowtie2 -k.')
    parser.add_argument('--score-min', default='', type=str,
                        help='--score-min for bowtie2 -k.')
    parser.add_argument('--log-level', default='INFO', 
                        choices=['NOTSET','DEBUG','INFO',
                            'WARNING','CRITICAL','ERROR','CRITICAL'],
                        help='Log level')
    args = parser.parse_args()

    # check if fastqs have correct dimension
    if args.paired_end and len(args.fastqs)!=2:
        raise ValueError('Need 2 fastqs for paired end.')
    if not args.paired_end and len(args.fastqs)!=1:
        raise ValueError('Need 1 fastq for single end.')

    log.setLevel(args.log_level)
    log.info(sys.argv)    
    return args

def get_read_length(fastq):
    # code extracted from Daniel Kim's ATAQC module
    # https://github.com/kundajelab/ataqc/blob/master/run_ataqc.py
    def getFileHandle(filename, mode="r"):
        if (re.search('.gz$',filename) or re.search('.gzip',filename)):
            if (mode=="r"):
                mode="rb";
            return gzip.open(filename,mode)
        else:
            return open(filename,mode)
    total_reads_to_consider = 1000000
    line_num = 0
    total_reads_considered = 0
    max_length = 0
    with getFileHandle(fastq, 'rb') as fp:
        for line in fp:
            if line_num % 4 == 1:
                if len(line.strip()) > max_length:
                    max_length = len(line.strip())
                total_reads_considered += 1
            if total_reads_considered >= total_reads_to_consider:
                break
            line_num += 1
    return int(max_length)

def make_read_length_file(fastq, out_dir):
    prefix = os.path.join(out_dir,
        os.path.basename(strip_ext_fastq(fastq)))
    txt = '{}.read_length.txt'.format(prefix)
    read_length = get_read_length(fastq)
    with open(txt,'w') as fp:
        fp.write(str(read_length))
    return txt

def bowtie2_se(fastq, ref_index_prefix, 
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
        ref_index_prefix,
        fastq,
        align_log,
        prefix)
    run_shell_cmd(cmd)
    return bam, align_log

def bowtie2_pe(fastq1, fastq2, ref_index_prefix, 
        multimapping, score_min, nth, out_dir):
    prefix = os.path.join(out_dir,
        os.path.basename(strip_ext_fastq(fastq1)))
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
        ref_index_prefix,
        fastq1,
        fastq2,
        align_log,
        prefix)
    run_shell_cmd(cmd)
    return bam, align_log

def chk_bowtie2_index(prefix):    
    index_1 = '{}.1.bt2'.format(prefix)
    if not os.path.exists(index_1):
        raise Exception("Bowtie2 index does not exists. "+
            "Prefix = {}".format(prefix))

def main():
    # read params
    args = parse_arguments()
    log.info('Initializing and making output directory...')

    # make out_dir
    mkdir_p(args.out_dir)

    # declare temp arrays
    temp_files = [] # files to deleted later at the end

    # generate read length file
    log.info('Generating read length file...')
    R1_read_length_file = make_read_length_file(
                            args.fastqs[0], args.out_dir)
    if args.paired_end:
        R2_read_length_file = make_read_length_file(
                            args.fastqs[1], args.out_dir)
    
    # if bowtie2 index is tarball then unpack it
    if args.bowtie2_index_prefix_or_tar.endswith('.tar'):
        log.info('Unpacking bowtie2 index tar...')
        tar = args.bowtie2_index_prefix_or_tar
        # untar
        untar(tar, args.out_dir)
        bowtie2_index_prefix = os.path.join(
            args.out_dir, os.path.basename(strip_ext_tar(tar)))
        temp_files.append('{}.*.bt2'.format(
            bowtie2_index_prefix))
    else:
        bowtie2_index_prefix = args.bowtie2_index_prefix_or_tar

    # check if bowties indices are unpacked on out_dir
    chk_bowtie2_index(bowtie2_index_prefix)

    # bowtie2
    log.info('Running bowtie2...')
    if args.paired_end:
        bam, align_log = bowtie2_pe(
            args.fastqs[0], args.fastqs[1], 
            bowtie2_index_prefix,
            args.multimapping, args.score_min, args.nth,
            args.out_dir)
    else:
        bam, align_log = bowtie2_se(
            args.fastqs[0], 
            bowtie2_index_prefix,
            args.multimapping, args.score_min, args.nth,
            args.out_dir)

    # initialize multithreading
    log.info('Initializing multithreading...')
    num_process = min(2,args.nth)
    log.info('Number of threads={}.'.format(num_process))
    pool = multiprocessing.Pool(num_process)
    
    # samtools index
    log.info('Running samtools index...')
    ret_val1 = pool.apply_async(
        samtools_index, (bam, args.out_dir))

    # samtools flagstat qc
    log.info('Running samtools flagstat...')
    ret_val2 = pool.apply_async(
        samtools_flagstat, (bam, args.out_dir))

    bai = ret_val1.get()
    flagstat_qc = ret_val2.get()

    # close multithreading
    pool.close()
    pool.join()

    # remove temporary/intermediate files
    log.info('Removing temporary files...')
    rm_f(temp_files)

    log.info('All done.')

if __name__=='__main__':
    main()