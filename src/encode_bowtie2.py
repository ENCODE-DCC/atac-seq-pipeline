#!/usr/bin/env python

# ENCODE DCC bowtie2 wrapper
# Author: Jin Lee (leepc12@gmail.com), Daniel Kim

import sys
import os
import re
import argparse
from encode_common_genomic import *

def parse_arguments():
    parser = argparse.ArgumentParser(prog='ENCODE DCC bowtie2 aligner.',
                                        description='')
    parser.add_argument('bowtie2_index_prefix_or_tar', type=str,
                        help='Path for prefix (or a tarball .tar) \
                            for reference bowtie2 index. \
                            Prefix must be like [PREFIX].1.bt2*. \
                            Tar ball must be packed without compression \
                            and directory by using command line \
                            "tar cvf [TAR] [TAR_PREFIX].*.bt2*".')
    parser.add_argument('fastqs', nargs='+', type=str,
                        help='List of FASTQs (R1 and R2). \
                            FASTQs must be compressed with gzip (with .gz).')
    parser.add_argument('--paired-end', action="store_true",
                        help='Paired-end FASTQs.')
    parser.add_argument('--multimapping', default=4, type=int,
                        help='Multimapping reads (for bowtie2 -k(m+1). '
                             'This will be incremented in an actual bowtie2 command line'
                             'e.g. --multimapping 3 will be bowtie2 -k 4')
    parser.add_argument('--nth', type=int, default=1,
                        help='Number of threads to parallelize.')
    parser.add_argument('--out-dir', default='', type=str,
                            help='Output directory.')
    parser.add_argument('--log-level', default='INFO', 
                        choices=['NOTSET','DEBUG','INFO',
                            'WARNING','CRITICAL','ERROR','CRITICAL'],
                        help='Log level')
    args = parser.parse_args()

    # check if fastqs have correct dimension
    if args.paired_end and len(args.fastqs)!=2:
        raise argparse.ArgumentTypeError('Need 2 fastqs for paired end.')
    if not args.paired_end and len(args.fastqs)!=1:
        raise argparse.ArgumentTypeError('Need 1 fastq for single end.')

    log.setLevel(args.log_level)
    log.info(sys.argv)    
    return args

# WDL glob() globs in an alphabetical order
# so R1 and R2 can be switched, which results in an
# unexpected behavior of a workflow.
# so we already prepended merge_fastqs_'end'_ (R1 or R2) 
# to the basename of original filename in 'trim_adapter' task.
# now it's time to strip it.
def strip_merge_fastqs_prefix(fastq):
    return re.sub(r'^merge\_fastqs\_R\d\_','',str(fastq))

def make_read_length_file(fastq, out_dir):
    basename = os.path.basename(strip_ext_fastq(fastq))
    prefix = os.path.join(out_dir,
        strip_merge_fastqs_prefix(basename))
    txt = '{}.read_length.txt'.format(prefix)
    read_length = get_read_length(fastq)
    with open(txt,'w') as fp:
        fp.write(str(read_length))
    return txt

def bowtie2_se(fastq, ref_index_prefix, 
        multimapping, nth, out_dir):
    basename = os.path.basename(strip_ext_fastq(fastq))
    prefix = os.path.join(out_dir,
        strip_merge_fastqs_prefix(basename))        
    bam = '{}.bam'.format(prefix)
    align_log = '{}.align.log'.format(prefix)

    cmd = 'bowtie2 {} --mm --threads {} -x {} -U {} 2> {} '
    cmd += '| samtools view -Su /dev/stdin | samtools sort - {}'
    cmd = cmd.format(
        '-k {}'.format(multimapping+1) if multimapping else '',
        nth,
        ref_index_prefix,
        fastq,
        align_log,
        prefix)
    run_shell_cmd(cmd)

    cmd2 = 'cat {}'.format(align_log)
    run_shell_cmd(cmd2)
    return bam, align_log

def bowtie2_pe(fastq1, fastq2, ref_index_prefix, 
        multimapping, nth, out_dir):
    basename = os.path.basename(strip_ext_fastq(fastq1))
    prefix = os.path.join(out_dir,
        strip_merge_fastqs_prefix(basename))
    bam = '{}.bam'.format(prefix)
    bai = '{}.bam.bai'.format(prefix)
    align_log = '{}.align.log'.format(prefix)

    cmd = 'bowtie2 {} -X2000 --mm --threads {} -x {} '
    cmd += '-1 {} -2 {} 2>{} | '
    cmd += 'samtools view -Su /dev/stdin | samtools sort - {}'
    cmd = cmd.format(
        '-k {}'.format(multimapping+1) if multimapping else '',
        nth,
        ref_index_prefix,
        fastq1,
        fastq2,
        align_log,
        prefix)
    run_shell_cmd(cmd)

    cmd2 = 'cat {}'.format(align_log)
    run_shell_cmd(cmd2)
    return bam, align_log

def chk_bowtie2_index(prefix):    
    index_1 = '{}.1.bt2'.format(prefix)
    index_2 = '{}.1.bt2l'.format(prefix)
    if not (os.path.exists(index_1) or os.path.exists(index_2)):
        raise Exception("Bowtie2 index does not exists. "+
            "Prefix = {}".format(prefix))

def main():
    # read params
    args = parse_arguments()

    log.info('Initializing and making output directory...')
    mkdir_p(args.out_dir)

    # declare temp arrays
    temp_files = [] # files to deleted later at the end

    # generate read length file
    log.info('Generating read length file...')
    R1_read_length_file = make_read_length_file(
                            args.fastqs[0], args.out_dir)
    # if args.paired_end:
    #     R2_read_length_file = make_read_length_file(
    #                         args.fastqs[1], args.out_dir)
    
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
        temp_files.append('{}.*.bt2l'.format(
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
            args.multimapping, args.nth,
            args.out_dir)
    else:
        bam, align_log = bowtie2_se(
            args.fastqs[0], 
            bowtie2_index_prefix,
            args.multimapping, args.nth,
            args.out_dir)

    log.info('Running samtools index...')
    bai = samtools_index(bam, args.out_dir)

    log.info('Running samtools flagstat...')
    flagstat_qc = samtools_flagstat(bam, args.out_dir)

    log.info('Removing temporary files...')
    rm_f(temp_files)

    log.info('List all files in output directory...')
    ls_l(args.out_dir)

    log.info('All done.')

if __name__=='__main__':
    main()
