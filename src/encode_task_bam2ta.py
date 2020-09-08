#!/usr/bin/env python

# ENCODE DCC BAM 2 TAGALIGN wrapper
# Author: Jin Lee (leepc12@gmail.com)

import sys
import os
import argparse
from encode_lib_common import (
    assert_file_not_empty, log, ls_l, mkdir_p, rm_f, run_shell_cmd,
    strip_ext_bam, strip_ext_ta)
from encode_lib_genomic import (
    samtools_name_sort, subsample_ta_pe, subsample_ta_se)


def parse_arguments():
    parser = argparse.ArgumentParser(prog='ENCODE DCC BAM 2 TAGALIGN.',
                                     description='')
    parser.add_argument('bam', type=str,
                        help='Path for BAM file.')
    parser.add_argument('--disable-tn5-shift', action="store_true",
                        help='Disable TN5 shifting for DNase-Seq.')
    parser.add_argument('--mito-chr-name', default='chrM',
                        help='Mito chromosome name.')
    parser.add_argument('--subsample', type=int, default=0,
                        help='Subsample TAGALIGN. \
                        This affects all downstream analysis.')
    parser.add_argument('--paired-end', action="store_true",
                        help='Paired-end BAM')
    parser.add_argument('--out-dir', default='', type=str,
                        help='Output directory.')
    parser.add_argument('--nth', type=int, default=1,
                        help='Number of threads to parallelize.')
    parser.add_argument('--mem-gb', type=float,
                        help='Max. memory for samtools sort in GB. '
                        'It should be total memory for this task (not memory per thread).')
    parser.add_argument('--log-level', default='INFO',
                        choices=['NOTSET', 'DEBUG', 'INFO',
                                 'WARNING', 'CRITICAL', 'ERROR',
                                 'CRITICAL'],
                        help='Log level')
    args = parser.parse_args()

    log.setLevel(args.log_level)
    log.info(sys.argv)
    return args


def bam2ta_se(bam, out_dir):
    prefix = os.path.join(out_dir,
                          os.path.basename(strip_ext_bam(bam)))
    ta = '{}.tagAlign.gz'.format(prefix)

    cmd = 'bedtools bamtobed -i {} | '
    cmd += 'awk \'BEGIN{{OFS="\\t"}}{{$4="N";$5="1000";print $0}}\' | '
    cmd += 'gzip -nc > {}'
    cmd = cmd.format(
        bam,
        ta)
    run_shell_cmd(cmd)
    return ta


def bam2ta_pe(bam, nth, out_dir):
    prefix = os.path.join(out_dir,
                          os.path.basename(strip_ext_bam(bam)))
    ta = '{}.tagAlign.gz'.format(prefix)
    # intermediate files
    bedpe = '{}.bedpe.gz'.format(prefix)
    nmsrt_bam = samtools_name_sort(bam, nth, out_dir)

    cmd1 = 'LC_COLLATE=C bedtools bamtobed -bedpe -mate1 -i {} | '
    # cmd1 += 'sort -k1,1 -k2,2n -k3,3n | '
    cmd1 += 'gzip -nc > {}'
    cmd1 = cmd1.format(
        nmsrt_bam,
        bedpe)
    run_shell_cmd(cmd1)
    rm_f(nmsrt_bam)

    cmd2 = 'zcat -f {} | '
    cmd2 += 'awk \'BEGIN{{OFS="\\t"}}'
    cmd2 += '{{printf "%s\\t%s\\t%s\\tN\\t1000\\t%s\\n'
    cmd2 += '%s\\t%s\\t%s\\tN\\t1000\\t%s\\n",'
    cmd2 += '$1,$2,$3,$9,$4,$5,$6,$10}}\' | '
    cmd2 += 'gzip -nc > {}'
    cmd2 = cmd2.format(
        bedpe,
        ta)
    run_shell_cmd(cmd2)
    rm_f(bedpe)
    return ta


def tn5_shift_ta(ta, out_dir):
    prefix = os.path.join(out_dir,
                          os.path.basename(strip_ext_ta(ta)))
    shifted_ta = '{}.tn5.tagAlign.gz'.format(prefix)

    cmd = 'zcat -f {} | '
    cmd += 'awk \'BEGIN {{OFS = "\\t"}}'
    cmd += '{{ if ($6 == "+") {{$2 = $2 + 4}} '
    cmd += 'else if ($6 == "-") {{$3 = $3 - 5}} print $0}}\' | '
    cmd += 'gzip -nc > {}'
    cmd = cmd.format(
        ta,
        shifted_ta)
    run_shell_cmd(cmd)
    return shifted_ta


def main():
    # read params
    args = parse_arguments()

    log.info('Initializing and making output directory...')
    mkdir_p(args.out_dir)

    # declare temp arrays
    temp_files = []  # files to deleted later at the end

    log.info('Converting BAM to TAGALIGN...')
    if args.paired_end:
        ta = bam2ta_pe(args.bam, args.nth, args.out_dir)
    else:
        ta = bam2ta_se(args.bam, args.out_dir)

    if args.subsample:
        log.info('Subsampling TAGALIGN...')
        if args.paired_end:
            subsampled_ta = subsample_ta_pe(
                ta, args.subsample, False,
                args.mito_chr_name, False, args.out_dir)
        else:
            subsampled_ta = subsample_ta_se(
                ta, args.subsample, False,
                args.mito_chr_name, args.out_dir)
        temp_files.append(ta)
    else:
        subsampled_ta = ta

    if args.disable_tn5_shift:
        shifted_ta = subsampled_ta
    else:
        log.info("TN5-shifting TAGALIGN...")
        shifted_ta = tn5_shift_ta(subsampled_ta, args.out_dir)
        temp_files.append(subsampled_ta)

    log.info('Checking if output is empty...')
    assert_file_not_empty(shifted_ta)

    log.info('Removing temporary files...')
    rm_f(temp_files)

    log.info('List all files in output directory...')
    ls_l(args.out_dir)

    log.info('All done.')


if __name__ == '__main__':
    main()
