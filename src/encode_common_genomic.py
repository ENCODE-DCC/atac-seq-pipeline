#!/usr/bin/env python

# ENCODE DCC common functions
# Author: Jin Lee (leepc12@gmail.com)

import os
from encode_common import *

def samtools_index(bam, out_dir=''):
    bai = '{}.bai'.format(bam)
    cmd = 'samtools index {}'.format(bam)
    run_shell_cmd(cmd)
    if os.path.abspath(out_dir)!= \
        os.path.abspath(os.path.dirname(bam)):
        cmd2 = 'mv {} {}'.format(bai, out_dir)
        return os.path.join(out_dir, os.path.basename(bai))
    else:
        return bai

def sambamba_index(bam, nth, out_dir=''):
    bai = '{}.bai'.format(bam)
    cmd = 'sambamba index {} -t {}'.format(bam, nth)
    run_shell_cmd(cmd)
    if os.path.abspath(out_dir)!= \
        os.path.abspath(os.path.dirname(bam)):
        cmd2 = 'mv {} {}'.format(bai, out_dir)
        return os.path.join(out_dir, os.path.basename(bai))
    else:
        return bai

def samtools_flagstat(bam, out_dir):
    prefix = os.path.join(out_dir,
        os.path.basename(strip_ext_bam(bam)))
    flagstat_qc = '{}.flagstat.qc'.format(prefix)

    cmd = 'samtools flagstat {} > {}'.format(
        bam,
        flagstat_qc)
    run_shell_cmd(cmd)
    return flagstat_qc

def sambamba_flagstat(bam, nth, out_dir):
    prefix = os.path.join(out_dir,
        os.path.basename(strip_ext_bam(bam)))
    flagstat_qc = '{}.flagstat.qc'.format(prefix)

    cmd = 'sambamba flagstat {} -t {} > {}'.format(
        bam,
        nth,
        flagstat_qc)
    run_shell_cmd(cmd)
    return flagstat_qc

def samtools_sort(bam, nth, out_dir):
    prefix = os.path.join(out_dir,
        os.path.basename(strip_ext_bam(bam)))
    srt_bam = '{}.srt.bam'.format(prefix)

    cmd = 'samtools sort {} -o {} -T {} -@ {}'.format(
        bam,
        srt_bam,
        prefix,
        nth)
    run_shell_cmd(cmd)
    return srt_bam

def sambamba_sort(bam, nth, out_dir):
    prefix = os.path.join(out_dir,
        os.path.basename(strip_ext_bam(bam)))
    srt_bam = '{}.srt.bam'.format(prefix)

    cmd = 'sambamba sort {} -o {} -t {}'.format(
        bam,
        srt_bam,
        nth)
    run_shell_cmd(cmd)
    return srt_bam

def samtools_name_sort(bam, nth, out_dir):
    prefix = os.path.join(out_dir,
        os.path.basename(strip_ext_bam(bam)))
    nmsrt_bam = '{}.nmsrt.bam'.format(prefix)

    cmd = 'samtools sort -n {} -o {} -T {} -@ {}'.format(
        bam,
        nmsrt_bam,
        prefix,
        nth)
    run_shell_cmd(cmd)
    return nmsrt_bam

def sambamba_name_sort(bam, nth, out_dir):
    prefix = os.path.join(out_dir,
        os.path.basename(strip_ext_bam(bam)))
    nmsrt_bam = '{}.nmsrt.bam'.format(prefix)

    cmd = 'sambamba sort -n {} -o {} -t {}'.format(
        bam,
        nmsrt_bam,
        nth)
    run_shell_cmd(cmd)
    return nmsrt_bam

def subsample_ta_se(ta, subsample, non_mito, out_dir):
    prefix = os.path.join(out_dir,
        os.path.basename(strip_ext_ta(ta)))
    ta_subsampled = '{}.{}{}.tagAlign.gz'.format(
        prefix,
        'no_chrM.' if non_mito else '',
        human_readable_number(subsample))

    cmd = 'zcat -f {} | '
    if non_mito:
        cmd += 'grep -v "chrM" | '
    cmd += 'shuf -n {} --random-source={} | '
    cmd += 'gzip -nc > {}'
    cmd = cmd.format(
        ta,
        subsample,
        ta,
        ta_subsampled)
    run_shell_cmd(cmd)
    return ta_subsampled

def subsample_ta_pe(ta, subsample, non_mito, r1_only, out_dir):
    prefix = os.path.join(out_dir,
        os.path.basename(strip_ext_ta(ta)))
    ta_subsampled = '{}.{}{}{}.tagAlign.gz'.format(
        prefix,
        'no_chrM.' if non_mito else '',
        'R1.' if r1_only else '',
        human_readable_number(subsample))

    cmd = 'zcat -f {} | '
    if non_mito:
        cmd += 'grep -v "chrM" | '
    cmd += 'sed \'N;s/\\n/\\t/\' | '
    cmd += 'shuf -n {} --random-source={} | '
    cmd += 'awk \'BEGIN{{OFS="\\t"}} '
    if r1_only:
        cmd += '{{printf "%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n'
        cmd += '",$1,$2,$3,$4,$5,$6}}\' | '
    else:
        cmd += '{{printf "%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n'
        cmd += '%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n",'
        cmd += '$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}}\' | '
    cmd += 'gzip -nc > {}'
    cmd = cmd.format(
        ta,
        subsample,
        ta,
        ta_subsampled)
    run_shell_cmd(cmd)
    return ta_subsampled

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
