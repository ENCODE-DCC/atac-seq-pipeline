#!/usr/bin/env python

# ENCODE DCC common functions python script
# Author: Jin Lee (leepc12@gmail.com)

import sys
import os
import re
import csv
import gzip
import logging
import subprocess
import signal

logging.basicConfig(
    format='[%(asctime)s %(levelname)s] %(message)s')
log = logging.getLogger(__name__)

BIG_INT = 99999999

def strip_ext_fastq(fastq):
    return re.sub(r'\.(fastq|fq|Fastq|Fq)\.gz$','',
                    str(fastq))

def strip_ext_bam(bam):
    return re.sub(r'\.(bam|Bam)$','',
                    str(bam))

def strip_ext_tar(tar):
    return re.sub(r'\.tar$','',
                    str(tar))

def strip_ext_ta(ta):
    return re.sub(r'\.(tagAlign|TagAlign|ta|Ta)\.gz$','',
                    str(ta))

def strip_ext_bed(bed):
    return re.sub(r'\.(bed|Bed)\.gz$','',
                    str(bed))

def strip_ext_npeak(npeak):
    return re.sub(r'\.(narrowPeak|NarrowPeak)\.gz$','',
                    str(npeak))

def strip_ext_rpeak(rpeak):
    return re.sub(r'\.(regionPeak|RegionPeak)\.gz$','',
                    str(npeak))

def strip_ext_gpeak(gpeak):
    return re.sub(r'\.(gappedPeak|GappedPeak)\.gz$','',
                    str(npeak))

def strip_ext_bpeak(bpeak):
    return re.sub(r'\.(broadPeak|BroadPeak)\.gz$','',
                    str(npeak))

def strip_ext_bigwig(bw):
    return re.sub(r'\.(bigwig|bw)$','',
                    str(bw))

def human_readable_number(num):
    for unit in ['','K','M','G','T','P']:
        if abs(num) < 1000:
            return '{}{}'.format(num, unit)
        num /= 1000
    return '{}{}'.format(num, 'E')

def human_readable_filesize(num):
    for unit in ['','KB','MB','GB','TB','PB']:
        if abs(num) < 1024.0:
            return '{}{}'.format(num, unit)
        num /= 1024.0
    return '{}{}'.format(num, 'EB')

def read_tsv(tsv):
    result = []
    with open(tsv,'r') as fp:        
        for row in csv.reader(fp,delimiter='\t'):
            result.append(row)
    return result

def write_tsv(tsv, arr): # arr must be list of list of string
    with open(tsv,'w') as fp:
        for i, a in enumerate(arr):
            s = '\t'.join(a) + ('\n' if i<len(arr)-1 else '')
            fp.write(s)

def mkdir_p(dirname):
    if not os.path.exists(dirname):
        os.makedirs(dirname)

def untar(tar, out_dir):
    cmd = 'tar xvf {} -C {}'.format(
        tar,
        out_dir)
    run_shell_cmd(cmd)

def get_num_lines(file):
    cmd = 'zcat -f {} | wc -l'.format(file)
    return int(run_shell_cmd(cmd))

def rm_rf(files):
    if files:
        if type(files)==list:
            run_shell_cmd('rm -rf {}'.format(' '.join(files)))
        else:
            run_shell_cmd('rm -rf {}'.format(files))

def run_shell_cmd(cmd): 
    try:
        p = subprocess.Popen(cmd, shell=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            universal_newlines=True,
            preexec_fn=os.setsid)
        log.info('run_shell_cmd: PID={}, CMD={}'.format(p.pid, cmd))
        ret = ''
        while True:
            line = p.stdout.readline()
            if line=='' and p.poll() is not None:
                break
            # log.debug('PID={}: {}'.format(p.pid,line.strip('\n')))
            print('PID={}: {}'.format(p.pid,line.strip('\n')))
            ret += line
            # sys.stdout.flush()
        p.communicate() # wait here
        if p.returncode > 0:
            raise subprocess.CalledProcessError(
                p.returncode, cmd)
        return ret.strip('\n')
    except:
        # kill all children processes when interrupted
        pgid = os.getpgid(p.pid)
        log.exception('Unknown exception caught. '+ \
            'Killing process group {}...'.format(pgid))
        # os.system("kill -{} -{}".format(signal.SIGKILL,pgid))
        os.killpg(os.getpgid(p.pid), signal.SIGKILL)
        p.terminate()

def samtools_index(bam, out_dir):
    prefix = os.path.join(out_dir,
        os.path.basename(strip_ext_bam(bam)))
    bai = '{}.bam.bai'.format(prefix)

    cmd = 'samtools index {}'.format(bam)
    run_shell_cmd(cmd)
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

def subsample_ta_se(ta, subsample, non_mito, out_dir):
    prefix = os.path.join(out_dir,
        os.path.basename(strip_ext_ta(ta)))
    subsampled_ta = '{}.{}.tagAlign.gz'.format(
        prefix, human_readable_number(subsample))

    cmd = 'zcat -f {} | '
    if non_mito:
        cmd += 'grep -v "chrM" | '
    cmd += 'shuf -n {} --random-source={} | '
    cmd += 'gzip -nc > {}'
    cmd = cmd.format(
        ta,
        subsample,
        ta,
        subsampled_ta)
    run_shell_cmd(cmd)
    return subsampled_ta

def subsample_ta_pe(ta, subsample, non_mito, out_dir):
    prefix = os.path.join(out_dir,
        os.path.basename(strip_ext_ta(ta)))
    ta = '{}.{}.tagAlign.gz'.format(
        prefix, human_readable_number(subsample))

    cmd = 'zcat -f {} | '
    if non_mito:
        cmd += 'grep -v "chrM" | '
    cmd += 'sed \'N;s/\\n/\\\t/\' | '
    cmd += 'shuf -n {} --random-source={} | '
    cmd += 'awk \'BEGIN{{OFS="\\t"}} | '
    cmd += '{{printf "%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n'
    cmd += '%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n",'
    cmd += '$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}}\' | '
    cmd += 'gzip -nc > {}'
    cmd = cmd.format(
        ta,
        subsample,
        ta,
        subsampled_ta)
    run_shell_cmd(cmd)
    return
