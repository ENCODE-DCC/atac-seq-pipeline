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

def strip_ext_npeak(npeak):
    return re.sub(r'\.(narrowPeak|NarrowPeak)\.gz$','',
                    str(npeak))

def strip_ext_bigwig(bw):
    return re.sub(r'\.(bigwig|bw)$','',
                    str(bw))

def read_tsv(tsv):
    result = []
    with open(tsv,'r') as fp:        
        for row in csv.reader(fp,delimiter='\t'):
            result.append(row)
    return result

def write_tsv(tsv, arr): # arr must be list of list of string
    with open(tsv,'w') as fp:
        for a in arr:
            fp.write('\t'.join(a)+'\n')            

def mkdir_p(dirname):
    if not os.path.exists(dirname):
        os.makedirs(dirname)

def untar(tar, out_dir):
    cmd = 'tar xvf {} -C {}'.format(
        tar,
        out_dir)
    run_shell_cmd(cmd)

def run_shell_cmd(cmd): 
    try:
        p = subprocess.Popen(cmd, shell=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            universal_newlines=True,
            preexec_fn=os.setsid)
        log.info('run_shell_cmd: PID={}, CMD={}'.format(p.pid, cmd))
        while True:
            line = p.stdout.readline()
            if line=='' and p.poll() is not None:
                break
            # log.debug('PID={}: {}'.format(p.pid,line.strip('\n')))
            print('PID={}: {}'.format(p.pid,line.strip('\n')))
            # sys.stdout.flush()        
        p.communicate() # wait here
        if p.returncode > 0:
            raise subprocess.CalledProcessError(
                p.returncode, cmd)
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
    cmd = 'samtools sort {} -o {} -@{}'.format(
        bam,
        srt_bam,        
        nth)
    run_shell_cmd(cmd)
    return srt_bam
