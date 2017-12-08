#!/usr/bin/env python

# ENCODE DCC common functions
# Author: Jin Lee (leepc12@gmail.com)

import sys
import os
import re
import csv
import gzip
import logging
import subprocess
import math
import signal
import time

logging.basicConfig(
    format='[%(asctime)s %(levelname)s] %(message)s',
    stream=sys.stdout)
log = logging.getLogger(__name__)

BIG_INT = 99999999

# string/system functions

def strip_ext_fastq(fastq):
    return re.sub(r'\.(fastq|fq|Fastq|Fq)\.gz$','',
                    str(fastq))

def strip_ext_bam(bam):
    return re.sub(r'\.(bam|Bam)$','',str(bam))

def strip_ext_tar(tar):
    return re.sub(r'\.tar$','',str(tar))

def strip_ext_ta(ta):
    return re.sub(r'\.(tagAlign|TagAlign|ta|Ta)\.gz$','',
                    str(ta))

def strip_ext_bed(bed):
    return re.sub(r'\.(bed|Bed)\.gz$','',str(bed))

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

def strip_ext_gz(f):
    return re.sub(r'\.gz$','',str(f))

def strip_ext(f, ext=''):
    if ext=='': ext = get_ext(f)
    return re.sub(r'\.({}|{}\.gz)$'.format(
        ext,ext),'',str(f))

def get_ext(f):    
    f_wo_gz = re.sub(r'\.gz$','',str(f))    
    return f_wo_gz.split('.')[-1]

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
            if len(row)==0:
                row=['']
            result.append(row)
    return result

def write_tsv(tsv, arr): # arr must be list of list of something
    with open(tsv,'w') as fp:
        for i, arr2 in enumerate(arr):
            s = '\t'.join([str(a) for a in arr2]) \
                + ('\n' if i<len(arr)-1 else '')
            fp.write(s)

def mkdir_p(dirname):    
    if dirname and not os.path.exists(dirname):
        os.makedirs(dirname)

def untar(tar, out_dir):
    cmd = 'tar xvf {} -C {}'.format(
        tar,
        out_dir if out_dir else '.')
    run_shell_cmd(cmd)

def gunzip(f, suffix, out_dir):
    if not f.endswith('.gz'):
        raise Exception('Cannot gunzip a file without .gz extension.')
    gunzipped = os.path.join(out_dir,
        os.path.basename(strip_ext_gz(f)))
    if suffix:
        gunzipped += '.{}'.format(suffix)
    # cmd = 'gzip -cd {} > {}'.format(f, gunzipped)
    cmd = 'zcat -f {} > {}'.format(f, gunzipped)
    run_shell_cmd(cmd)
    return gunzipped

def ls_l(d):
    cmd = 'ls -l {}'.format(d)
    run_shell_cmd(cmd)

def rm_f(files):
    if files:
        if type(files)==list:
            run_shell_cmd('rm -f {}'.format(' '.join(files)))
        else:
            run_shell_cmd('rm -f {}'.format(files))

def touch(f):
    run_shell_cmd('touch {}'.format(f))

def get_num_lines(f):
    cmd = 'zcat -f {} | wc -l'.format(f)
    return int(run_shell_cmd(cmd))

def assert_file_not_empty(f):
    if get_num_lines(f)==0:
        raise Exception('File is empty. {}'.format(f))

def write_txt(f,s):
    with open(f,'w') as fp:
        if type(s)!=list: arr = [s]
        else: arr = s
        for a in arr: fp.write(str(a)+'\n')

def hard_link(f, link):  # hard-link 'f' to 'link'
    # UNIX only
    if os.path.abspath(f)==os.path.abspath(link):
        raise Exception('Trying to hard-link itself. {}'.format(f))    
    os.link(f, link)
    return link

def make_hard_link(f, out_dir): # hard-link 'f' to 'out_dir'/'f'
    # UNIX only
    if os.path.dirname(f)==os.path.dirname(out_dir):
        raise Exception('Trying to hard-link itself. {}'.format(f))
    linked = os.path.join(out_dir, os.path.basename(f))
    rm_f(linked)
    os.link(f, linked)
    return linked

def make_empty_file(filename, out_dir):
    f = os.path.join(out_dir, os.path.basename(filename))
    touch(f)
    return f

def now():
    return time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())

def pdf2png(pdf, out_dir):
    prefix = os.path.join(out_dir,
        os.path.basename(strip_ext(pdf)))
    png = '{}.png'.format(prefix)
    
    cmd = 'gs -dFirstPage=1 -dLastPage=1 -dTextAlphaBits=4 '
    cmd += '-dGraphicsAlphaBits=4 -r110x110 -dUseCropBox -dQUIET '
    cmd += '-dSAFER -dBATCH -dNOPAUSE -dNOPROMPT -sDEVICE=png16m '
    cmd += '-sOutputFile={} -r144 {}'
    cmd = cmd.format(
        png,
        pdf)
    run_shell_cmd(cmd)
    return png

def run_shell_cmd(cmd): 
    try:
        p = subprocess.Popen(cmd, shell=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            universal_newlines=True,
            preexec_fn=os.setsid)
        pid = p.pid
        pgid = os.getpgid(pid)
        log.info('run_shell_cmd: PID={}, CMD={}'.format(pid, cmd))
        ret = ''
        while True:
            line = p.stdout.readline()
            if line=='' and p.poll() is not None:
                break
            # log.debug('PID={}: {}'.format(pid,line.strip('\n')))
            print('PID={}: {}'.format(pid,line.strip('\n')))
            ret += line
        p.communicate() # wait here
        if p.returncode > 0:
            raise subprocess.CalledProcessError(
                p.returncode, cmd)
        return ret.strip('\n')
    except:
        # kill all child processes
        log.exception('Unknown exception caught. '+ \
            'Killing process group {}...'.format(pgid))
        os.killpg(pgid, signal.SIGKILL)
        p.terminate()
        raise Exception('Unknown exception caught. PID={}'.format(pid))

# math

def nCr(n,r): # combination
    return math.factorial(n)/math.factorial(r)/math.factorial(n-r)

def infer_n_from_nC2(nC2): # calculate n from nC2
    if nC2:
        n=2
        while(nCr(n,2)!=nC2):
            n += 1
            if n > 99:
                raise argparse.ArgumentTypeError(
                'Cannot infer n from nC2. '+ \
                'wrong number of peakfiles '+ \
                'in command line arg. (--peaks)?')
        return n
    else:
        return 1

def infer_pair_label_from_idx(n, idx, prefix='rep'):
    cnt = 0
    for i in range(n):
        for j in range(i+1,n):
            if idx==cnt:
                return '{}{}-{}{}'.format(
                    prefix, i+1, prefix, j+1)
            cnt += 1
    raise argparse.ArgumentTypeError(
        'Cannot infer rep_id from n and idx.')
