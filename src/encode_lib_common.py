#!/usr/bin/env python

# ENCODE DCC common functions
# Author: Jin Lee (leepc12@gmail.com)

import sys
import os
import re
import csv
import logging
import subprocess
import math
import signal
import time
import argparse

logging.basicConfig(
    format='[%(asctime)s %(levelname)s] %(message)s',
    stream=sys.stdout)
log = logging.getLogger(__name__)

BIG_INT = 99999999

# string/system functions


def strip_ext_fastq(fastq):
    return re.sub(r'\.(fastq|fq|Fastq|Fq)\.gz$', '',
                  str(fastq))


def strip_ext_bam(bam):
    return re.sub(r'\.(bam|Bam)$', '', str(bam))


def strip_ext_tar(tar):
    return re.sub(r'\.tar$', '', str(tar))


def strip_ext_ta(ta):
    return re.sub(r'\.(tagAlign|TagAlign|ta|Ta)\.gz$', '',
                  str(ta))


def strip_ext_bed(bed):
    return re.sub(r'\.(bed|Bed)\.gz$', '', str(bed))


def strip_ext_npeak(npeak):
    return re.sub(r'\.(narrowPeak|NarrowPeak)\.gz$', '',
                  str(npeak))


def strip_ext_rpeak(rpeak):
    return re.sub(r'\.(regionPeak|RegionPeak)\.gz$', '',
                  str(rpeak))


def strip_ext_gpeak(gpeak):
    return re.sub(r'\.(gappedPeak|GappedPeak)\.gz$', '',
                  str(gpeak))


def strip_ext_bpeak(bpeak):
    return re.sub(r'\.(broadPeak|BroadPeak)\.gz$', '',
                  str(bpeak))


def get_peak_type(peak):
    if strip_ext_npeak(peak) != peak:
        return 'narrowPeak'
    elif strip_ext_rpeak(peak) != peak:
        return 'regionPeak'
    elif strip_ext_bpeak(peak) != peak:
        return 'broadPeak'
    elif strip_ext_gpeak(peak) != peak:
        return 'gappedPeak'
    else:
        raise Exception(
            'Unsupported peak type for stripping extension {}'.format(peak))


def strip_ext_peak(peak):  # returns a tuple (peak_type, stripped_filename)
    peak_type = get_peak_type(peak)
    if peak_type == 'narrowPeak':
        return strip_ext_npeak(peak)
    elif peak_type == 'regionPeak':
        return strip_ext_rpeak(peak)
    elif peak_type == 'broadPeak':
        return strip_ext_bpeak(peak)
    elif peak_type == 'gappedPeak':
        return strip_ext_gpeak(peak)
    else:
        raise Exception(
            'Unsupported peak type for stripping '
            'extension {}'.format(peak))


def strip_ext_bigwig(bw):
    return re.sub(r'\.(bigwig|bw)$', '',
                  str(bw))


def strip_ext_gz(f):
    return re.sub(r'\.gz$', '', str(f))


def strip_ext(f, ext=''):
    if ext == '':
        ext = get_ext(f)
    return re.sub(r'\.({}|{}\.gz)$'.format(
        ext, ext), '', str(f))


def get_ext(f):
    f_wo_gz = re.sub(r'\.gz$', '', str(f))
    return f_wo_gz.split('.')[-1]


def human_readable_number(num):
    for unit in ['', 'K', 'M', 'G', 'T', 'P']:
        if abs(num) < 1000:
            return '{}{}'.format(num, unit)
        num = int(num/1000.0)
    return '{}{}'.format(num, 'E')


def human_readable_filesize(num):
    for unit in ['', 'KB', 'MB', 'GB', 'TB', 'PB']:
        if abs(num) < 1024.0:
            return '{}{}'.format(num, unit)
        num = int(num/1024.0)
    return '{}{}'.format(num, 'EB')


def read_tsv(tsv):
    result = []
    with open(tsv, 'r') as fp:
        for row in csv.reader(fp, delimiter='\t'):
            if len(row) == 0:
                row = ['']
            result.append(row)
    return result


def write_tsv(tsv, arr):  # arr must be list of list of something
    with open(tsv, 'w') as fp:
        for i, arr2 in enumerate(arr):
            s = '\t'.join([str(a) for a in arr2]) \
                + ('\n' if i < len(arr)-1 else '')
            fp.write(s)


def mkdir_p(dirname):
    if dirname and not os.path.exists(dirname):
        os.makedirs(dirname)


def untar(tar, out_dir):
    cmd = 'tar xvf {} --no-same-owner -C {}'.format(
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
        if type(files) == list:
            run_shell_cmd('rm -f {}'.format(' '.join(files)))
        else:
            run_shell_cmd('rm -f {}'.format(files))


def touch(f):
    run_shell_cmd('touch {}'.format(f))


def get_num_lines(f):
    cmd = 'zcat -f {} | wc -l'.format(f)
    return int(run_shell_cmd(cmd))


def assert_file_not_empty(f, help=''):
    if not os.path.exists(f):
        raise Exception('File does not exist ({}). Help: {}'.format(f, help))
    elif get_num_lines(f) == 0:
        raise Exception('File is empty ({}). Help: {}'.format(f, help))


def write_txt(f, s):
    with open(f, 'w') as fp:
        if type(s) != list:
            arr = [s]
        else:
            arr = s
        for a in arr:
            fp.write(str(a)+'\n')


def hard_link(f, link):  # hard-link 'f' to 'link'
    # UNIX only
    if os.path.abspath(f) == os.path.abspath(link):
        raise Exception('Trying to hard-link itself. {}'.format(f))
    os.link(f, link)
    return link


def make_hard_link(f, out_dir):  # hard-link 'f' to 'out_dir'/'f'
    # UNIX only
    if os.path.dirname(f) == os.path.dirname(out_dir):
        raise Exception('Trying to hard-link itself. {}'.format(f))
    linked = os.path.join(out_dir, os.path.basename(f))
    rm_f(linked)
    os.link(f, linked)
    return linked


def soft_link(f, link):  # soft-link 'f' to 'link'
    # UNIX only
    if os.path.abspath(f) == os.path.abspath(link):
        raise Exception('Trying to soft-link itself. {}'.format(f))
    os.symlink(f, link)
    return link


def make_soft_link(f, out_dir):  # soft-link 'f' to 'out_dir'/'f'
    # UNIX only
    if os.path.dirname(f) == os.path.dirname(out_dir):
        raise Exception('Trying to soft-link itself. {}'.format(f))
    linked = os.path.join(out_dir, os.path.basename(f))
    rm_f(linked)
    os.symlink(f, linked)
    return linked


def copy_f_to_f(f, dest):  # copy 'f' to 'out_dir'/'f'
    if os.path.abspath(f) == os.path.abspath(dest):
        raise Exception('Trying to copy to itself. {}'.format(f))
    cmd = 'cp -f {} {}'.format(f, dest)
    run_shell_cmd(cmd)
    return dest


def copy_f_to_dir(f, out_dir):  # copy 'f' to 'out_dir'/'f'
    out_dir = os.path.abspath(out_dir)
    if not os.path.isdir(out_dir):
        raise Exception('Invalid destination directory {}.'.format(out_dir))
    dest = os.path.join(out_dir, os.path.basename(f))
    return copy_f_to_f(f, dest)


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
    p = subprocess.Popen(
        ['/bin/bash', '-o', 'pipefail'],  # to catch error in pipe
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        universal_newlines=True,
        preexec_fn=os.setsid)  # to make a new process with a new PGID
    pid = p.pid
    pgid = os.getpgid(pid)
    log.info('run_shell_cmd: PID={}, PGID={}, CMD={}'.format(pid, pgid, cmd))
    stdout, stderr = p.communicate(cmd)
    rc = p.returncode
    err_str = 'PID={}, PGID={}, RC={}\nSTDERR={}\nSTDOUT={}'.format(
        pid, pgid, rc, stderr.strip(), stdout.strip())
    if rc:
        # kill all child processes
        try:
            os.killpg(pgid, signal.SIGKILL)
        except:
            pass
        finally:
            raise Exception(err_str)
    else:
        log.info(err_str)
    return stdout.strip('\n')

# math


def nCr(n, r):  # combination
    return math.factorial(n)/math.factorial(r)/math.factorial(n-r)


def infer_n_from_nC2(nC2):  # calculate n from nC2
    if nC2:
        n = 2
        while(nCr(n, 2) != nC2):
            n += 1
            if n > 99:
                raise argparse.ArgumentTypeError(
                    'Cannot infer n from nC2. ' +
                    'wrong number of peakfiles ' +
                    'in command line arg. (--peaks)?')
        return n
    else:
        return 1


def infer_pair_label_from_idx(n, idx, prefix='rep'):
    cnt = 0
    for i in range(n):
        for j in range(i+1, n):
            if idx == cnt:
                return '{}{}_vs_{}{}'.format(
                    prefix, i+1, prefix, j+1)
            cnt += 1
    raise argparse.ArgumentTypeError(
        'Cannot infer rep_id from n and idx.')


def which(executable):
    return run_shell_cmd('which {}'.format(executable))
