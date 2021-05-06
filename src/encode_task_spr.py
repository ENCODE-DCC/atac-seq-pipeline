#!/usr/bin/env python

# ENCODE DCC pseudo replicator wrapper
# Author: Jin Lee (leepc12@gmail.com)

import sys
import os
import argparse
from encode_lib_common import (
    assert_file_not_empty, get_num_lines, log, ls_l, mkdir_p, rm_f,
    run_shell_cmd, strip_ext_ta)


def parse_arguments():
    parser = argparse.ArgumentParser(
        prog='ENCODE DCC pseudo replicator.')
    parser.add_argument('ta', type=str,
                        help='Path for TAGALIGN file.')
    parser.add_argument('--paired-end', action="store_true",
                        help='Paired-end TAGALIGN.')
    parser.add_argument('--pseudoreplication-random-seed',
                        type=int, default=0,
                        help='Set it to 0 to use file\'s size (in bytes) as random seed.'
                             'Otherwise this seed will be used for GNU shuf --random-source=sha256(seed).'
                             'It is useful when random seed based on input file size does not work.')
    parser.add_argument('--out-dir', default='', type=str,
                        help='Output directory.')
    parser.add_argument('--log-level', default='INFO',
                        choices=['NOTSET', 'DEBUG', 'INFO',
                                 'WARNING', 'CRITICAL', 'ERROR',
                                 'CRITICAL'],
                        help='Log level')
    args = parser.parse_args()

    log.setLevel(args.log_level)
    log.info(sys.argv)
    return args


def spr_se(ta, pseudoreplication_random_seed, out_dir):
    prefix = os.path.join(out_dir,
                          os.path.basename(strip_ext_ta(ta)))
    tmp_pr1 = '{}.00'.format(prefix)
    tmp_pr2 = '{}.01'.format(prefix)
    ta_pr1 = '{}.pr1.tagAlign.gz'.format(prefix)
    ta_pr2 = '{}.pr2.tagAlign.gz'.format(prefix)
    nlines = int((get_num_lines(ta)+1)/2)

    if pseudoreplication_random_seed == 0:
        random_seed = os.path.getsize(ta)
        log.info(
            'Using input file\'s size {seed} as random seed for pseudoreplication.'.format(seed=seed)
        )
    else:
        random_seed = pseudoreplication_random_seed
        log.info(
            'Using a fixed integer {seed} as random seed for pseudoreplication.'.format(seed=seed)
        )

    # bash-only
    run_shell_cmd(
        'zcat {ta} | shuf --random-source=<(openssl enc '
        '-aes-256-ctr -pass pass:{random_seed} '
        '-nosalt </dev/zero 2>/dev/null) | '
        'split -d -l {nlines} - {prefix}.'.format(
            ta=ta,
            random_seed=random_seed,
            nlines=nlines,
            prefix=prefix,
        )
    )

    run_shell_cmd('gzip -nc {tmp_pr1} > {ta_pr1}'.format(tmp_pr1=tmp_pr1, ta_pr1=ta_pr1))
    run_shell_cmd('gzip -nc {tmp_pr2} > {ta_pr2}'.format(tmp_pr2=tmp_pr2, ta_pr2=ta_pr2))

    rm_f([tmp_pr1, tmp_pr2])
    return ta_pr1, ta_pr2


def spr_pe(ta, pseudoreplication_random_seed, out_dir):
    prefix = os.path.join(out_dir,
                          os.path.basename(strip_ext_ta(ta)))
    tmp_pr1 = '{}.00'.format(prefix)
    tmp_pr2 = '{}.01'.format(prefix)
    ta_pr1 = '{}.pr1.tagAlign.gz'.format(prefix)
    ta_pr2 = '{}.pr2.tagAlign.gz'.format(prefix)
    nlines = int((get_num_lines(ta)/2+1)/2)

    if pseudoreplication_random_seed == 0:
        random_seed = os.path.getsize(ta)
        log.info(
            'Using input file\'s size {seed} as random seed for pseudoreplication.'.format(seed=seed)
        )
    else:
        random_seed = pseudoreplication_random_seed
        log.info(
            'Using a fixed integer {seed} as random seed for pseudoreplication.'.format(seed=seed)
        )

    # bash-only
    run_shell_cmd(
        'zcat -f {ta} | sed \'N;s/\\n/\\t/\' | '
        'shuf --random-source=<(openssl enc -aes-256-ctr '
        '-pass pass:${random_seed} -nosalt </dev/zero 2>/dev/null) | '
        'split -d -l {nlines} - {prefix}.'.format(
            ta=ta,
            random_seed=random_seed,
            nlines=nlines,
            prefix=prefix,
        )
    )

    run_shell_cmd(
        'zcat -f {tmp_pr1} | '
        'awk \'BEGIN{{OFS="\\t"}} '
        '{{printf "%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n'
        '%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n",'
        '$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}}\' | '
        'gzip -nc > {ta_pr1}'.format(
            tmp_pr1=tmp_pr1,
            ta_pr1=ta_pr1,
        )
    )

    run_shell_cmd(
        'zcat -f {tmp_pr2} | '
        'awk \'BEGIN{{OFS="\\t"}} '
        '{{printf "%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n'
        '%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n",'
        '$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}}\' | '
        'gzip -nc > {ta_pr2}'.format(
            tmp_pr2=tmp_pr2,
            ta_pr2=ta_pr2,
        )
    )

    rm_f([tmp_pr1, tmp_pr2])
    return ta_pr1, ta_pr2


def main():
    # read params
    args = parse_arguments()
    log.info('Initializing and making output directory...')
    mkdir_p(args.out_dir)

    log.info('Making self-pseudo replicates...')
    if args.paired_end:
        ta_pr1, ta_pr2 = spr_pe(
            args.ta, args.pseudoreplication_random_seed, args.out_dir,
        )
    else:
        ta_pr1, ta_pr2 = spr_se(
            args.ta, args.pseudoreplication_random_seed, args.out_dir,
        )

    log.info('List all files in output directory...')
    ls_l(args.out_dir)

    log.info('Checking if output is empty...')
    assert_file_not_empty(ta_pr1)
    assert_file_not_empty(ta_pr2)

    log.info('All done.')


if __name__ == '__main__':
    main()
