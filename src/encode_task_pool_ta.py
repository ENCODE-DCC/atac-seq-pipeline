#!/usr/bin/env python

# ENCODE DCC TAGALIGN pooler wrapper
# Author: Jin Lee (leepc12@gmail.com)

import sys
import os
import argparse
from encode_lib_common import (
    log, ls_l, make_hard_link, mkdir_p, run_shell_cmd, strip_ext_ta)


def parse_arguments():
    parser = argparse.ArgumentParser(prog='ENCODE DCC TAGALIGN pooler.',
                                     description='')
    parser.add_argument('tas', nargs='+', type=str,
                        help='List of TAGALIGNs to be pooled.')
    parser.add_argument('--prefix', type=str,
                        help='Basename prefix.')
    parser.add_argument('--out-dir', default='', type=str,
                        help='Output directory.')
    parser.add_argument('--col',
                        help='Number of columns to keep in a pooled TAGALIGN. '
                             'Keep all columns if not defined.')
    parser.add_argument('--log-level', default='INFO',
                        choices=['NOTSET', 'DEBUG', 'INFO',
                                 'WARNING', 'CRITICAL', 'ERROR',
                                 'CRITICAL'],
                        help='Log level')
    args = parser.parse_args()

    log.setLevel(args.log_level)
    log.info(sys.argv)
    return args


def pool_ta(tas, col, basename_prefix, out_dir):
    if len(tas) > 1:
        if basename_prefix is None:
            prefix = os.path.join(out_dir,'basename_prefix')
        else:
            prefix = os.path.join(out_dir,
                              os.path.basename(strip_ext_ta(tas[0])))
        pooled_ta = '{}.pooled.tagAlign.gz'.format(prefix)

        cmd = 'zcat -f {} | '
        if col is not None:
            cmd += 'cut -f 1-{} | '.format(col)
        cmd += 'gzip -nc > {}'
        cmd = cmd.format(
            ' '.join(tas),
            pooled_ta)
        run_shell_cmd(cmd)
        return pooled_ta
    else:
        raise Exception('Needs at least two TAs (or BEDs) to be pooled.')


def main():
    # read params
    args = parse_arguments()

    log.info('Initializing and making output directory...')
    mkdir_p(args.out_dir)

    log.info('Pooling TAGALIGNs...')
    pool_ta(args.tas, args.col, args.prefix, args.out_dir)

    log.info('List all files in output directory...')
    ls_l(args.out_dir)

    log.info('All done.')


if __name__ == '__main__':
    main()
