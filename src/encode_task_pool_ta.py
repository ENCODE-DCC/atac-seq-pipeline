#!/usr/bin/env python

# ENCODE DCC TAGALIGN pooler wrapper
# Author: Jin Lee (leepc12@gmail.com)

import sys
import os
import argparse
from encode_lib_common import *

def parse_arguments():
    parser = argparse.ArgumentParser(prog='ENCODE DCC TAGALIGN pooler.',
                                        description='')
    parser.add_argument('tas', nargs='+', type=str,
                        help='List of TAGALIGNs to be pooled.')
    parser.add_argument('--out-dir', default='', type=str,
                            help='Output directory.')
    parser.add_argument('--log-level', default='INFO', 
                        choices=['NOTSET','DEBUG','INFO',
                            'WARNING','CRITICAL','ERROR','CRITICAL'],
                        help='Log level')
    args = parser.parse_args()

    log.setLevel(args.log_level)
    log.info(sys.argv)
    return args

def pool_ta(tas, out_dir):
    if len(tas)>1:
        prefix = os.path.join(out_dir,
            os.path.basename(strip_ext_ta(tas[0])))
        pooled_ta = '{}.pooled.tagAlign.gz'.format(prefix)

        cmd = 'zcat -f {} | gzip -nc > {}'
        cmd = cmd.format(
            ' '.join(tas),
            pooled_ta)
        run_shell_cmd(cmd)
        return pooled_ta
    else:
        return make_hard_link(tas[0], out_dir)

def main():
    # read params
    args = parse_arguments()

    log.info('Initializing and making output directory...')
    mkdir_p(args.out_dir)

    log.info('Pooling TAGALIGNs...')
    pool_ta(args.tas, args.out_dir)

    log.info('List all files in output directory...')
    ls_l(args.out_dir)

    log.info('All done.')

if __name__=='__main__':
    main()