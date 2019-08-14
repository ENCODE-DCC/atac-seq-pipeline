#!/usr/bin/env python

# ENCODE frac mito
# Author: Jin Lee (leepc12@gmail.com)

import sys
import os
import argparse
from encode_lib_common import *
from encode_lib_log_parser import parse_flagstat_qc

def parse_arguments():
    parser = argparse.ArgumentParser(
        prog='ENCODE frac mito',
        description='Calculates fraction of mito reads')
    parser.add_argument('non_mito_samstat', type=str,
                        help='Path for SAMstats log file')
    parser.add_argument('mito_samstat', type=str,
                        help='Path for SAMstats log file (mito only)')
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

def frac_mito(non_mito_samstat, mito_samstat, out_dir):
    prefix = os.path.join(
        out_dir,
        os.path.basename(strip_ext(non_mito_samstat,
                                   'non_mito.samstats.qc')))
    frac_mito_qc = '{}.frac_mito.qc'.format(prefix)

    non_mito_samstat_dict = parse_flagstat_qc(non_mito_samstat)
    mito_samstat_dict = parse_flagstat_qc(mito_samstat)

    if 'mapped' in non_mito_samstat_dict:
        # backward compatibility (old key name was 'total')
        key_mapped = 'mapped'
    elif 'mapped_reads' in non_mito_samstat_dict:
        key_mapped = 'mapped_reads'
    Rn = non_mito_samstat_dict[key_mapped]

    if 'mapped' in mito_samstat_dict:
        # backward compatibility (old key name was 'total')
        key_mapped = 'mapped'
    elif 'mapped_reads' in mito_samstat_dict:
        key_mapped = 'mapped_reads'
    Rm = mito_samstat_dict[key_mapped]

    frac = float(Rm)/float(Rn + Rm)
    with open(frac_mito_qc, 'w') as fp:
        fp.write('non_mito_reads\t{}\n'.format(Rn))
        fp.write('mito_reads\t{}\n'.format(Rm))
        fp.write('frac_mito_reads\t{}\n'.format(frac))

    return frac_mito_qc

def main():
    # read params
    args = parse_arguments()
    log.info('Initializing and making output directory...')
    mkdir_p(args.out_dir)

    frac_mito_qc = frac_mito(args.non_mito_samstat,
                             args.mito_samstat,
                             args.out_dir)

    log.info('List all files in output directory...')
    ls_l(args.out_dir)

    log.info('All done.')

if __name__=='__main__':
    main()