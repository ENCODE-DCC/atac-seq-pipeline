#!/usr/bin/env python

# ENCODE DCC choose control wrapper
# Author: Jin Lee (leepc12@gmail.com)

import sys
import os
import argparse
from encode_lib_common import (
    copy_f_to_f, get_num_lines, log, ls_l, mkdir_p, write_txt)


def parse_arguments():
    parser = argparse.ArgumentParser(
        prog='ENCODE DCC Choose control.',
        description='Choose appropriate control for each IP replicate.'
                    'ctl_for_repN.tagAlign.gz will be generated for each '
                    'IP replicate on --out-dir. '
                    'This outputs a file with integers '
                    '(chosen control index for each replicate per line).')
    parser.add_argument('--tas', type=str, nargs='+', required=True,
                        help='List of experiment TAG-ALIGN per IP replicate.')
    parser.add_argument('--ctl-tas', type=str, nargs='+', required=True,
                        help='List of control TAG-ALIGN per IP replicate.')
    parser.add_argument('--ta-pooled', type=str, nargs='*',
                        help='Pooled experiment TAG-ALIGN.')
    parser.add_argument('--ctl-ta-pooled', type=str, nargs='*',
                        help='Pooled control TAG-ALIGN.')
    parser.add_argument('--ctl-depth-ratio', type=float, required=True,
                        help='Control depth ratio.')
    parser.add_argument('--always-use-pooled-ctl', action="store_true",
                        help='Always use pooled control for all IP '
                             'replicates.')
    parser.add_argument('--out-tsv-basename', default='chosen_ctl.tsv', type=str,
                        help='Output TSV basename '
                             '(will be written on directory --out-dir). '
                             'This TSV file has chosen control index per line '
                             '(for each replicate).')
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


def main():
    # read params
    args = parse_arguments()
    log.info('Initializing and making output directory...')

    # make out_dir (root of all outputs)
    mkdir_p(args.out_dir)

    # reproducibility QC
    log.info('Choosing appropriate control for each IP replicate...')
    ctl_ta_idx = [0]*len(args.tas)
    if len(args.ctl_tas) == 1:
        # if only one control, use it for all replicates
        pass
    elif args.always_use_pooled_ctl:
        # if --always-use-pooled-ctl, then always use pooled control
        ctl_ta_idx = [-1]*len(args.tas)
    else:
        # if multiple controls,
        # check # of lines in replicate/control tagaligns and
        # apply ctl_depth_ratio

        # num lines in tagaligns
        nlines = [get_num_lines(ta) for ta in args.tas]
        # num lines in control tagaligns
        nlines_ctl = [get_num_lines(ctl_ta) for ctl_ta in args.ctl_tas]

        # check every num lines in every pair of control tagaligns
        # if ratio of two entries in any pair > ctl_depth_ratio then
        # use pooled control for all
        use_pooled_ctl = False
        for i in range(len(nlines_ctl)):
            for j in range(i+1, len(nlines_ctl)):
                if nlines_ctl[i]/float(nlines_ctl[j]) > \
                        args.ctl_depth_ratio or \
                        nlines_ctl[j]/float(nlines_ctl[i]) > \
                        args.ctl_depth_ratio:
                    use_pooled_ctl = True
                    log.info(
                        'Number of reads in controls differ by a factor of {}.'
                        'Using pooled controls.'.format(
                            args.ctl_depth_ratio))
                    break

        if use_pooled_ctl:
            # use pooled control for all exp replicates
            ctl_ta_idx = [-1]*len(args.tas)
        else:
            for i in range(len(args.tas)):
                if i > len(args.ctl_tas)-1:
                    ctl_ta_idx[i] = -1  # use pooled control
                elif nlines_ctl[i] < nlines[i]:
                    log.info(
                        'Fewer reads in control {} than experiment replicate '
                        '{}. Using pooled control for replicate {}.'.format(
                            i+1, i+1, i+1))
                    ctl_ta_idx[i] = -1  # use pooled control
                else:
                    ctl_ta_idx[i] = i

    log.info('Writing idx.txt...')
    out_txt = os.path.join(args.out_dir, args.out_tsv_basename)
    write_txt(out_txt, ctl_ta_idx)

    log.info('List all files in output directory...')
    ls_l(args.out_dir)

    log.info('All done.')


if __name__ == '__main__':
    main()
