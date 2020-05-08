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
                        help='Control depth ratio (between any two controls).')
    parser.add_argument('--ctl-depth-limit', type=int, default=200000000,
                        help='Control depth limit. If read depth of chosen control is '
                             'over this limit then such control should be subsampled.')
    parser.add_argument('--exp-ctl-depth-ratio-limit', type=float, default=5.0,
                        help='Exp vs. control depth ratio limit. ')
    parser.add_argument('--always-use-pooled-ctl', action="store_true",
                        help='Always use pooled control for all IP '
                             'replicates.')
    parser.add_argument('--out-tsv-basename', default='chosen_ctl.tsv', type=str,
                        help='Output TSV basename '
                             '(will be written on directory --out-dir). '
                             'This TSV file has chosen control index '
                             'per line (for each exp replicate).')
    parser.add_argument('--out-tsv-subsample-basename', default='chosen_ctl_subsample.tsv', type=str,
                        help='Output TSV subsample basename '
                             '(will be written on directory --out-dir). '
                             'This TSV file has number of reads to subsample control '
                             'per line (for each exp replicate). '
                             '0 means no subsampling for control.')
    parser.add_argument('--out-txt-subsample-pooled-basename', default='chosen_ctl_subsample_pooled.txt', type=str,
                        help='Output TXT subsample basename for pooled control'
                             '(will be written on directory --out-dir). '
                             'This TXT file has a single line for '
                             'number of reads to subsample pooled control control'
                             '0 means no subsampling for control.')
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
    num_rep = len(args.tas)
    num_ctl = len(args.ctl_tas)

    # num lines in tagaligns
    depths = [get_num_lines(ta) for ta in args.tas]
    # num lines in control tagaligns
    depths_ctl = [get_num_lines(ctl_ta) for ctl_ta in args.ctl_tas]
    depth_rep_pooled = sum(depths)
    depth_ctl_pooled = sum(depths_ctl)

    # make them dicts including -1 key (meaning pooled one)
    depths = dict(enumerate(depths))
    depths_ctl = dict(enumerate(depths_ctl))

    depths[-1] = depth_rep_pooled
    depths_ctl[-1] = depth_ctl_pooled

    ctl_ta_idx = [0]*num_rep
    if num_ctl == 1:
        # if only one control, use it for all replicates
        pass
    elif args.always_use_pooled_ctl:
        # if --always-use-pooled-ctl, then always use pooled control
        ctl_ta_idx = [-1]*num_rep
    else:
        # if multiple controls,
        # check # of lines in replicate/control tagaligns and
        # apply ctl_depth_ratio

        # make depths dicts including pooled ones

        # check every num lines in every pair of control tagaligns
        # if ratio of two entries in any pair > ctl_depth_ratio then
        # use pooled control for all
        use_pooled_ctl = False
        for i in range(num_ctl):
            for j in range(i+1, num_ctl):
                if depths_ctl[i]/float(depths_ctl[j]) > \
                        args.ctl_depth_ratio or \
                        depths_ctl[j]/float(depths_ctl[i]) > \
                        args.ctl_depth_ratio:
                    use_pooled_ctl = True
                    log.info(
                        'Number of reads in controls differ by a factor of {}.'
                        'Using pooled controls.'.format(
                            args.ctl_depth_ratio))
                    break

        if use_pooled_ctl:
            # use pooled control for all exp replicates
            ctl_ta_idx = [-1]*num_rep
        else:
            for i in range(num_rep):
                if i > num_ctl-1:
                    ctl_ta_idx[i] = -1  # use pooled control
                elif depths_ctl[i] < depths[i]:
                    log.info(
                        'Fewer reads in control {} than experiment replicate '
                        '{}. Using pooled control for replicate {}.'.format(
                            i+1, i+1, i+1))
                    ctl_ta_idx[i] = -1  # use pooled control
                else:
                    ctl_ta_idx[i] = i

    ctl_ta_subsample = [0] * num_rep
    ctl_ta_subsampled_pooled = 0
    if args.exp_ctl_depth_ratio_limit or args.ctl_depth_limit:
        # subsampling chosen control for each replicate
        for rep in range(num_rep):
            chosen_ctl = ctl_ta_idx[rep]
            depth = depths[rep]
            depth_ctl = depths_ctl[chosen_ctl]
            limit = int(max(depth * args.exp_ctl_depth_ratio_limit, args.ctl_depth_limit))
            if depth_ctl > limit:
                ctl_ta_subsample[rep] = limit

        # subsampling pooled control for pooled replicate
        limit = int(max(depth_rep_pooled * args.exp_ctl_depth_ratio_limit, args.ctl_depth_limit))
        if depth_ctl_pooled > limit:
            ctl_ta_subsampled_pooled = limit

    # for each replicate check
    log.info('Writing idx.txt...')
    out_txt = os.path.join(args.out_dir, args.out_tsv_basename)
    write_txt(out_txt, ctl_ta_idx)

    log.info('Writing subsample txt...')
    out_subsample_txt = os.path.join(args.out_dir, args.out_tsv_subsample_basename)
    write_txt(out_subsample_txt, ctl_ta_subsample)

    log.info('Writing subsample_pooled txt...')
    out_subsample_pooled_txt = os.path.join(args.out_dir, args.out_txt_subsample_pooled_basename)
    write_txt(out_subsample_pooled_txt, ctl_ta_subsampled_pooled)

    log.info('List all files in output directory...')
    ls_l(args.out_dir)

    log.info('All done.')


if __name__ == '__main__':
    main()
