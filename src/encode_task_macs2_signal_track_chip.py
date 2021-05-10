#!/usr/bin/env python

# ENCODE DCC MACS2 signal track wrapper
# Author: Jin Lee (leepc12@gmail.com)

import sys
import os
import argparse
from encode_lib_common import (
    get_num_lines, log, ls_l, mkdir_p, rm_f, run_shell_cmd, strip_ext_ta,
    get_gnu_sort_param,
)
from encode_lib_genomic import (
    subsample_ta_se, subsample_ta_pe,
)


def parse_arguments():
    parser = argparse.ArgumentParser(prog='ENCODE DCC MACS2 signal track',
                                     description='')
    parser.add_argument('tas', type=str, nargs='+',
                        help='Path for TAGALIGN file (first) and control TAGALIGN file (second; optional).')
    parser.add_argument('--fraglen', type=int, required=True,
                        help='Fragment length.')
    parser.add_argument('--shift', type=int, default=0,
                        help='macs2 callpeak --shift.')
    parser.add_argument('--chrsz', type=str,
                        help='2-col chromosome sizes file.')
    parser.add_argument('--gensz', type=str,
                        help='Genome size (sum of entries in 2nd column of \
                            chr. sizes file, or hs for human, ms for mouse).')
    parser.add_argument('--pval-thresh', default=0.01, type=float,
                        help='P-Value threshold.')
    parser.add_argument('--ctl-subsample', default=0, type=int,
                        help='Subsample control to this read depth '
                             '(0: no subsampling).')
    parser.add_argument('--ctl-paired-end', action="store_true",
                        help='Paired-end control TA.')
    parser.add_argument('--mem-gb', type=float, default=4.0,
                        help='Max. memory for this job in GB. '
                        'This will be used to determine GNU sort -S (defaulting to 0.5 of this value). '
                        'It should be total memory for this task (not memory per thread).')
    parser.add_argument('--out-dir', default='', type=str,
                        help='Output directory.')
    parser.add_argument('--log-level', default='INFO',
                        choices=['NOTSET', 'DEBUG', 'INFO',
                                 'WARNING', 'CRITICAL', 'ERROR',
                                 'CRITICAL'],
                        help='Log level')
    args = parser.parse_args()
    if len(args.tas) == 1:
        args.tas.append('')
    log.setLevel(args.log_level)
    log.info(sys.argv)
    return args


def macs2_signal_track(ta, ctl_ta, chrsz, gensz, pval_thresh, shift, fraglen,
                       ctl_subsample, ctl_paired_end, mem_gb, out_dir):
    basename_ta = os.path.basename(strip_ext_ta(ta))
    if ctl_ta:
        if ctl_subsample:
            if ctl_paired_end:
                ctl_ta = subsample_ta_pe(
                    ctl_ta, ctl_subsample,
                    non_mito=False, mito_chr_name=None, r1_only=False,
                    out_dir=out_dir)
            else:
                ctl_ta = subsample_ta_se(
                    ctl_ta, ctl_subsample,
                    non_mito=False, mito_chr_name=None,
                    out_dir=out_dir)
        basename_ctl_ta = os.path.basename(strip_ext_ta(ctl_ta))
        basename_prefix = '{}_x_{}'.format(basename_ta, basename_ctl_ta)
        if len(basename_prefix) > 200:  # UNIX cannot have len(filename) > 255
            basename_prefix = '{}_x_control'.format(basename_ta)
    else:
        basename_prefix = basename_ta
    prefix = os.path.join(out_dir, basename_prefix)
    fc_bigwig = '{}.fc.signal.bigwig'.format(prefix)
    pval_bigwig = '{}.pval.signal.bigwig'.format(prefix)
    # temporary files
    fc_bedgraph = '{}.fc.signal.bedgraph'.format(prefix)
    fc_bedgraph_srt = '{}.fc.signal.srt.bedgraph'.format(prefix)
    pval_bedgraph = '{}.pval.signal.bedgraph'.format(prefix)
    pval_bedgraph_srt = '{}.pval.signal.srt.bedgraph'.format(prefix)

    temp_files = []

    run_shell_cmd(
        ' macs2 callpeak '
        '-t {ta} {ctl_param} -f BED -n {prefix} -g {gensz} -p {pval_thresh} '
        '--nomodel --shift {shiftsize} --extsize {extsize} --keep-dup all -B --SPMR'.format(
            ta=ta,
            ctl_param='-c {ctl_ta}'.format(ctl_ta=ctl_ta) if ctl_ta else '',
            prefix=prefix,
            gensz=gensz,
            pval_thresh=pval_thresh,
            shiftsize=0,
            extsize=fraglen,
        )
    )

    run_shell_cmd(
        'macs2 bdgcmp -t "{prefix}_treat_pileup.bdg" '
        '-c "{prefix}_control_lambda.bdg" '
        '--o-prefix "{prefix}" -m FE '.format(
            prefix=prefix,
        )
    )

    run_shell_cmd(
        'bedtools slop -i "{prefix}_FE.bdg" -g {chrsz} -b 0 | '
        'awk \'{{if ($3 != -1) print $0}}\' |'
        'bedClip stdin {chrsz} {fc_bedgraph}'.format(
            prefix=prefix,
            chrsz=chrsz,
            fc_bedgraph=fc_bedgraph,
        )
    )

    # sort and remove any overlapping regions in bedgraph by comparing two lines in a row
    run_shell_cmd(
        'LC_COLLATE=C sort -k1,1 -k2,2n {sort_param} {fc_bedgraph} | '
        'awk \'BEGIN{{OFS="\\t"}}{{if (NR==1 || NR>1 && (prev_chr!=$1 || '
        'prev_chr==$1 && prev_chr_e<=$2)) '
        '{{print $0}}; prev_chr=$1; prev_chr_e=$3;}}\' > {fc_bedgraph_srt}'.format(
            sort_param=get_gnu_sort_param(mem_gb * 1024 ** 3, ratio=0.5),
            fc_bedgraph=fc_bedgraph,
            fc_bedgraph_srt=fc_bedgraph_srt
        )
    )
    rm_f(fc_bedgraph)

    run_shell_cmd(
        'bedGraphToBigWig {fc_bedgraph_srt} {chrsz} {fc_bigwig}'.format(
            fc_bedgraph_srt=fc_bedgraph_srt,
            chrsz=chrsz,
            fc_bigwig=fc_bigwig,
        )
    )
    rm_f(fc_bedgraph_srt)

    # sval counts the number of tags per million in the (compressed) BED file
    sval = float(get_num_lines(ta))/1000000.0

    run_shell_cmd(
        'macs2 bdgcmp -t "{prefix}_treat_pileup.bdg" '
        '-c "{prefix}_control_lambda.bdg" '
        '--o-prefix {prefix} -m ppois -S {sval}'.format(
            prefix=prefix,
            sval=sval,
        )
    )

    run_shell_cmd(
        'bedtools slop -i "{prefix}_ppois.bdg" -g {chrsz} -b 0 | '
        'awk \'{{if ($3 != -1) print $0}}\' |'
        'bedClip stdin {chrsz} {pval_bedgraph}'.format(
            prefix=prefix,
            chrsz=chrsz,
            pval_bedgraph=pval_bedgraph,
        )
    )

    # sort and remove any overlapping regions in bedgraph by comparing two lines in a row
    run_shell_cmd(
        'LC_COLLATE=C sort -k1,1 -k2,2n {sort_param} {pval_bedgraph} | '
        'awk \'BEGIN{{OFS="\\t"}}{{if (NR==1 || NR>1 && (prev_chr!=$1 || '
        'prev_chr==$1 && prev_chr_e<=$2)) '
        '{{print $0}}; prev_chr=$1; prev_chr_e=$3;}}\' > {pval_bedgraph_srt}'.format(
            sort_param=get_gnu_sort_param(mem_gb * 1024 ** 3, ratio=0.5),
            pval_bedgraph=pval_bedgraph,
            pval_bedgraph_srt=pval_bedgraph_srt,
        )
    )
    rm_f(pval_bedgraph)

    run_shell_cmd(
        'bedGraphToBigWig {pval_bedgraph_srt} {chrsz} {pval_bigwig}'.format(
            pval_bedgraph_srt=pval_bedgraph_srt,
            chrsz=chrsz,
            pval_bigwig=pval_bigwig
        )
    )
    rm_f(pval_bedgraph_srt)

    # remove temporary files
    temp_files.append("{prefix}_*".format(prefix=prefix))
    rm_f(temp_files)

    return fc_bigwig, pval_bigwig


def main():
    # read params
    args = parse_arguments()

    log.info('Initializing and making output directory...')
    mkdir_p(args.out_dir)

    log.info('Calling peaks and generating signal tracks with MACS2...')
    fc_bigwig, pval_bigwig = macs2_signal_track(
        args.tas[0], args.tas[1], args.chrsz, args.gensz, args.pval_thresh,
        args.shift, args.fraglen, args.ctl_subsample, args.ctl_paired_end,
        args.mem_gb, args.out_dir)

    log.info('List all files in output directory...')
    ls_l(args.out_dir)

    log.info('All done.')


if __name__ == '__main__':
    main()
