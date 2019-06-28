#!/usr/bin/env python

# ENCODE annot_enrich (fraction of reads in annotated regions) wrapper
# Author: Daniel Kim, Jin Lee (leepc12@gmail.com)

import sys
import os
import argparse
import warnings
from encode_lib_common import run_shell_cmd, strip_ext_ta, ls_l, get_num_lines, log, logging, mkdir_p, rm_f
warnings.filterwarnings("ignore")

def parse_arguments():
    parser = argparse.ArgumentParser(prog='ENCODE annot_enrich (fraction of reads in annotated regions)')
    parser.add_argument('--ta', type=str, help='TAG-ALIGN file (from task bam2ta).')
    parser.add_argument('--dnase', type=str, help='DNase definition bed file.')
    parser.add_argument('--blacklist', type=str, help='Blacklist bed file.')
    parser.add_argument('--prom', type=str, help='Promoter definition bed file.')
    parser.add_argument('--enh', type=str, help='Enhancer definition bed file.')
    parser.add_argument('--out-dir', default='', type=str, help='Output directory.')
    parser.add_argument('--log-level', default='INFO', help='Log level',
                        choices=['NOTSET','DEBUG','INFO','WARNING','CRITICAL','ERROR','CRITICAL'])
    args = parser.parse_args()
    log.setLevel(args.log_level)
    log.info(sys.argv)
    return args

def get_fract_reads_in_regions(reads_bed, regions_bed):
    """Function that takes in bed file of reads and bed file of regions and
    gets fraction of reads sitting in said regions
    """
    # uses new run_shell_cmd
    cmd = "bedtools sort -i {}  | "
    cmd += "bedtools merge -i stdin | "
    cmd += "bedtools intersect -u -nonamecheck -a {} -b stdin | "
    cmd += "wc -l"
    #cmd += "bedtools intersect -c -nonamecheck -a stdin -b {} | "
    #cmd += "awk '{{ sum+=$4 }} END {{ print sum }}'"
    cmd = cmd.format(regions_bed, reads_bed)
    intersect_read_count = int(run_shell_cmd(cmd))
    total_read_count = get_num_lines(reads_bed)
    fract_reads = float(intersect_read_count) / total_read_count

    return intersect_read_count, fract_reads

def main():
    # read params
    args = parse_arguments()
    FINAL_BED = args.ta
    OUTPUT_PREFIX = os.path.join(
        args.out_dir,
        os.path.basename(strip_ext_ta(FINAL_BED)))

    DNASE = args.dnase if args.dnase and os.path.basename(args.dnase)!='null' else ''
    BLACKLIST = args.blacklist if args.blacklist and os.path.basename(args.blacklist)!='null' else ''
    PROM = args.prom if args.prom and os.path.basename(args.prom)!='null' else ''
    ENH = args.enh if args.enh and os.path.basename(args.enh)!='null' else ''

    result = []
    # Dnase regions
    if DNASE:
        reads_dnase, fract_dnase = get_fract_reads_in_regions(FINAL_BED, DNASE)
        result.append(('fraction_of_reads_in_universal_DHS_regions', str(reads_dnase), str(fract_dnase)))

    # Blacklist regions
    if BLACKLIST:
        reads_blacklist, \
            fract_blacklist = get_fract_reads_in_regions(FINAL_BED, BLACKLIST)
        result.append(('fraction_of_reads_in_blacklist_regions', str(reads_blacklist), str(fract_blacklist)))

    # Prom regions
    if PROM:
        reads_prom, fract_prom = get_fract_reads_in_regions(FINAL_BED, PROM)
        result.append(('fraction_of_reads_in_promoter_regions', str(reads_prom), str(fract_prom)))

    # Enh regions
    if ENH:
        reads_enh, fract_enh = get_fract_reads_in_regions(FINAL_BED, ENH)
        result.append(('fraction_of_reads_in_enhancer_regions', str(reads_enh), str(fract_enh)))

    annot_enrich_qc = OUTPUT_PREFIX + '.annot_enrich.qc'
    with open(annot_enrich_qc, 'w') as fp:
        for line in result:
            fp.write('\t'.join(line) + '\n')

    log.info('List all files in output directory...')
    ls_l(args.out_dir)

    log.info('All done.')

if __name__=='__main__':
    main()
