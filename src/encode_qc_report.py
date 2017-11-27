#!/usr/bin/env python

# ENCODE DCC reporting module wrapper
# Author: Jin Lee (leepc12@gmail.com)

import sys
import os
import json
import base64
import argparse
from encode_common import *
from encode_common_log_parser import *
from encode_common_html import *
from collections import OrderedDict

def parse_arguments():
    parser = argparse.ArgumentParser(
                        prog='ENCODE DCC generate HTML report and QC JSON.',
                        description='')
    parser.add_argument('--name', type=str, default='Untitled',
                        help='Name of sample.')
    parser.add_argument('--desc', type=str, default='No description',
                        help='Description for sample.')
    parser.add_argument('--paired-end', action="store_true",
                        help='Paired-end sample.')
    parser.add_argument('--pipeline-type', type=str, required=True,
                        choices=['atac','dnase','tf','histone'],
                        help='Pipeline type.')
    parser.add_argument('--peak-caller', type=str, required=True,
                        help='Description for sample.')
    parser.add_argument('--idr-thresh', type=float, required=True,
                        help='IDR threshold.')
    parser.add_argument('--flagstat-qcs', type=str, nargs='*',
                        help='List of flagstat QC (raw BAM) files per replicate.')
    parser.add_argument('--nodup-flagstat-qcs', type=str, nargs='*',
                        help='List of flagstat QC (filtered BAM) files per replicate.')
    parser.add_argument('--dup-qcs', type=str, nargs='*',
                        help='List of dup QC files per replicate.')
    parser.add_argument('--pbc-qcs', type=str, nargs='*',
                        help='List of PBC QC files per replicate.')
    parser.add_argument('--xcor-plots', type=str, nargs='*',
                        help='List of cross-correlation QC plot files per replicate.')
    parser.add_argument('--xcor-scores', type=str, nargs='*',
                        help='List of cross-correlation QC score files per replicate.')
    parser.add_argument('--idr-plots', type=str, nargs='*',
                        help='List of IDR plot files per a pair of two replicates.')
    parser.add_argument('--idr-plots-pr', type=str, nargs='*',
                        help='List of IDR plot files per replicate.')
    parser.add_argument('--idr-plot-ppr', type=str, nargs='*',
                        help='IDR plot file for pooled pseudo replicate.')
    parser.add_argument('--frip-qcs', type=str, nargs='*',
                        help='List of FRiP score files per replicate.')
    parser.add_argument('--frip-qcs-pr1', type=str, nargs='*',
                        help='List of FRiP score files for 1st pseudo replicates per replicate.')
    parser.add_argument('--frip-qcs-pr2', type=str, nargs='*',
                        help='List of FRiP score files for 2nd pseudo replicates per replicate.')
    parser.add_argument('--frip-qc-pooled', type=str, nargs='*',
                        help='FRiP score file for pooled replicates.')
    parser.add_argument('--frip-qc-ppr1', type=str, nargs='*',
                        help='FRiP score file for 1st pooled pseudo replicates.')
    parser.add_argument('--frip-qc-ppr2', type=str, nargs='*',
                        help='FRiP score file for 2nd pooled pseudo replicates.')
    parser.add_argument('--frip-idr-qcs', type=str, nargs='*',
                        help='List of IDR FRiP score files per a pair of two replicates.')
    parser.add_argument('--frip-idr-qcs-pr', type=str, nargs='*',
                        help='List of IDR FRiP score files for pseudo replicates per replicate.')
    parser.add_argument('--frip-idr-qc-ppr', type=str, nargs='*',
                        help='IDR FRiP score file for pooled pseudo replicates.')
    parser.add_argument('--frip-overlap-qcs', type=str, nargs='*',
                        help='List of overlapping peak FRiP score files \
                            per a pair of two replicates.')
    parser.add_argument('--frip-overlap-qcs-pr', type=str, nargs='*',
                        help='List of overlapping peak FRiP score files \
                            for pseudo replicates per replicate.')
    parser.add_argument('--frip-overlap-qc-ppr', type=str, nargs='*',
                        help='Overlapping peak FRiP score file \
                            for pooled pseudo replicates.')
    parser.add_argument('--idr-reproducibility-qc', type=str, nargs='*',
                        help='IDR reproducibility QC file.')
    parser.add_argument('--overlap-reproducibility-qc', type=str, nargs='*',
                        help='Overlapping peak reproducibility QC file.')
    parser.add_argument('--ataqc-log', type=str, nargs='*',
                        help='ATAQC *_qc.txt.')
    parser.add_argument('--out-qc-html', default='qc.html', type=str,
                        help='Output QC report HTML file.')
    parser.add_argument('--out-qc-json', default='qc.json', type=str,
                        help='Output QC JSON file.')
    parser.add_argument('--log-level', default='INFO',
                        choices=['NOTSET','DEBUG','INFO',
                            'WARNING','CRITICAL','ERROR','CRITICAL'],
                        help='Log level')
    args = parser.parse_args()

    log.setLevel(args.log_level)
    log.info(sys.argv)
    return args

def get_full_name(pipeline_type):
    if pipeline_type=='atac':
        return 'ATAC-Seq'
    elif pipeline_type=='dnase':
        return 'DNase-Seq'
    elif pipeline_type=='tf':
        return 'TF ChIP-Seq'
    elif pipeline_type=='histone':
        return 'Histone ChIP-Seq'
    return pipeline_type

def main():
    # read params
    args = parse_arguments()

    log.info('Initializing...')
    html = ''
    json_all = OrderedDict()
    json_all['name']=args.name
    json_all['desc']=args.desc

    html += html_heading(1, args.name)
    html += html_paragraph(args.desc)
    html += html_paragraph('Report generated at {}'.format(now()))
    html += html_paragraph('Pipeline type: {}'.format(get_full_name(args.pipeline_type)))
    html += html_paragraph('Peak caller: {}'.format(args.peak_caller.upper()))

    log.info('Parsing QC logs...')
    if args.flagstat_qcs:
        html += html_heading(2, 'Flagstat QC (raw BAM)')
        json_objs = [parse_flagstat_qc(qc) for qc in args.flagstat_qcs]
        html += html_vert_table_multi_rep(json_objs, args.paired_end)
        json_all['flagstat_qc'] = json_objs

    if args.nodup_flagstat_qcs:
        html += html_heading(2, 'Flagstat QC (filtered BAM)')
        json_objs = [parse_flagstat_qc(qc) for qc in args.nodup_flagstat_qcs]
        html += html_vert_table_multi_rep(json_objs, args.paired_end)
        json_all['nodup_flagstat_qc'] = json_objs

    if args.dup_qcs:
        # check if file is empty (when filter.no_dup_removal is on)
        # if empty then skip
        if get_num_lines(args.dup_qcs[0]):
            html += html_heading(2, 'MarkDuplicate QC')
            json_objs = [parse_dup_qc(qc) for qc in args.dup_qcs]
            html += html_vert_table_multi_rep(json_objs, args.paired_end)
            json_all['dup_qc'] = json_objs

    if args.pbc_qcs:
        # check if file is empty (when filter.no_dup_removal is on)
        # if empty then skip
        if get_num_lines(args.pbc_qcs[0]):
            html += html_heading(2, 'Library complexity QC')
            json_objs = [parse_pbc_qc(qc) for qc in args.pbc_qcs]
            html += html_vert_table_multi_rep(json_objs, args.paired_end)
            html += html_help_pbc()
            json_all['pbc_qc'] = json_objs

    if args.xcor_plots:
        html += html_heading(2, 'Enrichment (strand cross-correlation measures) QC')
        json_objs = [parse_xcor_score(qc) for qc in args.xcor_scores]
        html += html_vert_table_multi_rep(json_objs, args.paired_end)
        html += html_help_xcor()
        for i, xcor_plot in enumerate(args.xcor_plots):
            html += html_embedded_png(xcor_plot, 'rep{}'.format(i+1), 60)
        json_all['xcor_score'] = json_objs

    # frip (Enrichment QC) for raw peaks
    if args.frip_qcs or args.frip_qcs_pr1 or args.frip_qcs_pr2 \
        or args.frip_qc_pooled or args.frip_qc_ppr1 or args.frip_qc_ppr2:
        html += html_heading(2, 'Enrichment QC (Fraction of reads in {} raw peaks)'.format(
            args.peak_caller.upper()))
        json_objs_frip = []
        row_header_frip = []
        true_rep_labels = ['rep{}'.format(i+1) for i, qc in enumerate(args.frip_qcs)]
        rep_pr1_labels = ['rep{}-pr1'.format(i+1) for i, qc in enumerate(args.frip_qcs_pr1)]
        rep_pr2_labels = ['rep{}-pr2'.format(i+1) for i, qc in enumerate(args.frip_qcs_pr2)]

        if args.frip_qcs:
            json_objs = [parse_frip_qc(qc) for qc in args.frip_qcs]
            json_objs_frip.extend(json_objs)
            row_header_frip.extend(true_rep_labels)
        if args.frip_qcs_pr1:
            json_objs = [parse_frip_qc(qc) for qc in args.frip_qcs_pr1]
            json_objs_frip.extend(json_objs)
            row_header_frip.extend(rep_pr1_labels)
        if args.frip_qcs_pr2:
            json_objs = [parse_frip_qc(qc) for qc in args.frip_qcs_pr2]
            json_objs_frip.extend(json_objs)
            row_header_frip.extend(rep_pr2_labels)
        if args.frip_qc_pooled:
            json_obj = parse_frip_qc(args.frip_qc_pooled[0])
            json_objs_frip.append(json_obj)
            row_header_frip.append('pooled')
        if args.frip_qc_ppr1:
            json_obj = parse_frip_qc(args.frip_qc_ppr1[0])
            json_objs_frip.append(json_obj)
            row_header_frip.append('ppr1')
        if args.frip_qc_ppr2:
            json_obj = parse_frip_qc(args.frip_qc_ppr2[0])
            json_objs_frip.append(json_obj)
            row_header_frip.append('ppr2')
        json_all['frip_qc'] = OrderedDict(
            zip(row_header_frip,json_objs_frip))
        html += html_vert_table_multi_rep(json_objs_frip,args.paired_end,row_header_frip)
        html += html_help_FRiP(args.peak_caller)

    # reproducibility_qc for naive-overlap and IDR
    if args.overlap_reproducibility_qc or args.idr_reproducibility_qc:
        html += html_heading(2, 'Reproducibility QC and Peak Detection Statistics')
        json_objs_reproducibility_qc = []
        row_header_reproducibility_qc = []
        if args.overlap_reproducibility_qc:
            json_obj = parse_reproducibility_qc(args.overlap_reproducibility_qc[0])
            json_all['overlap_reproducibility_qc'] = json_obj
            json_objs_reproducibility_qc.append(json_obj)
            row_header_reproducibility_qc.append('overlap')
        if args.idr_reproducibility_qc:
            json_obj = parse_reproducibility_qc(args.idr_reproducibility_qc[0])
            json_all['idr_reproducibility_qc'] = json_obj
            json_objs_reproducibility_qc.append(json_obj)
            row_header_reproducibility_qc.append('IDR')
        html += html_vert_table_multi_rep(
            json_objs_reproducibility_qc, 
            args.paired_end,
            row_header_reproducibility_qc)
        if args.overlap_reproducibility_qc:
            html += html_help_overlap()
        if args.idr_reproducibility_qc:
            html += html_help_idr(args.idr_thresh)

    # frip (Enrichment QC) for overlap
    if args.frip_overlap_qcs or args.frip_overlap_qcs_pr or args.frip_overlap_qc_ppr:
        html += html_heading(2, 'Enrichment QC (Fraction of reads in overlapping peaks)')
        json_objs_frip_overlap = []
        row_header_frip_overlap = []
        if args.frip_overlap_qcs:
            # infer num_rep from size of --frip-overlap-qcs
            num_rep = infer_n_from_nC2(len(args.frip_overlap_qcs))
            json_objs = [parse_frip_qc(qc) for qc in args.frip_overlap_qcs]
            json_objs_frip_overlap.extend(json_objs)
            row_header_frip_overlap.extend(
                [infer_pair_label_from_idx(num_rep, i) for i in range(len(args.frip_overlap_qcs))])
        if args.frip_overlap_qcs_pr:
            json_objs = [parse_frip_qc(qc) for qc in args.frip_overlap_qcs_pr]
            json_objs_frip_overlap.extend(json_objs)
            row_header_frip_overlap.extend(
                ['rep{}-pr'.format(i+1) for i in range(len(args.frip_overlap_qcs_pr))])
        if args.frip_overlap_qc_ppr:
            json_obj = parse_frip_qc(args.frip_overlap_qc_ppr[0])
            json_objs_frip_overlap.append(json_obj)
            row_header_frip_overlap.append('ppr')
        json_all['overlap_frip_qc'] = OrderedDict(zip(row_header_frip_overlap,json_objs_frip_overlap))
        html += html_vert_table_multi_rep(
            json_objs_frip_overlap,
            args.paired_end,
            row_header_frip_overlap)
        html += html_help_overlap_FRiP()

    # frip (Enrichment QC) for IDR
    if args.frip_idr_qcs or args.frip_idr_qcs_pr or args.frip_idr_qc_ppr:
        html += html_heading(2, 'Enrichment QC (Fraction of reads in IDR peaks)')
        json_objs_frip_idr = []
        row_header_frip_idr = []
        if args.frip_idr_qcs:
            # infer num_rep from size of --frip-idr-qcs
            num_rep = infer_n_from_nC2(len(args.frip_idr_qcs))
            json_objs = [parse_frip_qc(qc) for qc in args.frip_idr_qcs]
            json_objs_frip_idr.extend(json_objs)
            row_header_frip_idr.extend(
                [infer_pair_label_from_idx(num_rep, i) for i in range(len(args.frip_idr_qcs))])
        if args.frip_idr_qcs_pr:
            json_objs = [parse_frip_qc(qc) for qc in args.frip_idr_qcs_pr]
            json_objs_frip_idr.extend(json_objs)
            row_header_frip_idr.extend(
                ['rep{}-pr'.format(i+1) for i in range(len(args.frip_idr_qcs_pr))])
        if args.frip_idr_qc_ppr:
            json_obj = parse_frip_qc(args.frip_idr_qc_ppr[0])
            json_objs_frip_idr.append(json_obj)
            row_header_frip_idr.append('ppr')
        json_all['idr_frip_qc'] = OrderedDict(zip(row_header_frip_idr,json_objs_frip_idr))
        html += html_vert_table_multi_rep(
            json_objs_frip_idr,
            args.paired_end,
            row_header_frip_idr)
        html += html_help_idr_FRiP()

    # IDR plots
    if args.idr_plots or args.idr_plots_pr or args.idr_plot_ppr:
        html += html_heading(2, 'IDR (Irreproducible Discovery Rate) plots')
        if args.idr_plots:
            # infer num_rep from size of --idr-plots
            num_rep = infer_n_from_nC2(len(args.idr_plots))        
            # embed PNGs into HTML
            for i, png in enumerate(args.idr_plots):
                pair_label = infer_pair_label_from_idx(num_rep, i)
                html += html_embedded_png(png, pair_label)
        if args.idr_plots_pr:
            for i, png in enumerate(args.idr_plots_pr):
                html += html_embedded_png(png, 'rep{}-pr'.format(i+1))
        if args.idr_plot_ppr:
            png = args.idr_plot_ppr[0]
            html += html_embedded_png(png, 'ppr')

    if html:
        log.info('Creating HTML report...')
        write_txt(args.out_qc_html, html)

    log.info('Write JSON file...')
    write_txt(args.out_qc_json, json.dumps(json_all, indent=4))
    # b = json.loads(a, object_pairs_hook=OrderedDict)

    log.info('All done.')

if __name__=='__main__':
    main()
