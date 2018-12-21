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

ATAQC_HTML_HEAD = \
'''
<head>
<title>QC report</title>
<style>
  .qc_table {
      font-family:"Lucida Sans Unicode", "Lucida Grande", Sans-Serif;
      font-size:12px;
      width:480px;
      text-align:left;
      border-collapse:collapse;
      margin:20px;
  }

  .qc_table th{
      font-size:14px;
      font-weight:normal;
      background:#8c1515;
      border-top:4px solid #700000;
      border-bottom:1px solid #fff;
      color:white;
      padding:8px;
  }

  .qc_table td{
      background:#f2f1eb;
      border-bottom:1px solid #fff;
      color:black;
      border-top:1px solid transparent;
      padding:8px;
  }

  .qc_table .fail{
      color:#ff0000;
      font-weight:bold;
  }
  </style>
</head>
'''

def parse_arguments():
    parser = argparse.ArgumentParser(
                        prog='ENCODE DCC generate HTML report and QC JSON.',
                        description='')
    parser.add_argument('--title', type=str, default='Untitled',
                        help='Title of sample.')
    parser.add_argument('--desc', type=str, default='No description',
                        help='Description for sample.')
    parser.add_argument('--genome', type=str,
                        help='Reference genome.')
    parser.add_argument('--pipeline-ver', type=str,
                        help='Pipeline version.')
    parser.add_argument('--multimapping', default=0, type=int,
                        help='Multimapping reads.')
    parser.add_argument('--paired-end', action="store_true",
                        help='Paired-end sample.')
    parser.add_argument('--pipeline-type', type=str, required=True,
                        choices=['atac','dnase','tf','histone'],
                        help='Pipeline type.')
    parser.add_argument('--peak-caller', type=str, required=True,
                        help='Description for sample.')
    parser.add_argument('--macs2-cap-num-peak', default=0, type=int,
                        help='Capping number of peaks by taking top N peaks for MACS.')
    parser.add_argument('--spp-cap-num-peak', default=0, type=int,
                        help='Capping number of peaks by taking top N peaks for SPP.')
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
    parser.add_argument('--ctl-flagstat-qcs', type=str, nargs='*',
                        help='List of flagstat QC (raw BAM) files per control.')
    parser.add_argument('--ctl-nodup-flagstat-qcs', type=str, nargs='*',
                        help='List of flagstat QC (filtered BAM) files per control.')
    parser.add_argument('--ctl-dup-qcs', type=str, nargs='*',
                        help='List of dup QC files per control.')
    parser.add_argument('--ctl-pbc-qcs', type=str, nargs='*',
                        help='List of PBC QC files per control.')
    parser.add_argument('--xcor-plots', type=str, nargs='*',
                        help='List of cross-correlation QC plot files per replicate.')
    parser.add_argument('--xcor-scores', type=str, nargs='*',
                        help='List of cross-correlation QC score files per replicate.')
    parser.add_argument('--jsd-plot', type=str, nargs='*',
                        help='Fingerprint JSD plot.')
    parser.add_argument('--jsd-qcs', type=str, nargs='*',
                        help='List of JSD qc files.')
    parser.add_argument('--idr-plots', type=str, nargs='*',
                        help='List of IDR plot files per a pair of two replicates.')
    parser.add_argument('--idr-plots-pr', type=str, nargs='*',
                        help='List of IDR plot files per replicate.')
    parser.add_argument('--idr-plot-ppr', type=str, nargs='*',
                        help='IDR plot file for pooled pseudo replicate.')

    parser.add_argument('--frip-macs2-qcs', type=str, nargs='*',
                        help='List of macs2 FRiP score files per replicate.')
    parser.add_argument('--frip-macs2-qcs-pr1', type=str, nargs='*',
                        help='List of macs2 FRiP score files for 1st pseudo replicates per replicate.')
    parser.add_argument('--frip-macs2-qcs-pr2', type=str, nargs='*',
                        help='List of macs2 FRiP score files for 2nd pseudo replicates per replicate.')
    parser.add_argument('--frip-macs2-qc-pooled', type=str, nargs='*',
                        help='macs2 FRiP score file for pooled replicates.')
    parser.add_argument('--frip-macs2-qc-ppr1', type=str, nargs='*',
                        help='macs2 FRiP score file for 1st pooled pseudo replicates.')
    parser.add_argument('--frip-macs2-qc-ppr2', type=str, nargs='*',
                        help='macs2 FRiP score file for 2nd pooled pseudo replicates.')
    parser.add_argument('--frip-spp-qcs', type=str, nargs='*',
                        help='List of spp FRiP score files per replicate.')
    parser.add_argument('--frip-spp-qcs-pr1', type=str, nargs='*',
                        help='List of spp FRiP score files for 1st pseudo replicates per replicate.')
    parser.add_argument('--frip-spp-qcs-pr2', type=str, nargs='*',
                        help='List of spp FRiP score files for 2nd pseudo replicates per replicate.')
    parser.add_argument('--frip-spp-qc-pooled', type=str, nargs='*',
                        help='spp FRiP score file for pooled replicates.')
    parser.add_argument('--frip-spp-qc-ppr1', type=str, nargs='*',
                        help='spp FRiP score file for 1st pooled pseudo replicates.')
    parser.add_argument('--frip-spp-qc-ppr2', type=str, nargs='*',
                        help='spp FRiP score file for 2nd pooled pseudo replicates.')

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
    parser.add_argument('--ataqc-txts', type=str, nargs='*',
                        help='ATAQC QC metrics JSON files *_qc.txt.')
    parser.add_argument('--ataqc-htmls', type=str, nargs='*',
                        help='ATAQC HTML reports *_qc.html.')
    parser.add_argument('--out-qc-html', default='qc.html', type=str,
                        help='Output QC report HTML file.')
    parser.add_argument('--out-qc-json', default='qc.json', type=str,
                        help='Output QC JSON file.')
    parser.add_argument('--qc-json-ref', type=str,
                        help='Reference QC JSON file to be compared to output QC JSON (developer\'s purpose only).')
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

def get_rep_labels(arr):
    return (['rep'+str(i+1) for i, a in enumerate(arr)])

def get_ctl_labels(arr):
    return (['ctl'+str(i+1) for i, a in enumerate(arr)])

def main():
    # read params
    args = parse_arguments()

    log.info('Initializing...')
    html = ATAQC_HTML_HEAD

    json_all = OrderedDict()
    json_all['general'] = OrderedDict()
    json_all['general']['date'] = now()
    json_all['general']['pipeline_ver'] = args.pipeline_ver
    json_all['general']['peak_caller'] = args.peak_caller
    json_all['general']['genome'] = args.genome
    json_all['general']['description'] = args.desc
    json_all['general']['title'] = args.title
    json_all['general']['paired_end'] = args.paired_end
    
    html += html_heading(1, args.title)
    html += html_paragraph(args.desc)
    html += html_paragraph('Pipeline version: {}'.format(args.pipeline_ver))
    html += html_paragraph('Report generated at {}'.format(now()))
    html += html_paragraph('Paired-end: {}'.format(args.paired_end))
    html += html_paragraph('Pipeline type: {}'.format(get_full_name(args.pipeline_type)))
    html += html_paragraph('Genome: {}'.format(args.genome))
    html += html_paragraph('Peak caller: {}'.format(args.peak_caller.upper()))

    # order of each chapter in report
    html_align = ''
    html_peakcall = ''
    html_enrich = ''
    html_etc = ''
    html_ataqc = ''

    log.info('Parsing QC logs...')

    if args.flagstat_qcs or args.ctl_flagstat_qcs:
        html_align += html_heading(2, 'Flagstat (raw BAM)')
        json_objs_all = []
        row_header = []
        if args.flagstat_qcs:
            json_objs = [parse_flagstat_qc(qc) for qc in args.flagstat_qcs]
            json_all['flagstat_qc'] = json_objs
            json_objs_all.extend(json_objs)
            row_header.extend(get_rep_labels(json_objs))
        if args.ctl_flagstat_qcs:
            json_objs = [parse_flagstat_qc(qc) for qc in args.ctl_flagstat_qcs]
            json_all['ctl_flagstat_qc'] = json_objs
            json_objs_all.extend(json_objs)
            row_header.extend(get_ctl_labels(json_objs))

        html_align += html_vert_table_multi_rep(json_objs_all, args.paired_end, row_header)

    if args.dup_qcs and get_num_lines(args.dup_qcs[0]) \
        or args.ctl_dup_qcs and get_num_lines(args.ctl_dup_qcs[0]):
        html_align += html_heading(2, 'Marking duplicates (filtered BAM)')
        html_align += html_help_filter(args.multimapping, args.paired_end)
        # check if file is empty (when filter.no_dup_removal is on)
        # if empty then skip
        json_objs_all = []
        row_header = []
        if args.dup_qcs and get_num_lines(args.dup_qcs[0]):
            json_objs = [parse_dup_qc(qc) for qc in args.dup_qcs]
            json_all['dup_qc'] = json_objs
            json_objs_all.extend(json_objs)
            row_header.extend(get_rep_labels(json_objs))
        if args.ctl_dup_qcs and get_num_lines(args.ctl_dup_qcs[0]):
            json_objs = [parse_dup_qc(qc) for qc in args.ctl_dup_qcs]
            json_all['ctl_dup_qc'] = json_objs
            json_objs_all.extend(json_objs)
            row_header.extend(get_ctl_labels(json_objs))
        html_align += html_vert_table_multi_rep(json_objs_all, args.paired_end, row_header)

    if args.pbc_qcs and get_num_lines(args.pbc_qcs[0]) \
        or args.ctl_pbc_qcs and get_num_lines(args.ctl_pbc_qcs[0]):
        html_align += html_heading(2, 'Library complexity (filtered non-mito BAM)')
        # check if file is empty (when filter.no_dup_removal is on)
        # if empty then skip
        json_objs_all = []
        row_header = []
        if args.pbc_qcs and get_num_lines(args.pbc_qcs[0]):
            json_objs = [parse_pbc_qc(qc) for qc in args.pbc_qcs]
            json_all['pbc_qc'] = json_objs
            json_objs_all.extend(json_objs)
            row_header.extend(get_rep_labels(json_objs))
        if args.ctl_pbc_qcs and get_num_lines(args.ctl_pbc_qcs[0]):
            json_objs = [parse_pbc_qc(qc) for qc in args.ctl_pbc_qcs]
            json_all['ctl_pbc_qc'] = json_objs
            json_objs_all.extend(json_objs)
            row_header.extend(get_ctl_labels(json_objs))
        html_align += html_vert_table_multi_rep(json_objs_all, args.paired_end, row_header)
        html_align += html_help_pbc()

    if args.nodup_flagstat_qcs or args.ctl_nodup_flagstat_qcs:
        html_align += html_heading(2, 'Flagstat (filtered/deduped BAM)')
        html_align += html_paragraph('Filtered and duplicates removed')
        json_objs_all = []
        row_header = []
        if args.nodup_flagstat_qcs:
            json_objs = [parse_flagstat_qc(qc) for qc in args.nodup_flagstat_qcs]
            json_all['nodup_flagstat_qc'] = json_objs
            json_objs_all.extend(json_objs)
            row_header.extend(get_rep_labels(json_objs))
        if args.ctl_nodup_flagstat_qcs:
            json_objs = [parse_flagstat_qc(qc) for qc in args.ctl_nodup_flagstat_qcs]
            json_all['ctl_nodup_flagstat_qc'] = json_objs
            json_objs_all.extend(json_objs)
            row_header.extend(get_ctl_labels(json_objs))
        html_align += html_vert_table_multi_rep(json_objs_all, args.paired_end, row_header)

    # IDR plots
    if args.idr_plots or args.idr_plots_pr or args.idr_plot_ppr:
        html_peakcall += html_heading(2, 'IDR (Irreproducible Discovery Rate) plots')
        if args.idr_plots:
            # infer num_rep from size of --idr-plots
            num_rep = infer_n_from_nC2(len(args.idr_plots))        
            # embed PNGs into HTML
            for i, png in enumerate(args.idr_plots):
                pair_label = infer_pair_label_from_idx(num_rep, i)
                html_peakcall += html_embedded_png(png, pair_label)
        if args.idr_plots_pr:
            for i, png in enumerate(args.idr_plots_pr):
                html_peakcall += html_embedded_png(png, 'rep{}-pr'.format(i+1))
        if args.idr_plot_ppr:
            png = args.idr_plot_ppr[0]
            html_peakcall += html_embedded_png(png, 'ppr')

    # reproducibility_qc for naive-overlap and IDR
    if args.overlap_reproducibility_qc or args.idr_reproducibility_qc:
        html_peakcall += html_heading(2, 'Reproducibility QC and peak detection statistics')
        if args.peak_caller=='macs2':
            html_peakcall += html_help_macs2(args.macs2_cap_num_peak)
        elif args.peak_caller=='spp':
            html_peakcall += html_help_spp(args.spp_cap_num_peak)
        else:
            raise Exception("Unsupported peak_caller")
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
        html_peakcall += html_vert_table_multi_rep(
            json_objs_reproducibility_qc, 
            args.paired_end,
            row_header_reproducibility_qc)
        if args.overlap_reproducibility_qc:
            html_peakcall += html_help_overlap()
        if args.idr_reproducibility_qc:
            html_peakcall += html_help_idr(args.idr_thresh)

    if args.xcor_plots and args.xcor_scores:
        html_enrich += html_heading(2, 'Strand cross-correlation measures')
        json_objs = [parse_xcor_score(qc) for qc in args.xcor_scores]
        html_enrich += html_paragraph('Performed on subsampled reads ({})'.format(
            human_readable_number(json_objs[0]['num_reads'])))

        html_enrich += html_vert_table_multi_rep(json_objs, args.paired_end)
        html_enrich += html_help_xcor()
        for i, xcor_plot in enumerate(args.xcor_plots):
            html_enrich += html_embedded_png(xcor_plot, 'rep{}'.format(i+1), 60)
        json_all['xcor_score'] = json_objs

    # frip (Enrichment QC) for MACS2 raw peaks
    if args.frip_macs2_qcs or args.frip_macs2_qcs_pr1 or args.frip_macs2_qcs_pr2 \
        or args.frip_macs2_qc_pooled or args.frip_macs2_qc_ppr1 or args.frip_macs2_qc_ppr2:
        # html_enrich += html_heading(2, 'Fraction of reads in {} raw peaks'.format(
        #     args.peak_caller.upper()))
        # html_enrich += html_help_macs2(args.macs2_cap_num_peak)
        json_objs_frip = []
        row_header_frip = []
        true_rep_labels = ['rep{}'.format(i+1) for i, qc in enumerate(args.frip_macs2_qcs)]
        rep_pr1_labels = ['rep{}-pr1'.format(i+1) for i, qc in enumerate(args.frip_macs2_qcs_pr1)]
        rep_pr2_labels = ['rep{}-pr2'.format(i+1) for i, qc in enumerate(args.frip_macs2_qcs_pr2)]

        if args.frip_macs2_qcs:
            json_objs = [parse_frip_qc(qc) for qc in args.frip_macs2_qcs]
            json_objs_frip.extend(json_objs)
            row_header_frip.extend(true_rep_labels)
        if args.frip_macs2_qcs_pr1:
            json_objs = [parse_frip_qc(qc) for qc in args.frip_macs2_qcs_pr1]
            json_objs_frip.extend(json_objs)
            row_header_frip.extend(rep_pr1_labels)
        if args.frip_macs2_qcs_pr2:
            json_objs = [parse_frip_qc(qc) for qc in args.frip_macs2_qcs_pr2]
            json_objs_frip.extend(json_objs)
            row_header_frip.extend(rep_pr2_labels)
        if args.frip_macs2_qc_pooled:
            json_obj = parse_frip_qc(args.frip_macs2_qc_pooled[0])
            json_objs_frip.append(json_obj)
            row_header_frip.append('pooled')
        if args.frip_macs2_qc_ppr1:
            json_obj = parse_frip_qc(args.frip_macs2_qc_ppr1[0])
            json_objs_frip.append(json_obj)
            row_header_frip.append('ppr1')
        if args.frip_macs2_qc_ppr2:
            json_obj = parse_frip_qc(args.frip_macs2_qc_ppr2[0])
            json_objs_frip.append(json_obj)
            row_header_frip.append('ppr2')
        json_all['frip_macs2_qc'] = OrderedDict(
            zip(row_header_frip,json_objs_frip))
        # html_enrich += html_vert_table_multi_rep(json_objs_frip,args.paired_end,row_header_frip)
        # html_enrich += html_help_FRiP(args.peak_caller)

    # frip (Enrichment QC) for SPP raw peaks
    if args.frip_spp_qcs or args.frip_spp_qcs_pr1 or args.frip_spp_qcs_pr2 \
        or args.frip_spp_qc_pooled or args.frip_spp_qc_ppr1 or args.frip_spp_qc_ppr2:
        # html_enrich += html_heading(2, 'Fraction of reads in {} raw peaks'.format(
        #     args.peak_caller.upper()))
        # html_peakcall += html_help_spp(args.spp_cap_num_peak)
        json_objs_frip = []
        row_header_frip = []
        true_rep_labels = ['rep{}'.format(i+1) for i, qc in enumerate(args.frip_spp_qcs)]
        rep_pr1_labels = ['rep{}-pr1'.format(i+1) for i, qc in enumerate(args.frip_spp_qcs_pr1)]
        rep_pr2_labels = ['rep{}-pr2'.format(i+1) for i, qc in enumerate(args.frip_spp_qcs_pr2)]

        if args.frip_spp_qcs:
            json_objs = [parse_frip_qc(qc) for qc in args.frip_spp_qcs]
            json_objs_frip.extend(json_objs)
            row_header_frip.extend(true_rep_labels)
        if args.frip_spp_qcs_pr1:
            json_objs = [parse_frip_qc(qc) for qc in args.frip_spp_qcs_pr1]
            json_objs_frip.extend(json_objs)
            row_header_frip.extend(rep_pr1_labels)
        if args.frip_spp_qcs_pr2:
            json_objs = [parse_frip_qc(qc) for qc in args.frip_spp_qcs_pr2]
            json_objs_frip.extend(json_objs)
            row_header_frip.extend(rep_pr2_labels)
        if args.frip_spp_qc_pooled:
            json_obj = parse_frip_qc(args.frip_spp_qc_pooled[0])
            json_objs_frip.append(json_obj)
            row_header_frip.append('pooled')
        if args.frip_spp_qc_ppr1:
            json_obj = parse_frip_qc(args.frip_spp_qc_ppr1[0])
            json_objs_frip.append(json_obj)
            row_header_frip.append('ppr1')
        if args.frip_spp_qc_ppr2:
            json_obj = parse_frip_qc(args.frip_spp_qc_ppr2[0])
            json_objs_frip.append(json_obj)
            row_header_frip.append('ppr2')
        json_all['frip_spp_qc'] = OrderedDict(
            zip(row_header_frip,json_objs_frip))
        # html_enrich += html_vert_table_multi_rep(json_objs_frip,args.paired_end,row_header_frip)
        # html_enrich += html_help_FRiP(args.peak_caller)

    # frip (Enrichment QC) for overlap
    if args.frip_overlap_qcs or args.frip_overlap_qcs_pr or args.frip_overlap_qc_ppr:
        html_enrich += html_heading(2, 'Fraction of reads in overlapping peaks')
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
        html_enrich += html_vert_table_multi_rep(
            json_objs_frip_overlap,
            args.paired_end,
            row_header_frip_overlap)
        html_enrich += html_help_overlap_FRiP()

    # frip (Enrichment QC) for IDR
    if args.frip_idr_qcs or args.frip_idr_qcs_pr or args.frip_idr_qc_ppr:
        html_enrich += html_heading(2, 'Fraction of reads in IDR peaks')
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
        html_enrich += html_vert_table_multi_rep(
            json_objs_frip_idr,
            args.paired_end,
            row_header_frip_idr)
        html_enrich += html_help_idr_FRiP()

    if args.jsd_plot:
        html_etc += html_heading(2, 'Fingerprint and Jensen-Shannon distance')
        json_objs = [parse_jsd_qc(qc) for qc in args.jsd_qcs]
        html_etc += html_vert_table_multi_rep(json_objs, args.paired_end)
        html_etc += html_embedded_png(args.jsd_plot[0], '', 40)
        json_all['jsd_qc'] = json_objs

    # ATAQC
    if args.ataqc_txts and args.ataqc_htmls:
        json_objs = [parse_ataqc_txt(txt) for txt in args.ataqc_txts]
        json_all['ataqc'] = json_objs
        html_ataqc += html_heading(2, 'Summary table')
        html_ataqc += html_vert_table_multi_rep(json_objs)
        for i, ataqc_html in enumerate(args.ataqc_htmls):
            html_ataqc += html_heading(2, 'Replicate {}'.format(i+1))
            html_ataqc += html_parse_body_from_file(ataqc_html).replace('<h2>ATAqC</h2>','')
            html_ataqc += '\n<br>'

    if html_align:
        html_align = html_heading(1, 'Alignment')+'<hr>'+html_align
    if html_peakcall:
        html_peakcall = html_heading(1, 'Peak calling')+'<hr>'+html_peakcall
    if html_enrich:
        html_enrich = html_heading(1, 'Enrichment')+'<hr>'+html_enrich
    if html_etc:
        html_etc = html_heading(1, 'Other quality metrics')+'<hr>'+html_etc
    if html_ataqc:
        html_ataqc = html_heading(1, 'ATAQC')+'<hr>'+html_ataqc

    html += html_align + html_peakcall + html_enrich + html_etc + html_ataqc

    log.info('Creating HTML report...')
    write_txt(args.out_qc_html, html)

    log.info('Converting format of qc.json...')
    # convert format of json_all
    # old: [](zero-based array) 0=rep1, 1=rep2, ...
    # new: obj (.rep1, .rep2, ...)
    json_all_new_format = OrderedDict()
    for key in json_all:
        val = json_all[key]
        if type(val)==list:
            json_all_new_format[key] = OrderedDict()
            for i, rep in enumerate(val):
                json_all_new_format[key]["rep{}".format(i+1)] = val[i]
        else:
            json_all_new_format[key] = json_all[key]
    log.info('Write JSON file...')
    write_txt(args.out_qc_json, json.dumps(json_all_new_format, indent=4))
    # b = json.loads(a, object_pairs_hook=OrderedDict)

    log.info('Comparing JSON file with reference...')
    if args.qc_json_ref:
        # exclude general section from comparing 
        #  because it includes metadata like date, pipeline_ver, ...
        # we just want to compare quality metric values only
        json_all_new_format.pop("general")
        with open(args.qc_json_ref,'r') as fp:
            json_ref = json.load(fp, object_pairs_hook=OrderedDict)
            json_ref.pop("general")
            match_qc_json_ref = json_all_new_format==json_ref
    else:
        match_qc_json_ref = False
    run_shell_cmd('echo {} > qc_json_ref_match.txt'.format(match_qc_json_ref))

    log.info('All done.')

if __name__=='__main__':
    main()
