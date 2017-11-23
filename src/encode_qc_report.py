#!/usr/bin/env python

# ENCODE DCC reporting module wrapper
# Author: Jin Lee (leepc12@gmail.com)

import sys
import os
import json
import base64
import argparse
from encode_common import *
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
    parser.add_argument('--frip-idr-qcs', type=str, nargs='*',
                        help='List of IDR FRiP score files \
                            per a pair of two replicates.')
    parser.add_argument('--frip-idr-qcs-pr', type=str, nargs='*',
                        help='List of IDR FRiP score files \
                            for pseudo replicates per replicate.')
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

def get_long_keyname(key, paired_end=False):
    short_to_long = {
        'total' : 'Total',
        'total_qc_failed' : 'Total(QC-failed)',
        'duplicates' : 'Dupes',
        'duplicates_qc_failed' : 'Dupes(QC-failed)',
        'mapped' : 'Mapped',
        'mapped_qc_failed' : 'Mapped(QC-failed)',
        'mapped_pct' : '% Mapped',
        'paired' : 'Paired',
        'paired_qc_failed' : 'Paired(QC-failed)',
        'read1' : 'Read1',
        'read1_qc_failed' : 'Read1(QC-failed)',
        'read2' : 'Read2',
        'read2_qc_failed' : 'Read2(QC-failed)',
        'paired_properly' : 'Properly Paired',
        'paired_properly_qc_failed' : 'Properly Paired(QC-failed)',
        'paired_properly_pct' : '% Properly Paired',
        'with_itself' : 'With itself',
        'with_itself_qc_failed' : 'With itself(QC-failed)',
        'singletons' : 'Singletons',
        'singletons_qc_failed' : 'Singletons(QC-failed)',
        'singletons_pct' : '% Singleton',
        'diff_chroms' : 'Diff. Chroms',
        'diff_chroms_qc_failed' : 'Diff. Chroms (QC-failed)',
        'unpaired_reads' : 'Unpaired Reads',
        'paired_reads' : 'Paired Reads',
        'unmapped_reads' : 'Unmapped Reads',
        'unpaired_dupes' : 'Unpaired Dupes',
        'paired_dupes' : 'Paired Dupes',
        'paired_opt_dupes' : 'Paired Opt. Dupes',
        'dupes_pct' : '% Dupes/100',
        'num_reads' : 'Reads',
        'est_frag_len' : 'Est. Fragment Len.',
        'corr_est_frag_len' : 'Corr. Est. Fragment Len.',
        'phantom_peak' : 'Phantom Peak',
        'corr_phantom_peak' : 'Corr. Phantom Peak',
        'argmin_corr' : 'Argmin. Corr.',
        'min_corr' : 'Min. Corr.',
        'NSC' : 'NSC',
        'RSC' : 'RSC',
        'Nt' : 'Nt',
        'N1' : 'N1',
        'N2' : 'N2',
        'Np' : 'Np',
        'N_opt' : 'N optimal',
        'N_consv' : 'N conservative',
        'opt_set' : 'Optimal Set',
        'consv_set' : 'Conservative Set',
        'rescue_ratio' : 'Rescue Ratio',
        'self_consistency_ratio' : 'Self Consistency Ratio',
        'reproducibility_test' : 'Reproducibility Test',
        'FRiP' : 'Fraction of Reads in Peak',}

    short_to_long_se = {
        'total_read_pairs' : 'Total Reads',
        'distinct_read_pairs' : 'Distinct Reads',
        'one_read_pair' : 'One Read',
        'two_read_pair' : 'Two Reads',
        'NRF' : 'NRF = Distinct/Total',
        'PBC1' : 'PBC1 = OneRead/Distinct',
        'PBC2' : 'PBC2 = OneRead/TwoReads',}

    short_to_long_pe = {
        'total_read_pairs' : 'Total Read Pairs',
        'distinct_read_pairs' : 'Distinct Read Pairs',
        'one_read_pair' : 'One Read Pair',
        'two_read_pair' : 'Two Read Pairs',
        'NRF' : 'NRF = Distinct/Total',
        'PBC1' : 'PBC1 = OnePair/Distinct',
        'PBC2' : 'PBC2 = OnePair/TwoPair',}

    if key in short_to_long:
        return short_to_long[key]
    if paired_end and key in short_to_long_pe:
        return short_to_long_pe[key]
    if not paired_end and key in short_to_long_se:
        return short_to_long_se[key]
    return key

def html_heading(lvl, label):
    html = '<h{lvl}>{label}</h{lvl}>\n'
    return html.format(lvl=lvl, label=label)

def html_paragraph(txt):
    html = '<p>{txt}</p>\n'
    return html.format(txt=txt)

def html_embedded_png(png, caption, size_pct=100):
    html = '''
    <figure style="display:inline-block">
      <img src="data:image/png;base64,{encoded}" alt="{caption}" height="{size_pct}%"/>
      <figcaption style="text-align:center">{caption}</figcaption>
    </figure>
    '''
    encoded = base64.b64encode(open(png, 'rb').read()).encode('ascii')
    return html.format(size_pct=size_pct, encoded=encoded, caption=caption)

def html_vert_table_multi_rep(json_objs, paired_end=False, row_header=[]): # json_objs=list of OrderedDict
    html = '<table border="1" style="border-collapse:collapse">{header}{content}</table><br>\n'
    # make row header list
    if row_header:
        row_header = [' '] + row_header
    else:
        row_header = [' '] + ['rep'+str(i+1) 
        for i, json_obj in enumerate(json_objs)]
    # make column header list
    col_header = json_objs[0].keys()
    # make row header
    header = '<tr><th bgcolor="#EEEEEE">'+'</th><th bgcolor="#EEEEEE">'.join(row_header)+'</th></tr>\n'
    # contents
    content = ''
    for row in col_header:
        content += '<tr><th bgcolor="#EEEEEE" style="text-align:left">{}</th><td>{}</td></tr>\n'.format(
            get_long_keyname(row, paired_end),
            '</td><td>'.join(
                [str(json_obj[row]) for json_obj in json_objs]))
    return html.format(header=header, content=content)

def html_horz_table(json_obj, paired_end=False):
    html = '<table border="1" style="border-collapse:collapse">{header}{content}</table><br>\n'

    # make row header
    header = '<tr><th bgcolor="#EEEEEE">'+\
        '</th><th bgcolor="#EEEEEE">'.join(
            [get_long_keyname(key, paired_end) for key in json_obj])+'</th></tr>\n'
    # contents
    content = '<tr><td>'+\
        '</th><th>'.join([json_obj[col] for col in json_obj])+'</td></tr>\n'
    return html.format(header=header, content=content)

def html_help_pbc():
    html = """
    <div id='help-pbc'><p>NRF (non redundant fraction) <br>
    PBC1 (PCR Bottleneck coefficient 1) <br>
    PBC2 (PCR Bottleneck coefficient 2) <br>
    PBC1 is the primary measure. Provisionally <br>
    <ul>
    <li>0-0.5 is severe bottlenecking</li>
    <li>0.5-0.8 is moderate bottlenecking </li>
    <li>0.8-0.9 is mild bottlenecking </li>
    <li>0.9-1.0 is no bottlenecking </li>
    </ul></p></div><br>
    """
    return html

def html_help_xcor():
    html = """
    <div id='help-xcor'><p>
    NOTE1: For SE datasets, reads from replicates are randomly subsampled.<br>
    NOTE2: For PE datasets, the first end of each read-pair is selected and the reads are then randomly subsampled.<br>
    <ul>
    <li>Normalized strand cross-correlation coefficient (NSC) = col9 in outFile </li>
    <li>Relative strand cross-correlation coefficient (RSC) = col10 in outFile </li>
    <li>Estimated fragment length = col3 in outFile, take the top value </li>
    </ul></p></div><br>
    """
    return html

def html_help_idr(idr_thresh):
    html = """
    <div id='help-idr'><p>IDR (Irreproducible Discovery Rate) peaks<br>
    <ul>
    <li>N1: Replicate 1 self-consistent IDR {idr_thresh} peaks (comparing two pseudoreplicates generated by subsampling Rep1 reads) </li>
    <li>N2: Replicate 2 self-consistent IDR {idr_thresh} peaks (comparing two pseudoreplicates generated by subsampling Rep2 reads) </li>
    <li>Nt: True Replicate consistent IDR {idr_thresh} peaks (comparing true replicates Rep1 vs Rep2 ) </li>
    <li>Np: Pooled-pseudoreplicate consistent IDR {idr_thresh} peaks (comparing two pseudoreplicates generated by subsampling pooled reads from Rep1 and Rep2 ) </li>
    <li>Self-consistency Ratio: max(N1,N2) / min (N1,N2) </li>
    <li>Rescue Ratio: max(Np,Nt) / min (Np,Nt) </li>
    <li>Reproducibility Test: If Self-consistency Ratio >2 AND Rescue Ratio > 2, then 'Fail' else 'Pass' </li>
    </ul></p></div><br>
    """.format(idr_thresh=idr_thresh)
    return html

def html_help_idr_FRiP():
    html = """
    <div id='help-idr-FRiP'><p>
    <ul>
    <li>ppr: IDR peaks comparing pooled pseudo replicates </li>
    <li>rep1-pr: IDR peaks comparing pseudoreplicates from replicate 1 </li>
    <li>rep2-pr: IDR peaks comparing pseudoreplicates from replicate 2 </li>
    <li>repi-repj: IDR peaks comparing true replicates (rep i vs. rep j) </li>
    </ul></p></div><br>
    """
    return html

def html_help_overlap():
    html = """
    <div id='help-overlap'><p>Overlapping peaks<br>
    <ul>
    <li>N1: Replicate 1 self-consistent overlapping peaks (comparing two pseudoreplicates generated by subsampling Rep1 reads) </li>
    <li>N2: Replicate 2 self-consistent overlapping peaks (comparing two pseudoreplicates generated by subsampling Rep2 reads) </li>
    <li>Nt: True Replicate consisten overlapping peaks (comparing true replicates Rep1 vs Rep2 ) </li>
    <li>Np: Pooled-pseudoreplicate consistent overlapping peaks (comparing two pseudoreplicates generated by subsampling pooled reads from Rep1 and Rep2 ) </li>
    <li>Self-consistency Ratio: max(N1,N2) / min (N1,N2) </li>
    <li>Rescue Ratio: max(Np,Nt) / min (Np,Nt) </li>
    <li>Reproducibility Test: If Self-consistency Ratio >2 AND Rescue Ratio > 2, then 'Fail' else 'Pass' </li>
    </ul></p></div><br>
    """
    return html

def html_help_overlap_FRiP():
    html = """
    <div id='help-overlap-FRiP'><p>
    <ul>
    <li>ppr: Overlapping peaks comparing pooled pseudo replicates </li>
    <li>rep1-pr: Overlapping peaks comparing pseudoreplicates from replicate 1 </li>
    <li>rep2-pr: Overlapping peaks comparing pseudoreplicates from replicate 2 </li>
    <li>repi-repj: Overlapping peaks comparing true replicates (rep i vs. rep j) </li>
    </ul></p></div><br>
    """
    return html

def parse_bowtie2_align_log(txt):
    result = OrderedDict()
    return result

def parse_flagstat_qc(txt):
    result = OrderedDict()
    total = ''
    total_qc_failed = ''
    duplicates = ''
    duplicates_qc_failed = ''
    mapped = ''
    mapped_qc_failed = ''
    mapped_pct = ''
    paired = ''
    paired_qc_failed = ''
    read1 = ''
    read1_qc_failed = ''
    read2 = ''
    read2_qc_failed = ''
    paired_properly = ''
    paired_properly_qc_failed = ''
    paired_properly_pct = ''
    with_itself = ''
    with_itself_qc_failed = ''
    singletons = ''
    singletons_qc_failed = ''
    singletons_pct = ''
    diff_chroms = ''
    diff_chroms_qc_failed = ''

    with open(txt, 'r') as f:
        for line in f:    
            if ' in total ' in line:
                tmp1 = line.split(' in total ')
                line1 = tmp1[0]
                tmp1 = line1.split(' + ')
                total = tmp1[0]
                total_qc_failed = tmp1[1]
            if ' duplicates' in line:
                tmp2 = line.split(' duplicates')
                line2 = tmp2[0]
                tmp2 = line2.split(' + ')
                duplicates = tmp2[0]
                duplicates_qc_failed = tmp2[1]
            if ' mapped (' in line:
                tmp3 = line.split(' mapped (')               
                line3_1 = tmp3[0]
                tmp3_1 = line3_1.split(' + ')
                mapped = tmp3_1[0]
                mapped_qc_failed = tmp3_1[1]
                line3_2 = tmp3[1]
                tmp3_2 = line3_2.split(':')
                mapped_pct = tmp3_2[0].replace('%','')
            if ' paired in sequencing' in line:
                tmp2 = line.split(' paired in sequencing')
                line2 = tmp2[0]
                tmp2 = line2.split(' + ')
                paired = tmp2[0]
                paired_qc_failed = tmp2[1]
            if ' read1' in line:
                tmp2 = line.split(' read1')
                line2 = tmp2[0]
                tmp2 = line2.split(' + ')
                read1 = tmp2[0]
                read1_qc_failed = tmp2[1]
            if ' read2' in line:
                tmp2 = line.split(' read2')
                line2 = tmp2[0]
                tmp2 = line2.split(' + ')
                read2 = tmp2[0]
                read2_qc_failed = tmp2[1]
            if ' properly paired (' in line:
                tmp3 = line.split(' properly paired (')              
                line3_1 = tmp3[0]
                tmp3_1 = line3_1.split(' + ')
                paired_properly = tmp3_1[0]
                paired_properly_qc_failed = tmp3_1[1]
                line3_2 = tmp3[1]
                tmp3_2 = line3_2.split(':')
                paired_properly_pct = tmp3_2[0].replace('%','')
            if ' with itself and mate mapped' in line:
                tmp3 = line.split(' with itself and mate mapped')
                line3_1 = tmp3[0]
                tmp3_1 = line3_1.split(' + ')
                with_itself = tmp3_1[0]
                with_itself_qc_failed = tmp3_1[1]
            if ' singletons (' in line:
                tmp3 = line.split(' singletons (')               
                line3_1 = tmp3[0]
                tmp3_1 = line3_1.split(' + ')
                singletons = tmp3_1[0]
                singletons_qc_failed = tmp3_1[1]
                line3_2 = tmp3[1]
                tmp3_2 = line3_2.split(':')
                singletons_pct = tmp3_2[0].replace('%','')       
            if ' with mate mapped to a different chr' in line:
                tmp3 = line.split(' with mate mapped to a different chr')
                line3_1 = tmp3[0]
                tmp3_1 = line3_1.split(' + ')
                diff_chroms = tmp3_1[0]
                diff_chroms_qc_failed = tmp3_1[1]
    if total:
        result['total'] = int(total)
    if total_qc_failed:
        result['total_qc_failed'] = int(total_qc_failed)
    if duplicates:
        result['duplicates'] = int(duplicates)
    if duplicates_qc_failed:
        result['duplicates_qc_failed'] = int(duplicates_qc_failed)
    if mapped:
        result['mapped'] = int(mapped)
    if mapped_qc_failed:
        result['mapped_qc_failed'] = int(mapped_qc_failed)
    if mapped_pct:
        result['mapped_pct'] = float(mapped_pct)
    if paired:
        result['paired'] = int(paired)
    if paired_qc_failed:
        result['paired_qc_failed'] = int(paired_qc_failed)
    if read1:
        result['read1'] = int(read1)
    if read1_qc_failed:
        result['read1_qc_failed'] = int(read1_qc_failed)
    if read2:
        result['read2'] = int(read2)
    if read2_qc_failed:
        result['read2_qc_failed'] = int(read2_qc_failed)
    if paired_properly:
        result['paired_properly'] = int(paired_properly)
    if paired_properly_qc_failed:
        result['paired_properly_qc_failed'] = int(paired_properly_qc_failed)
    if paired_properly_pct and paired_properly_pct!='N/A':
        result['paired_properly_pct'] = float(paired_properly_pct)
    if with_itself:
        result['with_itself'] = int(with_itself)
    if with_itself_qc_failed:
        result['with_itself_qc_failed'] = int(with_itself_qc_failed)
    if singletons:
        result['singletons'] = int(singletons)
    if singletons_qc_failed:
        result['singletons_qc_failed'] = int(singletons_qc_failed)
    if singletons_pct and singletons_pct!='N/A':
        result['singletons_pct'] = float(singletons_pct)
    if diff_chroms:
        result['diff_chroms'] = int(diff_chroms)
    if diff_chroms_qc_failed:
        result['diff_chroms_qc_failed'] = int(diff_chroms_qc_failed)
    return result

def parse_dup_qc(txt):
    result = OrderedDict()
    paired_reads = ''
    unpaired_reads = ''
    unmapped_reads = ''
    unpaired_dupes = ''
    paired_dupes = ''
    paired_opt_dupes = ''
    dupes_pct = ''

    picard_log_found = False
    # picard markdup
    with open(txt, 'r') as f:
        header = '' # if 'UNPAIRED_READS_EXAMINED' in header
        content = ''
        for line in f:
            if header:
                content = line
                picard_log_found = True
                break
            if 'UNPAIRED_READS_EXAMINED' in line:
                header = line
    if picard_log_found:
        header_items = header.split('\t')
        content_items = content.split('\t')
        m = dict(zip(header_items, content_items))
        unpaired_reads = m['UNPAIRED_READS_EXAMINED']
        paired_reads = m['READ_PAIRS_EXAMINED']
        unmapped_reads = m['UNMAPPED_READS']
        unpaired_dupes = m['UNPAIRED_READ_DUPLICATES']
        paired_dupes = m['READ_PAIR_DUPLICATES']
        paired_opt_dupes = m['READ_PAIR_OPTICAL_DUPLICATES']
        dupes_pct = m['PERCENT_DUPLICATION']
    else:
        # sambamba markdup
        with open(txt, 'r') as f:
            for line in f:    
                if ' end pairs' in line:
                    tmp1 = line.strip().split(' ')
                    paired_reads = tmp1[1]
                if ' single ends ' in line:
                    tmp1 = line.strip().split(' ')
                    unpaired_reads = tmp1[1]
                    unmapped_reads = tmp1[6]
                if 'found ' in line:
                    tmp1 = line.strip().split(' ')
                    if paired_reads == '0':
                        unpaired_dupes = tmp1[1] # SE
                        paired_dupes = 0
                    else:
                        unpaired_dupes = 0
                        paired_dupes = str(int(tmp1[1])/2) # PE
                if paired_reads == '0': # SE
                    dupes_pct = '{0:.2f}'.format(
                                float(unpaired_dupes)/float(unpaired_reads))
                else:
                    dupes_pct = '{0:.2f}'.format(
                                float(paired_dupes)/float(paired_reads))
    if unpaired_reads:
        result['unpaired_reads'] = int(unpaired_reads)
    if paired_reads:
        result['paired_reads'] = int(paired_reads)
    if unmapped_reads:
        result['unmapped_reads'] = int(unmapped_reads)
    if unpaired_dupes:
        result['unpaired_dupes'] = int(unpaired_dupes)
    if paired_dupes:
        result['paired_dupes'] = int(paired_dupes)
    if paired_opt_dupes:
        result['paired_opt_dupes'] = int(paired_opt_dupes)
    if dupes_pct:
        result['dupes_pct'] = float(dupes_pct)
    return result

def parse_pbc_qc(txt):
    result = OrderedDict()
    with open(txt, 'r') as f:
        for line in f:
            arr = line.strip().split('\t')
            break
    result['total_read_pairs'] = int(arr[0])
    result['distinct_read_pairs'] = int(arr[1])
    result['one_read_pair'] = int(arr[2])
    result['two_read_pair'] = int(arr[3])
    result['NRF'] = float(arr[4])
    result['PBC1'] = float(arr[5])
    result['PBC2'] = float(arr[6])
    return result

def parse_xcor_score(txt):
    result = OrderedDict()
    with open(txt, 'r') as f:
        arr = f.readlines()[0].strip().split('\t')
    result['num_reads'] = int(arr[1])
    result['est_frag_len'] = int(arr[2])
    result['corr_est_frag_len'] = float(arr[3])
    result['phantom_peak'] = int(arr[4])
    result['corr_phantom_peak'] = float(arr[5])
    result['argmin_corr'] = int(arr[6])
    result['min_corr'] = float(arr[7])
    result['NSC'] = float(arr[8])
    result['RSC'] = float(arr[9])
    return result

def parse_reproducibility_qc(txt):
    with open(txt, 'r') as f:
        lines = f.readlines()
        header = lines[0].strip()
        content = lines[1].strip()
        result = OrderedDict(
            zip(header.split('\t'), content.split('\t')))
        for key in result:
            if key.startswith('N'):
                result[key] = int(result[key])
            if key.endswith('_ratio'):
                result[key] = float(result[key])
    return result

def parse_frip_qc(txt):
    result = OrderedDict()
    with open(txt, 'r') as f:
        frip = f.readlines()[0].strip()
    result['FRiP'] = float(frip)
    return result

def parse_multi_col_txt(txt): # to read ATAQC log
    result = OrderedDict()
    with open(txt, 'r') as f:
        lines = f.readlines()    
    for line in lines:
        line = line.strip().replace(' reads; of these:','')
        arr = line.split('\t')
        for j in range(1,len(arr)):
            header = arr[0] 
            header += '' if len(arr)==2 else '-{}'.format(j)
            content = arr[j]
            result[header] = content
    return result

def main():
    # read params
    args = parse_arguments()

    log.info('Initializing...')
    html = ''
    json_all = OrderedDict()

    if args.name:
        html += html_heading(1, args.name)
    if args.desc:
        html += html_paragraph(args.desc)
    html += html_paragraph('Report generated at {}'.format(now()))

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
                [infer_pair_label_from_idx(num_rep, i)
                    for i in range(len(args.frip_overlap_qcs))])
        if args.frip_overlap_qcs_pr:
            json_objs = [parse_frip_qc(qc) for qc in args.frip_overlap_qcs_pr]
            json_objs_frip_overlap.extend(json_objs)
            row_header_frip_overlap.extend(
                ['rep{}-pr'.format(i+1) 
                    for i in range(len(args.frip_overlap_qcs_pr))])
        if args.frip_overlap_qc_ppr:
            json_obj = parse_frip_qc(args.frip_overlap_qc_ppr[0])
            json_objs_frip_overlap.append(json_obj)
            row_header_frip_overlap.append('ppr')
        json_all['overlap_frip_qc'] = OrderedDict(
            zip(row_header_frip_overlap,json_objs_frip_overlap))
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
                [infer_pair_label_from_idx(num_rep, i)
                    for i in range(len(args.frip_idr_qcs))])
        if args.frip_idr_qcs_pr:
            json_objs = [parse_frip_qc(qc) for qc in args.frip_idr_qcs_pr]
            json_objs_frip_idr.extend(json_objs)
            row_header_frip_idr.extend(
                ['rep{}-pr'.format(i+1) 
                    for i in range(len(args.frip_idr_qcs_pr))])
        if args.frip_idr_qc_ppr:
            json_obj = parse_frip_qc(args.frip_idr_qc_ppr[0])
            json_objs_frip_idr.append(json_obj)
            row_header_frip_idr.append('ppr')
        json_all['idr_frip_qc'] = OrderedDict(
            zip(row_header_frip_idr,json_objs_frip_idr))
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

    if json_all:
        log.info('Write JSON file...')        
        write_txt(args.out_qc_json, json.dumps(json_all, indent=4))
        # b = json.loads(a, object_pairs_hook=OrderedDict)

    log.info('All done.')

if __name__=='__main__':
    main()
