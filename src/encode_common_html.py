#!/usr/bin/env python

# ENCODE DCC reporting module wrapper
# Author: Jin Lee (leepc12@gmail.com)

import base64
import re
from encode_common_log_parser import get_long_keyname
from encode_common import *
from collections import OrderedDict

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
    encoded = base64.b64encode(open(png, 'rb').read()).decode("utf-8")
    return html.format(size_pct=size_pct, encoded=encoded, caption=caption)

def html_parse_body_from_file(f):
    with open(f,'r') as fp:
        # p = re.compile('<body[^>]*>((.|[\n\r])*)<\/body>')
        # return p.search(fp.read()).group()
        p = re.compile('<body[^>]*>((.|[\n\r])*)<\/body>')
        s = p.search(fp.read())
        if s:
            return s.group()
        else:
            return ''

def html_vert_table_multi_rep(json_objs, row_header=[]): # json_objs=list of OrderedDict
    html = '<table border="1" style="border-collapse:collapse">{header}{content}</table><br>\n'

    valid_row_header = []
    valid_json_objs = []
    col_header_dict = OrderedDict()
    for i, json_obj in enumerate(json_objs):
        if not json_obj:
            continue
        if row_header:
            valid_row_header.append(row_header[i])
        else:
            valid_row_header.append('rep'+str(i+1))
        valid_json_objs.append(json_obj)
        # find column header
        col_header_dict.update(json_obj)

    col_header = col_header_dict.keys()
    if not col_header:
        return ''
    # make row header
    header = '<tr><th bgcolor="#EEEEEE">'+'</th><th bgcolor="#EEEEEE">'.join([' ']+valid_row_header)+'</th></tr>\n'
    # contents
    content = ''
    for row in col_header:
        content += '<tr><th bgcolor="#EEEEEE" style="text-align:left">{}</th><td>{}</td></tr>\n'.format(
            get_long_keyname(row),
            '</td><td>'.join(
                [str_float_4_dec_pts(json_obj[row]) if row in json_obj else 'N/A'
                    for json_obj in valid_json_objs]))
    return html.format(header=header, content=content)

def html_help_filter(multimapping):
    html = """
    <div id='help-filter'>
    Filtered out (samtools view -F 1804):
    <ul>
    <li>read unmapped (0x4)</li>
    <li>mate unmapped (0x8, for paired-end)</li>"
    <li>not primary alignment (0x100)</li>
    <li>read fails platform/vendor quality checks (0x200)</li>
    <li>read is PCR or optical duplicate (0x400)</li>
    </ul></p></div><br>
    """
    return html

def html_help_pbc():
    html = """
    <div id='help-pbc'>
    Mitochondrial reads are filtered out.
    <p>NRF (non redundant fraction) <br>
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

def html_help_FRiP(peak_caller):
    html = """
    <div id='help-FRiP'><p>
    <ul>
    <li>repX: {peak_caller} peaks from true replicate X </li>
    <li>repX-prY: {peak_caller} peaks from Yth pseudoreplicates from replicate X </li>
    <li>pooled: {peak_caller} peaks from pooled true replicates </li>
    <li>ppr1: {peak_caller} peaks from 1st pooled pseudo replicates </li>
    <li>ppr2: {peak_caller} peaks from 2nd pooled pseudo replicates </li>
    </ul></p></div><br>
    """.format(peak_caller=peak_caller.upper())
    return html

def html_help_macs2(cap_num_peak):
    html = """
    <div id='help-macs2'><p>
    The number of peaks is capped at {} for peak-caller MACS2
    </p></div><br>
    """.format(human_readable_number(cap_num_peak))
    return html if cap_num_peak else ""

def html_help_spp(cap_num_peak):
    html = """
    <div id='help-spp'><p>
    The number of peaks is capped at {} for peak-caller SPP
    </p></div><br>
    """.format(human_readable_number(cap_num_peak))
    return html if cap_num_peak else ""

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

# for float, limit it to 4 decimal points
# join if list is given
def str_float_4_dec_pts(num):
    if type(num)==list:
        return ', '.join([str(n) for n in num])
    if type(num)==int:
        return str(num)
    elif type(num)==float:
        return '{0:.4f}'.format(num)
    else:
        return str(num)

