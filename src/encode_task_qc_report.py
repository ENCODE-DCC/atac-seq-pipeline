#!/usr/bin/env python

# ENCODE reporting module wrapper
# Author: Jin Lee (leepc12@gmail.com)

import sys
import json
import argparse
from collections import OrderedDict
from encode_lib_common import (
    infer_n_from_nC2, infer_pair_label_from_idx, log, now,
    run_shell_cmd, write_txt)
from encode_lib_log_parser import (
    MAP_KEY_DESC_ANNOT_ENRICH_QC, MAP_KEY_DESC_DUP_QC, MAP_KEY_DESC_FLAGSTAT_QC,
    MAP_KEY_DESC_FRAC_MITO_QC, MAP_KEY_DESC_FRIP_QC, MAP_KEY_DESC_JSD_QC,
    MAP_KEY_DESC_LIB_COMPLEXITY_QC, MAP_KEY_DESC_NUCLEOSOMAL_QC,
    MAP_KEY_DESC_NUM_PEAK_QC, MAP_KEY_DESC_PEAK_REGION_SIZE_QC,
    MAP_KEY_DESC_PICARD_EST_LIB_SIZE_QC, MAP_KEY_DESC_REPRODUCIBILITY_QC,
    MAP_KEY_DESC_TSS_ENRICH_QC, MAP_KEY_DESC_XCOR_SCORE)
from encode_lib_log_parser import (
    parse_annot_enrich_qc, parse_dup_qc, parse_flagstat_qc, parse_frac_mito_qc,
    parse_frip_qc, parse_jsd_qc, parse_lib_complexity_qc, parse_nucleosomal_qc,
    parse_num_peak_qc, parse_peak_region_size_qc, parse_picard_est_lib_size_qc,
    parse_reproducibility_qc, parse_tss_enrich_qc, parse_xcor_score)
from encode_lib_qc_category import QCCategory


def parse_arguments():
    parser = argparse.ArgumentParser(
        prog='ENCODE Final QC report/JSON generator.')
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
    parser.add_argument('--paired-end', action='store_true',
                        help='Endedness of all replicates.')
    parser.add_argument('--ctl-paired-end', action='store_true',
                        help='Endedness of all controls.')
    parser.add_argument('--paired-ends', type=str, nargs='*',
                        help='List of true/false for paired endedness of sample.')
    parser.add_argument('--ctl-paired-ends', type=str, nargs='*',
                        help='List of true/false for paired endedness of control.')
    parser.add_argument('--pipeline-type', type=str, required=True,
                        help='Pipeline type.')
    parser.add_argument('--aligner', type=str, required=True,
                        help='Aligner.')
    parser.add_argument('--peak-caller', type=str, required=True,
                        help='Peak caller.')
    parser.add_argument('--cap-num-peak', default=0, type=int,
                        help='Capping number of peaks by taking top N peaks.')
    parser.add_argument('--idr-thresh', type=float, required=True,
                        help='IDR threshold.')
    parser.add_argument('--pval-thresh', type=float,
                        help='pValue threshold for MACS2 peak caller.')
    parser.add_argument('--xcor-trim-bp', type=int,
                        help='FASTQs are trimmed to this for cross-correlation '
                             'analysis only.')
    parser.add_argument('--xcor-subsample-reads', type=int,
                        help='Subsampled TAG-ALIGNs to this depth for cross-correlation '
                             'analysis only.')
    parser.add_argument('--samstat-qcs', type=str, nargs='*',
                        help='List of samstat QC (raw BAM) files per replicate.')
    parser.add_argument('--nodup-samstat-qcs', type=str, nargs='*',
                        help='List of samstat QC (filtered BAM) files per replicate.')
    parser.add_argument('--dup-qcs', type=str, nargs='*',
                        help='List of dup QC files per replicate.')
    parser.add_argument('--lib-complexity-qcs', type=str, nargs='*',
                        help='List of lib complexity QC files per replicate.')
    parser.add_argument('--frac-mito-qcs', type=str, nargs='*',
                        help='List of fraction of mito QC files per replicate.')
    parser.add_argument('--ctl-samstat-qcs', type=str, nargs='*',
                        help='List of samstat QC (raw BAM) files per control.')
    parser.add_argument('--ctl-nodup-samstat-qcs', type=str, nargs='*',
                        help='List of samstat QC (filtered BAM) files per control.')
    parser.add_argument('--ctl-dup-qcs', type=str, nargs='*',
                        help='List of dup QC files per control.')
    parser.add_argument('--ctl-lib-complexity-qcs', type=str, nargs='*',
                        help='List of lib complexity QC files per control.')
    parser.add_argument('--ctl-frac-mito-qcs', type=str, nargs='*',
                        help='List of fraction of mito QC files per control.')
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

    parser.add_argument('--frip-qcs', type=str, nargs='*',
                        help='List of raw peak FRiP score files per replicate.')
    parser.add_argument('--frip-qcs-pr1', type=str, nargs='*',
                        help='List of raw peak FRiP score files for 1st pseudo replicates per replicate.')
    parser.add_argument('--frip-qcs-pr2', type=str, nargs='*',
                        help='List of raw peak FRiP score files for 2nd pseudo replicates per replicate.')
    parser.add_argument('--frip-qc-pooled', type=str, nargs='*',
                        help='Raw peak FRiP score file for pooled replicates.')
    parser.add_argument('--frip-qc-ppr1', type=str, nargs='*',
                        help='Raw peak FRiP score file for 1st pooled pseudo replicates.')
    parser.add_argument('--frip-qc-ppr2', type=str, nargs='*',
                        help='Raw peak FRiP score file for 2nd pooled pseudo replicates.')

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

    parser.add_argument('--annot-enrich-qcs', type=str, nargs='*',
                        help='List of annot_enrich QC files.')
    parser.add_argument('--tss-enrich-qcs', type=str, nargs='*',
                        help='List of TSS enrichment QC files.')
    parser.add_argument('--tss-large-plots', type=str, nargs='*',
                        help='List of TSS enrichment large plots.')
    parser.add_argument('--roadmap-compare-plots', type=str, nargs='*',
                        help='List of Roadmap-compare plots.')
    parser.add_argument('--fraglen-dist-plots', type=str, nargs='*',
                        help='List of fragment length distribution plots.')
    parser.add_argument('--fraglen-nucleosomal-qcs', type=str, nargs='*',
                        help='List of fragment length nucleosomal QC files.')
    parser.add_argument('--gc-plots', type=str, nargs='*',
                        help='List of GC bias plots.')
    parser.add_argument('--preseq-plots', type=str, nargs='*',
                        help='List of preseq plots.')
    parser.add_argument('--picard-est-lib-size-qcs', type=str, nargs='*',
                        help='List of Picard estimated library size QC files.')

    parser.add_argument('--peak-region-size-qcs', type=str, nargs='*',
                        help='List of peak region size QC files.')
    parser.add_argument('--peak-region-size-plots', type=str, nargs='*',
                        help='List of peak region size plot files.')
    parser.add_argument('--num-peak-qcs', type=str, nargs='*',
                        help='List of QC files with number of peaks.')

    parser.add_argument('--idr-opt-peak-region-size-qc', type=str, nargs='*',
                        help='IDR opt peak region size QC files.')
    parser.add_argument('--idr-opt-peak-region-size-plot', type=str, nargs='*',
                        help='IDR opt peak region size plot files.')
    parser.add_argument('--idr-opt-num-peak-qc', type=str, nargs='*',
                        help='QC file with number of peaks in IDR opt peak.')

    parser.add_argument('--overlap-opt-peak-region-size-qc', type=str, nargs='*',
                        help='overlap opt peak region size QC files.')
    parser.add_argument('--overlap-opt-peak-region-size-plot', type=str, nargs='*',
                        help='overlap opt peak region size plot files.')
    parser.add_argument('--overlap-opt-num-peak-qc', type=str, nargs='*',
                        help='QC file with number of peaks in overlap opt peak.')

    parser.add_argument('--out-qc-html', default='qc.html', type=str,
                        help='Output QC report HTML file.')
    parser.add_argument('--out-qc-json', default='qc.json', type=str,
                        help='Output QC JSON file.')
    parser.add_argument('--qc-json-ref', type=str,
                        help='Reference QC JSON file to be compared to output QC JSON (developer\'s purpose only).')
    parser.add_argument('--log-level', default='INFO',
                        choices=['NOTSET', 'DEBUG', 'INFO',
                                 'WARNING', 'CRITICAL', 'ERROR', 'CRITICAL'],
                        help='Log level')
    args = parser.parse_args()

    # parse with a special delimiter "_:_"
    for a in vars(args):
        value = getattr(args, a)
        if isinstance(value, list):
            setattr(args, a, split_entries_and_extend(value))

    if args.paired_ends is not None:
        for i, _ in enumerate(args.paired_ends):
            args.paired_ends[i] = str2bool(args.paired_ends[i])

    if args.ctl_paired_ends is not None:
        for i, _ in enumerate(args.ctl_paired_ends):
            args.ctl_paired_ends[i] = str2bool(args.ctl_paired_ends[i])

    log.setLevel(args.log_level)
    log.info(sys.argv)
    return args


def split_entries_and_extend(l, delim='_:_'):
    result = []
    for a in l:
        if isinstance(a, str):
            result.extend(a.split(delim))
        else:
            result.append(a)
    test_not_all_empty = [r for r in result if r]
    # if all empty then return []
    return result if test_not_all_empty else []


def str2bool(s):
    s = s.lower()
    if s not in ('false', 'true'):
        raise ValueError('Not a valid boolean string')
    return s == 'true'


def str_rep(i):
    return 'rep' + str(i + 1)


def str_ctl(i):
    return 'ctl' + str(i + 1)


MAP_KEY_DESC_GENERAL = {
    'date': 'Report generated at',
    'title': 'Title',
    'description': 'Description',
    'pipeline_ver': 'Pipeline version',
    'pipeline_type': 'Pipeline type',
    'aligner': 'Aligner',
    'peak_caller': 'Peak caller',
    'genome': 'Genome',
    'seq_endedness': 'Sequencing endedness'
}


def make_cat_root(args):
    cat_root = QCCategory(
        'root',
        html_head='<h1>QC Report</h1><hr>',
        map_key_desc=MAP_KEY_DESC_GENERAL
    )

    # add general information dict to cat_root
    d_general = dict([
        ('date', now()),
        ('title', args.title),
        ('description', args.desc),
        ('pipeline_ver', args.pipeline_ver),
        ('pipeline_type', args.pipeline_type),
        ('genome', args.genome),
        ('aligner', args.aligner),
        ('seq_endedness', {}),
        ('peak_caller', args.peak_caller),
    ])
    if args.paired_ends is not None:
        for i, paired_end in enumerate(args.paired_ends):
            d_general['seq_endedness']['rep{}'.format(i + 1)] = {'paired_end': paired_end}
    if args.ctl_paired_ends is not None:
        for i, paired_end in enumerate(args.ctl_paired_ends):
            d_general['seq_endedness']['ctl{}'.format(i + 1)] = {'paired_end': paired_end}
    cat_root.add_log(d_general, key='general')

    return cat_root


def make_cat_align(args, cat_root):
    cat_align = QCCategory(
        'align',
        html_head='<h1>Alignment quality metrics</h1><hr>',
        parent=cat_root
    )

    cat_align_samstat = QCCategory(
        'samstat',
        html_head='<h2>SAMstat (raw unfiltered BAM)</h2>',
        parser=parse_flagstat_qc,
        map_key_desc=MAP_KEY_DESC_FLAGSTAT_QC,
        parent=cat_align
    )
    if args.samstat_qcs:
        for i, qc in enumerate(args.samstat_qcs):
            if qc:
                cat_align_samstat.add_log(qc, key=str_rep(i))
    if args.ctl_samstat_qcs:
        for i, qc in enumerate(args.ctl_samstat_qcs):
            if qc:
                cat_align_samstat.add_log(qc, key=str_ctl(i))

    cat_align_dup = QCCategory(
        'dup',
        html_head='<h2>Marking duplicates (filtered BAM)</h2>',
        html_foot="""
            <div id='help-filter'>
            Filtered out (samtools view -F 1804):
            <ul>
            <li>read unmapped (0x4)</li>
            <li>mate unmapped (0x8, for paired-end)</li>
            <li>not primary alignment (0x100)</li>
            <li>read fails platform/vendor quality checks (0x200)</li>
            <li>read is PCR or optical duplicate (0x400)</li>
            </ul></p></div><br>
        """,
        parser=parse_dup_qc,
        map_key_desc=MAP_KEY_DESC_DUP_QC,
        parent=cat_align
    )
    if args.dup_qcs:
        for i, qc in enumerate(args.dup_qcs):
            if qc:
                cat_align_dup.add_log(qc, key=str_rep(i))
    if args.ctl_dup_qcs:
        for i, qc in enumerate(args.ctl_dup_qcs):
            if qc:
                cat_align_dup.add_log(qc, key=str_ctl(i))

    cat_align_frac_mito = QCCategory(
        'frac_mito',
        html_head='<h2>Fraction of mitochondrial reads (unfiltered BAM)</h2>',
        parser=parse_frac_mito_qc,
        map_key_desc=MAP_KEY_DESC_FRAC_MITO_QC,
        parent=cat_align
    )
    if args.frac_mito_qcs:
        for i, qc in enumerate(args.frac_mito_qcs):
            if qc:
                cat_align_frac_mito.add_log(qc, key=str_rep(i))
    if args.ctl_frac_mito_qcs:
        for i, qc in enumerate(args.ctl_frac_mito_qcs):
            if qc:
                cat_align_frac_mito.add_log(qc, key=str_ctl(i))

    cat_align_preseq = QCCategory(
        'preseq',
        html_foot="""
            <p>Preseq performs a yield prediction by subsampling the reads, calculating the
            number of distinct reads, and then extrapolating out to see where the
            expected number of distinct reads no longer increases. The confidence interval
            gives a gauge as to the validity of the yield predictions.</p>
        """,
        parent=cat_align
    )
    if args.preseq_plots:
        for i, plot in enumerate(args.preseq_plots):
            if plot:
                cat_align_preseq.add_plot(
                    plot, key=str_rep(i), size_pct=50)

    cat_align_nodup_samstat = QCCategory(
        'nodup_samstat',
        html_head='<h2>SAMstat (filtered/deduped BAM)</h2>',
        html_foot="""
            <p>Filtered and duplicates removed</p><br>
        """,
        parser=parse_flagstat_qc,
        map_key_desc=MAP_KEY_DESC_FLAGSTAT_QC,
        parent=cat_align
    )
    if args.nodup_samstat_qcs:
        for i, qc in enumerate(args.nodup_samstat_qcs):
            if qc:
                cat_align_nodup_samstat.add_log(qc, key=str_rep(i))
    if args.ctl_nodup_samstat_qcs:
        for i, qc in enumerate(args.ctl_nodup_samstat_qcs):
            if qc:
                cat_align_nodup_samstat.add_log(qc, key=str_ctl(i))

    cat_fraglen = QCCategory(
        'frag_len_stat',
        html_head='<h2>Fragment length statistics (filtered/deduped BAM)</h2>',
        html_foot="""
            <p>Open chromatin assays show distinct fragment length enrichments, as the cut
            sites are only in open chromatin and not in nucleosomes. As such, peaks
            representing different n-nucleosomal (ex mono-nucleosomal, di-nucleosomal)
            fragment lengths will arise. Good libraries will show these peaks in a
            fragment length distribution and will show specific peak ratios.</p><br>
            <ul>
            <li>NFR: Nucleosome free region
            </ul><br>
        """,
        parser=parse_nucleosomal_qc,
        map_key_desc=MAP_KEY_DESC_NUCLEOSOMAL_QC,
        parent=cat_align
    )

    if args.fraglen_nucleosomal_qcs:
        for i, qc in enumerate(args.fraglen_nucleosomal_qcs):
            if qc:
                cat_fraglen.add_log(qc, key=str_rep(i))
    if args.fraglen_dist_plots:
        for i, plot in enumerate(args.fraglen_dist_plots):
            if plot:
                cat_fraglen.add_plot(plot,
                                     key=str_rep(i), size_pct=50)

    cat_gc_bias = QCCategory(
        'gc_bias',
        html_head='<h2>Sequence quality metrics (filtered/deduped BAM)</h2>',
        html_foot="""
            <p>Open chromatin assays are known to have significant GC bias. Please take this
            into consideration as necessary.</p><br>
        """,
        parent=cat_align
    )
    if args.gc_plots:
        for i, plot in enumerate(args.gc_plots):
            if plot:
                cat_gc_bias.add_plot(plot,
                                     key=str_rep(i), size_pct=60)

    return cat_align


def make_cat_lib_complexity(args, cat_root):
    cat_lc = QCCategory(
        'lib_complexity',
        html_head='<h1>Library complexity quality metrics</h1><hr>',
        parent=cat_root
    )

    html_head_lib_complexity = """
        <div id='help-lib_complexity'>
        <p>Mitochondrial reads are filtered out by default.
        The non-redundant fraction (NRF) is the fraction of non-redundant mapped reads
        in a dataset; it is the ratio between the number of positions in the genome
        that uniquely mapped reads map to and the total number of uniquely mappable
        reads. The NRF should be > 0.8. The PBC1 is the ratio of genomic locations
        with EXACTLY one read pair over the genomic locations with AT LEAST one read
        pair. PBC1 is the primary measure, and the PBC1 should be close to 1.
        Provisionally 0-0.5 is severe bottlenecking, 0.5-0.8 is moderate bottlenecking,
        0.8-0.9 is mild bottlenecking, and 0.9-1.0 is no bottlenecking. The PBC2 is
        the ratio of genomic locations with EXACTLY one read pair over the genomic
        locations with EXACTLY two read pairs. The PBC2 should be significantly
        greater than 1. {pipeline_specific_info}
        </p><br>
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
    """.format(
        pipeline_specific_info="""See more details at
        <a href="https://www.encodeproject.org/chip-seq/transcription_factor/#standards">
        the ENCODE portal standard for ChIP-Seq pipeline</a>
        """ if args.pipeline_type in ('tf', 'histone') else ''
        )

    cat_lc_lib_complexity = QCCategory(
        'lib_complexity',
        html_head='<h2>Library complexity (filtered non-mito BAM)</h2>',
        html_foot=html_head_lib_complexity,
        parser=parse_lib_complexity_qc,
        map_key_desc=MAP_KEY_DESC_LIB_COMPLEXITY_QC,
        parent=cat_lc
    )
    if args.lib_complexity_qcs:
        for i, qc in enumerate(args.lib_complexity_qcs):
            if qc:
                cat_lc_lib_complexity.add_log(qc, key=str_rep(i))
    if args.ctl_lib_complexity_qcs:
        for i, qc in enumerate(args.ctl_lib_complexity_qcs):
            if qc:
                cat_lc_lib_complexity.add_log(qc, key=str_ctl(i))

    cat_lc_lib_size = QCCategory(
        'lib_size',
        parser=parse_picard_est_lib_size_qc,
        map_key_desc=MAP_KEY_DESC_PICARD_EST_LIB_SIZE_QC,
        parent=cat_lc
    )
    if args.picard_est_lib_size_qcs:
        for i, qc in enumerate(args.picard_est_lib_size_qcs):
            if qc:
                cat_lc_lib_size.add_log(qc, key=str_rep(i))

    return cat_lc


def make_cat_replication(args, cat_root):
    cat_replication = QCCategory(
        'replication',
        html_head='<h1>Replication quality metrics</h1><hr>',
        parent=cat_root
    )

    cat_idr = QCCategory(
        'idr',
        html_head='<h2>IDR (Irreproducible Discovery Rate) plots</h2>',
        parent=cat_replication,
    )
    if args.idr_plots:
        num_rep = infer_n_from_nC2(len(args.idr_plots))
        for i, plot in enumerate(args.idr_plots):
            if plot:
                cat_idr.add_plot(
                    plot,
                    key=infer_pair_label_from_idx(num_rep, i))
    if args.idr_plots_pr:
        for i, plot in enumerate(args.idr_plots_pr):
            if plot:
                cat_idr.add_plot(
                    plot,
                    key='rep{X}-pr1_vs_rep{X}-pr2'.format(X=i+1))
    if args.idr_plot_ppr:
        cat_idr.add_plot(args.idr_plot_ppr[0], key='pooled-pr1_vs_pooled-pr2')

    cat_reproducibility = QCCategory(
        'reproducibility',
        html_head='<h2>Reproducibility QC and peak detection statistics</h2>',
        html_foot="""
            <div id='help-reproducibility'><p>Reproducibility QC<br>
            <ul>
            <li>N1: Replicate 1 self-consistent peaks (comparing two pseudoreplicates generated by subsampling Rep1 reads) </li>
            <li>N2: Replicate 2 self-consistent peaks (comparing two pseudoreplicates generated by subsampling Rep2 reads) </li>
            <li>Ni: Replicate i self-consistent peaks (comparing two pseudoreplicates generated by subsampling RepX reads) </li>
            <li>Nt: True Replicate consistent peaks (comparing true replicates Rep1 vs Rep2) </li>
            <li>Np: Pooled-pseudoreplicate consistent peaks (comparing two pseudoreplicates generated by subsampling pooled reads from Rep1 and Rep2) </li>
            <li>Self-consistency Ratio: max(N1,N2) / min (N1,N2) </li>
            <li>Rescue Ratio: max(Np,Nt) / min (Np,Nt) </li>
            <li>Reproducibility Test: If Self-consistency Ratio >2 AND Rescue Ratio > 2, then 'Fail' else 'Pass' </li>
            </ul></p></div><br>
        """,
        parser=parse_reproducibility_qc,
        map_key_desc=MAP_KEY_DESC_REPRODUCIBILITY_QC,
        parent=cat_replication,
    )
    if args.overlap_reproducibility_qc:
        qc = args.overlap_reproducibility_qc[0]
        cat_reproducibility.add_log(qc, key='overlap')

    if args.idr_reproducibility_qc:
        qc = args.idr_reproducibility_qc[0]
        cat_reproducibility.add_log(qc, key='idr')

    if args.peak_caller == 'spp':
        extra_info = 'with FDR 0.01'
    elif args.peak_caller == 'macs2':
        extra_info = 'with p-val threshold {}'.format(args.pval_thresh)
    else:
        extra_info = ''

    cat_num_peak = QCCategory(
        'num_peaks',
        html_head='<h2>Number of raw peaks</h2>',
        html_foot="""
            Top {num_peak} raw peaks from {peak_caller} {extra_info}
        """.format(
            num_peak=args.cap_num_peak,
            peak_caller=args.peak_caller,
            extra_info=extra_info,
        ),
        parser=parse_num_peak_qc,
        map_key_desc=MAP_KEY_DESC_NUM_PEAK_QC,
        parent=cat_replication,
    )
    if args.num_peak_qcs:
        for i, qc in enumerate(args.num_peak_qcs):
            if qc:
                cat_num_peak.add_log(qc, key=str_rep(i))

    return cat_replication


def make_cat_peak_stat(args, cat_root):
    cat_peak_stat = QCCategory(
        'peak_stat',
        html_head='<h1>Peak calling statistics</h1><hr>',
        parent=cat_root
    )

    cat_peak_region_size = QCCategory(
        'peak_region_size',
        html_head='<h2>Peak region size</h2>',
        html_foot="""
        """,
        parser=parse_peak_region_size_qc,
        map_key_desc=MAP_KEY_DESC_PEAK_REGION_SIZE_QC,
        parent=cat_peak_stat,
    )
    if args.peak_region_size_qcs:
        for i, qc in enumerate(args.peak_region_size_qcs):
            if qc:
                cat_peak_region_size.add_log(qc, key=str_rep(i))
    if args.idr_opt_peak_region_size_qc:
        qc = args.idr_opt_peak_region_size_qc[0]
        cat_peak_region_size.add_log(qc, key='idr_opt')
    if args.overlap_opt_peak_region_size_qc:
        qc = args.overlap_opt_peak_region_size_qc[0]
        cat_peak_region_size.add_log(qc, key='overlap_opt')
    # plots
    if args.peak_region_size_plots:
        for i, plot in enumerate(args.peak_region_size_plots):
            if plot:
                cat_peak_region_size.add_plot(plot, key=str_rep(i),
                                              size_pct=35)
    if args.idr_opt_peak_region_size_plot:
        plot = args.idr_opt_peak_region_size_plot[0]
        cat_peak_region_size.add_plot(plot, key='idr_opt', size_pct=35)
    if args.overlap_opt_peak_region_size_plot:
        plot = args.overlap_opt_peak_region_size_plot[0]
        cat_peak_region_size.add_plot(plot, key='overlap_opt', size_pct=35)

    return cat_peak_stat


def make_cat_align_enrich(args, cat_root):
    cat_align_enrich = QCCategory(
        'align_enrich',
        html_head='<h1>Enrichment / Signal-to-noise ratio</h1><hr>',
        parent=cat_root
    )

    if args.pipeline_type in ('tf', 'histone'):
        html_head_xcor = '<h2>Strand cross-correlation measures (trimmed/filtered SE BAM)</h2>'
        html_foot_xcor = """
            <br><p>Performed on subsampled ({xcor_subsample_reads}) reads mapped from FASTQs that are trimmed to {xcor_trim_bp}.
            Such FASTQ trimming and subsampling reads are for cross-corrleation analysis only. 
            Untrimmed FASTQs are used for all the other analyses.</p>
            <div id='help-xcor'><p>
            NOTE1: For SE datasets, reads from replicates are randomly subsampled to {xcor_subsample_reads}.<br>
            NOTE2: For PE datasets, the first end (R1) of each read-pair is selected and trimmed to {xcor_trim_bp} the reads are then randomly subsampled to {xcor_subsample_reads}.<br>
        """.format(
            xcor_subsample_reads=args.xcor_subsample_reads,
            xcor_trim_bp=args.xcor_trim_bp,
        )
    else:
        html_head_xcor = '<h2>Strand cross-correlation measures (filtered BAM)</h2>'
        html_foot_xcor = """
            <br><p>Performed on subsampled ({xcor_subsample_reads}) reads.
            Such FASTQ trimming is for cross-corrleation analysis only.</p>
            <div id='help-xcor'><p>
        """.format(
            xcor_subsample_reads=args.xcor_subsample_reads
        )
    html_foot_xcor += """<ul>
        <li>Normalized strand cross-correlation coefficient (NSC) = col9 in outFile </li>
        <li>Relative strand cross-correlation coefficient (RSC) = col10 in outFile </li>
        <li>Estimated fragment length = col3 in outFile, take the top value </li>
        </ul></p></div><br>
    """

    cat_xcor = QCCategory(
        'xcor_score',
        html_head=html_head_xcor,
        html_foot=html_foot_xcor,
        parser=parse_xcor_score,
        map_key_desc=MAP_KEY_DESC_XCOR_SCORE,
        parent=cat_align_enrich,
    )
    if args.xcor_scores:
        for i, qc in enumerate(args.xcor_scores):
            if qc:
                cat_xcor.add_log(qc, key=str_rep(i))
    # trivial subcategory to show table legend before plots
    cat_xcor_plot = QCCategory(
        'xcor_plot',
        parent=cat_align_enrich,
    )
    if args.xcor_plots:
        for i, plot in enumerate(args.xcor_plots):
            if plot:
                cat_xcor_plot.add_plot(plot, key=str_rep(i), size_pct=60)

    cat_tss_enrich = QCCategory(
        'tss_enrich',
        html_head='<h2>TSS enrichment (filtered/deduped BAM)</h2>',
        html_foot="""
            <p>Open chromatin assays should show enrichment in open chromatin sites, such as
            TSS's. An average TSS enrichment in human (hg19) is above 6. A strong TSS enrichment is
            above 10. For other references please see https://www.encodeproject.org/atac-seq/</p><br>
        """,
        parser=parse_tss_enrich_qc,
        map_key_desc=MAP_KEY_DESC_TSS_ENRICH_QC,
        parent=cat_align_enrich
    )

    if args.tss_enrich_qcs:
        for i, qc in enumerate(args.tss_enrich_qcs):
            if qc:
                cat_tss_enrich.add_log(qc, key=str_rep(i))

    if args.tss_large_plots:
        for i, plot in enumerate(args.tss_large_plots):
            if plot:
                cat_tss_enrich.add_plot(plot, key=str_rep(i))

    cat_jsd = QCCategory(
        'jsd',
        html_head='<h2>Jensen-Shannon distance (filtered/deduped BAM)</h2>',
        parser=parse_jsd_qc,
        map_key_desc=MAP_KEY_DESC_JSD_QC,
        parent=cat_align_enrich
    )
    if args.jsd_plot:
        plot = args.jsd_plot[0]
        cat_jsd.add_plot(plot, size_pct=40)
    if args.jsd_qcs:
        for i, qc in enumerate(args.jsd_qcs):
            if qc:
                cat_jsd.add_log(qc, key=str_rep(i))

    return cat_align_enrich


def make_cat_peak_enrich(args, cat_root):
    cat_peak_enrich = QCCategory(
        'peak_enrich',
        html_head='<h1>Peak enrichment</h1><hr>',
        parent=cat_root
    )

    cat_frip = QCCategory(
        'frac_reads_in_peaks',
        html_head='<h2>Fraction of reads in peaks (FRiP)</h2>',
        html_foot="""
            <div id='help-FRiP'>
            For {peak_caller} raw peaks:<br>
            <p><ul>
            <li>repX: Peak from true replicate X </li>
            <li>repX-prY: Peak from Yth pseudoreplicates from replicate X </li>
            <li>pooled: Peak from pooled true replicates (pool of rep1, rep2, ...) </li>
            <li>pooled-pr1: Peak from 1st pooled pseudo replicate (pool of rep1-pr1, rep2-pr1, ...)</li>
            <li>pooled-pr2: Peak from 2nd pooled pseudo replicate (pool of rep1-pr2, rep2-pr2, ...)</li>
            </ul></p>
            <br>
            For overlap/IDR peaks:<br>
            <p><ul>
            <li>repX_vs_repY: Comparing two peaks from true replicates X and Y </li>
            <li>repX-pr1_vs_repX-pr2: Comparing two peaks from both pseudoreplicates from replicate X </li>
            <li>pooled-pr1_vs_pooled-pr2: Comparing two peaks from 1st and 2nd pooled pseudo replicates </li>
            </ul></p>
            </div>
        """.format(
            peak_caller=args.peak_caller),
        parent=cat_peak_enrich,
    )

    # raw peaks
    cat_frip_call_peak = QCCategory(
        args.peak_caller,
        html_head='<h3>FRiP for {} raw peaks</h3>'.format(args.peak_caller),
        parser=parse_frip_qc,
        map_key_desc=MAP_KEY_DESC_FRIP_QC,
        parent=cat_frip
    )
    if args.frip_qcs:
        for i, qc in enumerate(args.frip_qcs):
            if qc:
                cat_frip_call_peak.add_log(qc, key=str_rep(i))
    if args.frip_qcs_pr1:
        for i, qc in enumerate(args.frip_qcs_pr1):
            if qc:
                cat_frip_call_peak.add_log(qc, key=str_rep(i) + '-pr1')
    if args.frip_qcs_pr2:
        for i, qc in enumerate(args.frip_qcs_pr2):
            if qc:
                cat_frip_call_peak.add_log(qc, key=str_rep(i) + '-pr2')
    if args.frip_qc_pooled:
        cat_frip_call_peak.add_log(args.frip_qc_pooled[0], key='pooled')
    if args.frip_qc_ppr1:
        cat_frip_call_peak.add_log(args.frip_qc_ppr1[0], key='pooled-pr1')
    if args.frip_qc_ppr2:
        cat_frip_call_peak.add_log(args.frip_qc_ppr2[0], key='pooled-pr2')

    # overlap
    cat_frip_overlap = QCCategory(
        'overlap',
        html_head='<h3>FRiP for overlap peaks</h3>',
        parser=parse_frip_qc,
        map_key_desc=MAP_KEY_DESC_FRIP_QC,
        parent=cat_frip
    )
    if args.frip_overlap_qcs:
        num_rep = infer_n_from_nC2(len(args.frip_overlap_qcs))
        for i, qc in enumerate(args.frip_overlap_qcs):
            if qc:
                cat_frip_overlap.add_log(
                    qc, key=infer_pair_label_from_idx(num_rep, i))
    if args.frip_overlap_qcs_pr:
        for i, qc in enumerate(args.frip_overlap_qcs_pr):
            if qc:
                cat_frip_overlap.add_log(
                    qc, key='rep{X}-pr1_vs_rep{X}-pr2'.format(X=i + 1))
    if args.frip_overlap_qc_ppr:
        cat_frip_overlap.add_log(args.frip_overlap_qc_ppr[0],
                                 key='pooled-pr1_vs_pooled-pr2')

    # IDR
    cat_frip_idr = QCCategory(
        'idr',
        html_head='<h3>FRiP for IDR peaks</h3>',
        parser=parse_frip_qc,
        map_key_desc=MAP_KEY_DESC_FRIP_QC,
        parent=cat_frip
    )
    if args.frip_idr_qcs:
        num_rep = infer_n_from_nC2(len(args.frip_idr_qcs))
        for i, qc in enumerate(args.frip_idr_qcs):
            if qc:
                cat_frip_idr.add_log(
                    qc, key=infer_pair_label_from_idx(num_rep, i))
    if args.frip_idr_qcs_pr:
        for i, qc in enumerate(args.frip_idr_qcs_pr):
            if qc:
                cat_frip_idr.add_log(
                    qc, key='rep{X}-pr1_vs_rep{X}-pr2'.format(X=i+1))
    if args.frip_idr_qc_ppr:
        cat_frip_idr.add_log(args.frip_idr_qc_ppr[0],
                             key='pooled-pr1_vs_pooled-pr2')

    cat_annot_enrich = QCCategory(
        'frac_reads_in_annot',
        html_head='<h2>Annotated genomic region enrichment</h2>',
        html_foot="""
            <p>Signal to noise can be assessed by considering whether reads are falling into
            known open regions (such as DHS regions) or not. A high fraction of reads
            should fall into the universal (across cell type) DHS set. A small fraction
            should fall into the blacklist regions. A high set (though not all) should
            fall into the promoter regions. A high set (though not all) should fall into
            the enhancer regions. The promoter regions should not take up all reads, as
            it is known that there is a bias for promoters in open chromatin assays.</p><br>
        """,
        parser=parse_annot_enrich_qc,
        map_key_desc=MAP_KEY_DESC_ANNOT_ENRICH_QC,
        parent=cat_peak_enrich,
    )

    if args.annot_enrich_qcs:
        for i, qc in enumerate(args.annot_enrich_qcs):
            if qc:
                cat_annot_enrich.add_log(qc, key=str_rep(i))

    return cat_peak_enrich


def make_cat_etc(args, cat_root):
    cat_etc = QCCategory(
        'etc',
        html_head='<h1>Other quality metrics</h1><hr>',
        parent=cat_root
    )

    cat_roadmap = QCCategory(
        'roadmap',
        html_head='<h2>Comparison to Roadmap DNase</h2>',
        html_foot="""
            <p>This bar chart shows the correlation between the Roadmap DNase samples to
            your sample, when the signal in the universal DNase peak region sets are
            compared. The closer the sample is in signal distribution in the regions
            to your sample, the higher the correlation.</p><br>
        """,
        parent=cat_etc
    )
    if args.roadmap_compare_plots:
        for i, plot in enumerate(args.roadmap_compare_plots):
            if plot:
                cat_roadmap.add_plot(plot, key=str_rep(i), size_pct=50)

    return cat_etc


def main():
    # read params
    log.info('Parsing QC logs and reading QC plots...')
    args = parse_arguments()

    # make a root QCCategory
    cat_root = make_cat_root(args)

    # make QCCategory for each category
    make_cat_align(args, cat_root)
    make_cat_lib_complexity(args, cat_root)
    make_cat_replication(args, cat_root)
    make_cat_peak_stat(args, cat_root)
    make_cat_align_enrich(args, cat_root)
    make_cat_peak_enrich(args, cat_root)
    make_cat_etc(args, cat_root)

    log.info('Creating HTML report...')
    write_txt(args.out_qc_html, cat_root.to_html())

    log.info('Creating QC JSON file...')
    j = cat_root.to_dict()
    write_txt(args.out_qc_json, json.dumps(j, indent=4))

    if args.qc_json_ref:
        log.info('Comparing QC JSON file with reference...')
        # exclude general section from comparing
        # because it includes metadata like date, pipeline_ver, ...
        # we want to compare actual quality metrics only
        j.pop('general')
        # exclude JSD (last 3 columns are random)
        # JSD is tested in task level test.
        if 'align_enrich' in j and 'jsd' in j['align_enrich']:
            j['align_enrich'].pop('jsd')
        with open(args.qc_json_ref, 'r') as fp:
            j_ref = json.load(fp, object_pairs_hook=OrderedDict)
            if 'general' in j_ref:
                j_ref.pop("general")
            if 'align_enrich' in j_ref and 'jsd' in j_ref['align_enrich']:
                j_ref['align_enrich'].pop('jsd')
            match_qc_json_ref = j == j_ref
    else:
        match_qc_json_ref = False

    run_shell_cmd('echo {} > qc_json_ref_match.txt'.format(match_qc_json_ref))

    log.info('All done.')


if __name__ == '__main__':
    main()
