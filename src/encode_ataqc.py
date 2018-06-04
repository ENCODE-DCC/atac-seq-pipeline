#!/usr/bin/env python

# ENCODE DCC ATAQC wrapper
# Author: Daniel Kim, Jin Lee (leepc12@gmail.com)

import sys
import os
import re
import argparse
import multiprocessing
from encode_common_genomic import *
from run_ataqc import *
from encode_common_log_parser import parse_dup_qc, parse_flagstat_qc

def parse_arguments():
    parser = argparse.ArgumentParser(prog='ENCODE DCC ATAQC.', description='ATAQC')    
    parser.add_argument('--paired-end', action="store_true", help='Paired-end BAM.')
    parser.add_argument('--bowtie2-log', type=str, help='Read bowtie2 log file (from task bowtie2).')
    parser.add_argument('--read-len-log', type=str, help='Read length log file (from task bowtie2).')
    parser.add_argument('--bam', type=str, help='Raw BAM file.')    
    parser.add_argument('--flagstat-log', type=str, help='Flagstat log file for Raw BAM (from task bowtie2).')
    parser.add_argument('--nodup-bam', type=str, help='Raw BAM file (from task filter).')
    parser.add_argument('--nodup-flagstat-log', type=str, help='Flagstat log file for deduped BAM file (from task filter).')
    parser.add_argument('--pbc-log', type=str, help='PBC log file for deduped BAM file (from task filter).')
    parser.add_argument('--dup-log', type=str, help='Dup log file for deduped BAM file (from task filter).')
    parser.add_argument('--mito-dup-log', type=str, help='Mito dup log file (from task filter).')
    parser.add_argument('--ta', type=str, help='TAG-ALIGN file (from task bam2ta).')
    parser.add_argument('--bigwig', type=str, help='BIGWIG file (from task macs2).')
    parser.add_argument('--peak', type=str, help='Raw NARROWPEAK file (from task macs2).')
    parser.add_argument('--overlap-peak', type=str, help='Overlapping NARROWPEAK file (from task overlap).')
    parser.add_argument('--idr-peak', type=str, help='IDR NARROWPEAK file (from task idr).')
    parser.add_argument('--ref-fa', type=str, help='Reference fasta file.')
    parser.add_argument('--chrsz', type=str, help='2-col chromosome sizes file.')
    parser.add_argument('--tss-enrich', type=str, help='TSS enrichment definition bed file.')
    parser.add_argument('--dnase', type=str, help='DNase definition bed file.')
    parser.add_argument('--blacklist', type=str, help='Blacklist bed file.')
    parser.add_argument('--prom', type=str, help='Promoter definition bed file.')
    parser.add_argument('--enh', type=str, help='Enhancer definition bed file.')
    parser.add_argument('--reg2map', type=str, help='Reg2map file.')
    parser.add_argument('--reg2map-bed', type=str, help='Reg2map bed file.')
    parser.add_argument('--roadmap-meta', type=str, help='Roadmap metadata file.')
    parser.add_argument('--out-dir', default='', type=str, help='Output directory.')
    parser.add_argument('--log-level', default='INFO', help='Log level',
                        choices=['NOTSET','DEBUG','INFO','WARNING','CRITICAL','ERROR','CRITICAL'])
    args = parser.parse_args()
    log.setLevel(args.log_level)
    log.info(sys.argv)
    return args

def ataqc():
    # read params
    args = parse_arguments()

    # prefix
    if args.bam:
        prefix = strip_ext_bam(args.bam)
    elif args.nodup_bam:
        prefix = strip_ext_bam(args.bam)
    elif args.ta:
        prefix = strip_ext_ta(args.ta)
    elif args.peak:
        prefix = strip_ext(args.peak)
    elif args.bigwig:
        prefix = strip_ext_bigwig(args.bigwig)
    else:
        prefix = ''
    OUTPUT_PREFIX = os.path.join(args.out_dir, os.path.basename(prefix))
    # make index for bam and nodup_bam
    samtools_index(args.bam)
    samtools_index(args.nodup_bam)
    # experiment data files
    ALIGNMENT_LOG = args.bowtie2_log
    ALIGNED_BAM = args.bam
    COORDSORT_BAM = args.bam
    PBC_LOG = args.pbc_log
    DUP_LOG = args.dup_log
    FINAL_BAM = args.nodup_bam
    FINAL_BED = args.ta
    if args.idr_peak:
        PEAKS = args.peak
    else:
        PEAKS = args.overlap_peak
    NAIVE_OVERLAP_PEAKS = args.overlap_peak
    IDR_PEAKS = args.idr_peak
    BIGWIG = args.bigwig
    # genome ref data files
    GENOME = os.path.basename(args.ref_fa)
    NAME = ''
    REF = args.ref_fa
    CHROMSIZES = args.chrsz
    TSS = args.tss_enrich
    DNASE = args.dnase
    BLACKLIST = args.blacklist
    PROM = args.prom
    ENH = args.enh
    REG2MAP_BED = args.reg2map_bed if os.path.basename(args.reg2map_bed)!='null' else ''
    REG2MAP = args.reg2map
    ROADMAP_META = args.roadmap_meta

    log.info('Initializing and making output directory...')
    mkdir_p(args.out_dir)

    # Set up the log file and timing
    logging.basicConfig(filename='ataqc_test.log', level=logging.DEBUG)
    start = timeit.default_timer()

    # First check if paired/single
    # paired_status = determine_paired(FINAL_BAM)
    paired_status = "Paired-ended" if args.paired_end else "Single-ended"

    # Also get read length
    # read_len = get_read_length(FASTQ)
    if args.read_len_log:
        with open(args.read_len_log,'r') as fp:
            read_len = int(fp.read().strip())
    else:
        read_len = 0

    # Sequencing metrics: Bowtie1/2 alignment log, chrM, GC bias
    BOWTIE_STATS = get_bowtie_stats(ALIGNMENT_LOG)

    chr_m_reads, fraction_chr_m = get_chr_m(COORDSORT_BAM)

    gc_out, gc_plot, gc_summary = get_gc(FINAL_BAM,
                                         REF,
                                         OUTPUT_PREFIX)

    # Library complexity: Preseq results, NRF, PBC1, PBC2
    picard_est_library_size = get_picard_complexity_metrics(ALIGNED_BAM,
                                                            OUTPUT_PREFIX)
    preseq_data, preseq_log = run_preseq(ALIGNED_BAM, OUTPUT_PREFIX) # SORTED BAM
    #preseq_data = '/srv/scratch/dskim89/ataqc/results/2016-03-27.ENCODE_Hardison_Wold/wold/forebrain_e14-5.b1/forebrain_e14-5.b1.preseq.dat'
    #preseq_log = '/srv/scratch/dskim89/ataqc/results/2016-03-27.ENCODE_Hardison_Wold/wold/forebrain_e14-5.b1/forebrain_e14-5.b1.preseq.log'

    encode_lib_metrics = get_encode_complexity_measures(PBC_LOG)

    # Filtering metrics: duplicates, map quality
    num_mapq, fract_mapq = get_fract_mapq(ALIGNED_BAM)

    # if USE_SAMBAMBA_MARKDUP:
    #     read_dups, percent_dup = get_sambamba_dup_stats(DUP_LOG, paired_status)
    # else:    
    #     read_dups, percent_dup = get_picard_dup_stats(DUP_LOG, paired_status)
    dup_map = parse_dup_qc(DUP_LOG)
    if args.paired_end:
        read_dups = dup_map['paired_dupes']
    else:
        read_dups = dup_map['unpaired_dupes']
    percent_dup = dup_map['dupes_pct']

    # mito_dups, fract_dups_from_mito = get_mito_dups(ALIGNED_BAM,
    #                                                 OUTPUT_PREFIX,
    #                                                 paired_status,
    #                                                 use_sambamba=USE_SAMBAMBA_MARKDUP)
    with open(args.mito_dup_log,'r') as fp:
        # read mito_dup_log (TSV -> dict)
        mito_dup_map = {k: v for k,v in (map(str, line.split('\t')) for line in fp)}
        mito_dups = int(mito_dup_map['mito_dups'])
        total_dups = int(mito_dup_map['total_dups'])
        fract_dups_from_mito = mito_dups/float(total_dups)

    # [flagstat, mapped_count] = get_samtools_flagstat(ALIGNED_BAM)
    flagstat_map = parse_flagstat_qc(args.flagstat_log)
    with open(args.flagstat_log,'r') as fp:
        flagstat = fp.read()
    mapped_count = flagstat_map['mapped']

    # Final read statistics
    first_read_count, final_read_count, \
        fract_reads_left = get_final_read_count(ALIGNED_BAM,
                                                FINAL_BAM)

    # Insert size distribution - CAN'T GET THIS FOR SE FILES
    if paired_status == "Paired-ended":
        insert_data, insert_plot = get_insert_distribution(FINAL_BAM,
                                                           OUTPUT_PREFIX)
    else:
        insert_data = ''
        insert_plot = ''

    # Enrichments: V plot for enrichment
    tss_plot_file, tss_plot_large_file, tss_point_val = make_tss_plot(FINAL_BAM, # Use final to avoid duplicates
                                                                      TSS,
                                                                      OUTPUT_PREFIX,
                                                                      CHROMSIZES,
                                                                      read_len)

    # Signal to noise: reads in DHS regions vs not, reads falling
    # into blacklist regions
    reads_dnase, fract_dnase, reads_blacklist, fract_blacklist, \
        reads_prom, fract_prom, reads_enh, fract_enh, \
        reads_peaks, fract_peaks = get_signal_to_noise(FINAL_BED,
                                                       DNASE,
                                                       BLACKLIST,
                                                       PROM,
                                                       ENH,
                                                       PEAKS)

    # Also need to run n-nucleosome estimation
    if paired_status == 'Paired-ended':
        nucleosomal_qc = fragment_length_qc(read_picard_histogram(insert_data))
    else:
        nucleosomal_qc = ''

    # Peak metrics
    peak_counts = get_peak_counts(PEAKS, NAIVE_OVERLAP_PEAKS, IDR_PEAKS)
    raw_peak_summ, raw_peak_dist = get_region_size_metrics(PEAKS)
    naive_peak_summ, naive_peak_dist = get_region_size_metrics(NAIVE_OVERLAP_PEAKS)
    idr_peak_summ, idr_peak_dist = get_region_size_metrics(IDR_PEAKS)

    # Compare to roadmap
    roadmap_compare_plot = compare_to_roadmap(BIGWIG, REG2MAP_BED, REG2MAP,
                                              ROADMAP_META, OUTPUT_PREFIX)

    # Finally output the bar chart of reads
    read_count_data = [first_read_count, first_read_count*fract_mapq,
                       first_read_count*fract_mapq*(1-float(percent_dup)),
                       final_read_count]
    read_count_labels = ['Start', 'q>30', 'dups removed',
                         'chrM removed (final)']
    read_tracker_plot = track_reads(read_count_data, read_count_labels)

    # Take all this info and render the html file
    SAMPLE_INFO = OrderedDict([
        ('Sample', NAME),
        ('Genome', GENOME),
        ('Paired/Single-ended', paired_status),
        ('Read length', read_len),
    ])

    SUMMARY_STATS = OrderedDict([
        ('Read count from sequencer', first_read_count),
        ('Read count successfully aligned', mapped_count),
        ('Read count after filtering for mapping quality', num_mapq),
        ('Read count after removing duplicate reads',
            int(num_mapq - read_dups)),
        ('Read count after removing mitochondrial reads (final read count)',
            final_read_count),
    ])

    FILTERING_STATS = OrderedDict([
        ('Mapping quality > q30 (out of total)', (num_mapq, fract_mapq)),
        ('Duplicates (after filtering)', (read_dups, percent_dup)),
        ('Mitochondrial reads (out of total)', (chr_m_reads, fraction_chr_m)),
        ('Duplicates that are mitochondrial (out of all dups)',
            (mito_dups, fract_dups_from_mito)),
        ('Final reads (after all filters)', (final_read_count,
                                             fract_reads_left)),
    ])

    ENRICHMENT_PLOTS = {
        'tss': b64encode(open(tss_plot_large_file, 'rb').read())
    }

    ANNOT_ENRICHMENTS = OrderedDict([
        ('Fraction of reads in universal DHS regions', (reads_dnase,
                                                        fract_dnase)),
        ('Fraction of reads in blacklist regions', (reads_blacklist,
                                                    fract_blacklist)),
        ('Fraction of reads in promoter regions', (reads_prom, fract_prom)),
        ('Fraction of reads in enhancer regions', (reads_enh, fract_enh)),
        ('Fraction of reads in called peak regions', (reads_peaks,
                                                      fract_peaks)),
    ])

    SAMPLE = OrderedDict([
        ('Name', NAME),
        ('basic_info', SAMPLE_INFO),

        # Summary
        ('summary_stats', SUMMARY_STATS),
        ('read_tracker', read_tracker_plot),

        # Alignment statistics
        ('bowtie_stats', BOWTIE_STATS),
        ('samtools_flagstat', flagstat),

        # Filtering statistics
        ('filtering_stats', FILTERING_STATS),

        # Library complexity statistics
        ('encode_lib_complexity', encode_lib_metrics),
        ('picard_est_library_size', picard_est_library_size),
        ('yield_prediction', preseq_plot(preseq_data)),

        # Fragment length statistics
        ('fraglen_dist', fragment_length_plot(insert_data)),
        ('nucleosomal', nucleosomal_qc),

        # Peak metrics
        ('peak_counts', peak_counts),
        ('raw_peak_summ', raw_peak_summ),
        ('naive_peak_summ', naive_peak_summ),
        ('idr_peak_summ', idr_peak_summ),
        ('raw_peak_dist', raw_peak_dist),
        ('naive_peak_dist', naive_peak_dist),
        ('idr_peak_dist', idr_peak_dist),

        # GC
        ('gc_bias', plot_gc(gc_out)),

        # Annotation based statistics
        ('enrichment_plots', ENRICHMENT_PLOTS),
        ('TSS_enrichment', tss_point_val),
        ('annot_enrichments', ANNOT_ENRICHMENTS),

        # Roadmap plot
        ('roadmap_plot', roadmap_compare_plot),
    ])

    results = open('{0}_qc.html'.format(OUTPUT_PREFIX), 'w')
    results.write(html_template.render(sample=SAMPLE))
    results.close()

    # Also produce a text file of relevant stats (so that another module
    # can combine the stats) and put in using ordered dictionary
    textfile = open('{0}_qc.txt'.format(OUTPUT_PREFIX), 'w')
    for key, value in SAMPLE.iteritems():
        # Make sure to not get b64encode
        if isinstance(value, str) and (len(value) < 300) and (len(value) > 0):
            textfile.write('{0}\t{1}\n'.format(key, value))
        elif isinstance(value, int):
            textfile.write('{0}\t{1}\n'.format(key, value))
        elif isinstance(value, float):
            textfile.write('{0}\t{1}\n'.format(key, value))
        elif isinstance(value, OrderedDict):
            for dict_key, dict_value in value.iteritems():
                if isinstance(dict_value, tuple):
                    textfile.write('{0}'.format(dict_key))
                    for tuple_val in dict_value:
                        textfile.write('\t{0}'.format(tuple_val))
                    textfile.write('\n')
                else:
                    textfile.write('{0}\t{1}\n'.format(dict_key, dict_value))
        # QC tables go here
        elif isinstance(value, list):
            if 'bowtie' in value[0]: # Hack, fix this
                continue
            for result in value:
                textfile.write('{0}\t{1}\n'.format(result.metric,
                                                   result.message))
        else:
            pass
    textfile.close()

    # stop = timeit.default_timer()
    # print("Run time:", str(datetime.timedelta(seconds=int(stop - start))))

    log.info('List all files in output directory...')
    ls_l(args.out_dir)

    log.info('All done.')

if __name__=='__main__':
    ataqc()