#!/usr/bin/env python

# ENCODE DCC common functions
# Author: Jin Lee (leepc12@gmail.com)

import math
import os
import gzip
import re
import subprocess

from encode_lib_common import (
    get_gnu_sort_param,
    get_num_lines,
    get_peak_type,
    human_readable_number,
    rm_f,
    run_shell_cmd,
    strip_ext,
    strip_ext_bam,
    strip_ext_peak,
    strip_ext_ta,
    strip_ext_gz,
)


# https://github.com/samtools/samtools/blob/1.9/bam_sort.c#L70
DEFAULT_SAMTOOLS_MAX_MEM_MB_PER_THREAD = 768


def remove_chrs_from_bam(bam, chrs, chrsz, nth=1, out_dir=''):
    if len(chrs) == 0:
        raise ValueError('There must be at least one chromosome, zero found.')

    prefix = os.path.join(out_dir,
                          os.path.basename(strip_ext_bam(bam)))
    suffix = 'no_{}'.format('_'.join(chrs))
    final_bam = '{}.{}.bam'.format(prefix, suffix)
    tmp_chrsz = '{}.{}.tmp.chrsz'.format(prefix, suffix)

    # make a temp chrsz file
    cmd0 = 'zcat -f {chrsz} |'
    cmd0 += 'grep -v -P \'^({chrs})\\s\' | '
    cmd0 += 'awk \'BEGIN{{OFS="\\t"}} {{print $1,0,$2}}\' > {tmp_chrsz}'
    cmd0 = cmd0.format(
        chrsz=chrsz,
        chrs='|'.join(chrs),
        tmp_chrsz=tmp_chrsz)
    run_shell_cmd(cmd0)

    # remove chrs from BAM
    cmd1 = 'samtools view -b -L {tmp_chrsz} {bam} {res_param} > {final_bam}'
    cmd1 = cmd1.format(
        tmp_chrsz=tmp_chrsz,
        bam=bam,
        res_param=get_samtools_res_param('view', nth=nth),
        final_bam=final_bam)
    run_shell_cmd(cmd1)
    rm_f(tmp_chrsz)

    return final_bam


def samstat(bam, nth=1, mem_gb=None, out_dir=''):
    prefix = os.path.join(out_dir,
                          os.path.basename(strip_ext_bam(bam)))
    samstat_qc = '{}.samstats.qc'.format(prefix)

    run_shell_cmd(
        'samtools sort -n {bam} -T {prefix}.tmp {res_param} -O sam | '
        'SAMstats --sorted_sam_file - --outf {samstat_qc}'.format(
            bam=bam,
            prefix=prefix,
            res_param=get_samtools_res_param('sort', nth=nth, mem_gb=mem_gb),
            samstat_qc=samstat_qc,
        )
    )
    return samstat_qc


def samtools_index(bam, nth=1, out_dir=''):
    bai = '{}.bai'.format(bam)
    run_shell_cmd(
        'samtools index {bam} {res_param}'.format(
            bam=bam,
            res_param=get_samtools_res_param('index', nth=nth),
        )
    )
    if os.path.abspath(out_dir) != \
            os.path.abspath(os.path.dirname(bam)):
        return os.path.join(out_dir, os.path.basename(bai))
    else:
        return bai


def get_samtools_res_param(subcmd, nth=1, mem_gb=None):
    """Make resource parameters (-@, -m) for samtools.
        -@:
            This means number of total/additional threads.
        -m:
            This means memory per thread (for samtools sort only).
            It's clipped between 1 and DEFAULT_SAMTOOLS_MAX_MEM_MB_PER_THREAD MBs.

    Tested with samtools 1.9.
    Lower version of samtools work a bit differently.
    For such lower versions, -@ is number of threads.
    Run `samtools view --help` with your version and check if
    it is based on total or additional number of threads.

    Args:
        nth:
            Number of threads.
                - index: total threads.
                - all other sub-commands: additional threads.
        mem_gb:
            Total memory in GBs.
    """
    res_param = ''
    if subcmd == 'index':
        res_param += '-@ {num_total_threads} '.format(
            num_total_threads=nth,
        )
    else:
        res_param += '-@ {num_additional_threads} '.format(
            num_additional_threads=nth - 1,
        )

    if subcmd == 'sort':
        if nth and mem_gb:
            mem_mb_per_thread = min(
                math.floor(mem_gb * 1024.0 / nth),
                DEFAULT_SAMTOOLS_MAX_MEM_MB_PER_THREAD
            )
            res_param += '-m {mem}M '.format(mem=mem_mb_per_thread)

    return res_param


def samtools_sort(bam, nth=1, mem_gb=None, out_dir=''):
    prefix = os.path.join(out_dir,
                          os.path.basename(strip_ext_bam(bam)))
    srt_bam = '{}.srt.bam'.format(prefix)

    run_shell_cmd(
        'samtools sort {bam} -o {srt_bam} -T {prefix} {res_param}'.format(
            bam=bam,
            srt_bam=srt_bam,
            prefix=prefix,
            res_param=get_samtools_res_param('sort', nth=nth, mem_gb=mem_gb),
        )
    )
    return srt_bam


def samtools_name_sort(bam, nth=1, mem_gb=None, out_dir=''):
    prefix = os.path.join(out_dir,
                          os.path.basename(strip_ext_bam(bam)))
    nmsrt_bam = '{}.nmsrt.bam'.format(prefix)

    run_shell_cmd(
        'samtools sort -n {bam} -o {nmsrt_bam} -T {prefix} {res_param}'.format(
            bam=bam,
            nmsrt_bam=nmsrt_bam,
            prefix=prefix,
            res_param=get_samtools_res_param('sort', nth=nth, mem_gb=mem_gb),
        )
    )
    return nmsrt_bam


def bam_is_empty(bam, nth=1):
    cmd = 'samtools view -c {bam} {res_param}'.format(
        bam=bam,
        res_param=get_samtools_res_param('view', nth=nth),
    )
    return int(run_shell_cmd(cmd)) == 0


def locate_picard():
    try:
        cmd = 'which picard.jar'
        ret = run_shell_cmd(cmd)
        return ret
    except:
        try:
            # If picard.jar cannot be found, try with conda installed binary
            # This relies on that picard is correctly installed with a link
            # to the folder containing picard.jar
            cmd = 'which picard'
            picard = run_shell_cmd(cmd)
            ret = os.path.realpath(picard) + '.jar'
            if os.path.isfile(ret) and os.access(ret, os.R_OK):
                return ret
            else:
                msg = 'Potential bioconda installation of Picard tools'
                msg += ' located at:\n'
                msg += picard + '\n'
                msg += 'but the associated jar file:\n'
                msg += ret + '\n'
                msg += 'cannot be found.'
                raise Exception(msg)
        except:
            msg = 'Cannot find picard.jar or conda installation '\
                  'of Picard tools'
            raise Exception(msg)


def locate_trimmomatic():
    try:
        cmd = 'which trimmomatic.jar'
        ret = run_shell_cmd(cmd)
        return ret
    except:
        try:
            # If trimmomatic.jar cannot be found, try with conda installed binary
            # This relies on that trimmomatic is correctly installed with a link
            # to the folder containing trimmomatic.jar
            cmd = 'which trimmomatic'
            trimmomatic = run_shell_cmd(cmd)
            ret = os.path.realpath(trimmomatic) + '.jar'
            if os.path.isfile(ret) and os.access(ret, os.R_OK):
                return ret
            else:
                msg = 'Potential bioconda installation of trimmomatic'
                msg += ' located at:\n'
                msg += trimmomatic + '\n'
                msg += 'but the associated jar file:\n'
                msg += ret + '\n'
                msg += 'cannot be found.'
                raise Exception(msg)
        except:
            msg = 'Cannot find trimmomatic.jar or conda installation '\
                  'of trimmomatic'
            raise Exception(msg)


def subsample_ta_se(ta, subsample, non_mito, mito_chr_name, out_dir):
    prefix = os.path.join(out_dir,
                          os.path.basename(strip_ext_ta(ta)))
    ta_subsampled = '{}.{}{}tagAlign.gz'.format(
        prefix,
        'no_chrM.' if non_mito else '',
        '{}.'.format(human_readable_number(subsample)) if subsample > 0 else ''
    )

    # bash-only
    cmd = 'zcat -f {} | '
    if non_mito:
        # cmd += 'awk \'{{if ($1!="'+mito_chr_name+'") print $0}}\' | '
        cmd += 'grep -v \'^'+mito_chr_name+'\\b\' | '
    if subsample > 0:
        cmd += 'shuf -n {} --random-source=<(openssl enc -aes-256-ctr '
        cmd += '-pass pass:$(zcat -f {} | wc -c) -nosalt '
        cmd += '</dev/zero 2>/dev/null) | '
        cmd += 'gzip -nc > {}'
        cmd = cmd.format(
            ta,
            subsample,
            ta,
            ta_subsampled)
    else:
        cmd += 'gzip -nc > {}'
        cmd = cmd.format(
            ta,
            ta_subsampled)

    run_shell_cmd(cmd)
    return ta_subsampled


def subsample_ta_pe(ta, subsample, non_mito, mito_chr_name, r1_only, out_dir):
    prefix = os.path.join(out_dir,
                          os.path.basename(strip_ext_ta(ta)))
    if subsample % 2:
        raise ValueError(
            'Number of reads to subsample should be an even number '
            'for paired end TAG-ALIGN (BED) file. n={n}'.format(n=subsample))
    ta_subsampled = '{}.{}{}{}tagAlign.gz'.format(
        prefix,
        'no_chrM.' if non_mito else '',
        'R1.' if r1_only else '',
        '{}.'.format(human_readable_number(subsample)) if subsample > 0 else ''
    )
    ta_tmp = '{}.tagAlign.tmp'.format(prefix)

    cmd0 = 'zcat -f {} | '
    if non_mito:
        # cmd0 += 'awk \'{{if ($1!="'+mito_chr_name+'") print $0}}\' | '
        cmd0 += 'grep -v \'^'+mito_chr_name+'\\b\' | '
    cmd0 += 'sed \'N;s/\\n/\\t/\' '
    if subsample > 0:
        cmd0 += '| shuf -n {} --random-source=<(openssl enc -aes-256-ctr '
        cmd0 += '-pass pass:$(zcat -f {} | wc -c) -nosalt '
        cmd0 += '</dev/zero 2>/dev/null) > {}'
        cmd0 = cmd0.format(
            ta,
            int(subsample / 2),
            ta,
            ta_tmp)
    else:
        cmd0 += '> {}'
        cmd0 = cmd0.format(
            ta,
            ta_tmp)

    run_shell_cmd(cmd0)

    cmd = 'cat {} | '
    cmd += 'awk \'BEGIN{{OFS="\\t"}} '
    if r1_only:
        cmd += '{{printf "%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n'
        cmd += '",$1,$2,$3,$4,$5,$6}}\' | '
    else:
        cmd += '{{printf "%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n'
        cmd += '%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n",'
        cmd += '$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}}\' | '
    cmd += 'gzip -nc > {}'
    cmd = cmd.format(
        ta_tmp,
        ta_subsampled)
    run_shell_cmd(cmd)
    rm_f(ta_tmp)
    return ta_subsampled

# convert encode peak file to hammock (for Wash U browser track)


def peak_to_hammock(peak, mem_gb, out_dir):
    peak_type = get_peak_type(peak)
    prefix = os.path.join(out_dir, os.path.basename(
        strip_ext_peak(peak)))
    hammock = '{}.{}.hammock'.format(prefix, peak_type)
    hammock_tmp = '{}.tmp'.format(hammock)
    hammock_tmp2 = '{}.tmp2'.format(hammock)
    hammock_gz = '{}.gz'.format(hammock)
    hammock_gz_tbi = '{}.gz.tbi'.format(hammock)

    if get_num_lines(peak) == 0:
        run_shell_cmd(
            'zcat -f {peak} | gzip -nc > {hammock_gz}'.format(
                peak=peak,
                hammock_gz=hammock_gz,
            )
        )
        run_shell_cmd(
            'touch {hammock_gz_tbi}'.format(
                hammock_gz_tbi=hammock_gz_tbi,
            )
        )
    else:
        run_shell_cmd(
            'zcat -f {peak} | '
            'LC_COLLATE=C sort -k1,1V -k2,2n {sort_param} > {hammock_tmp}'.format(
                peak=peak,
                sort_param=get_gnu_sort_param(mem_gb * 1024 ** 3, ratio=0.5),
                hammock_tmp=hammock_tmp,
            )
        )

        with open(hammock_tmp, 'r') as fin, open(hammock_tmp2, 'w') as fout:
            id = 1
            for line in fin:
                lst = line.rstrip().split('\t')

                if peak_type == 'narrowPeak' or peak_type == 'regionPeak':
                    fout.write(
                        '{0[0]}\t{0[1]}\t{0[2]}\tscorelst:[{0[6]},{0[7]},'
                        '{0[8]}],id:{1},'.format(lst, id))
                    if len(lst[3]) > 1:
                        fout.write('name:"'+lst[3]+'",')
                    if lst[5] != '.':
                        fout.write('strand:"'+lst[5]+'",')
                    if lst[9] != '-1':
                        fout.write('sbstroke:['+lst[9]+']')
                elif peak_type == 'gappedPeak':
                    fout.write(
                        '{0[0]}\t{0[1]}\t{0[2]}\tscorelst:[{0[12]},{0[13]},'
                        '{0[14]}],id:{1},struct:{{thin:[[{0[1]},{0[2]}]],'
                        'thick:['.format(lst, id))
                    a = int(lst[1])
                    sizes = lst[10].split(',')
                    starts = lst[11].split(',')
                    for i in range(len(sizes)):
                        fout.write('[{0},{1}],'.format(
                            a+int(starts[i]), a+int(starts[i])+int(sizes[i])))
                    fout.write(']},')

                    if len(lst[3]) > 1:
                        fout.write('name:"'+lst[3]+'",')
                    if lst[5] != '.':
                        fout.write('strand:"'+lst[5]+'",')
                elif peak_type == 'broadPeak':
                    fout.write(
                        '{0[0]}\t{0[1]}\t{0[2]}\tscorelst:[{0[6]},{0[7]}],'
                        'id:{1},'.format(lst, id))
                    if len(lst[3]) > 1:
                        fout.write('name:"'+lst[3]+'",')
                    if lst[5] != '.':
                        fout.write('strand:"'+lst[5]+'",')
                else:
                    raise Exception("Unsupported peak_type {}".format(peak))
                id += 1

                fout.write('\n')

        run_shell_cmd(
            'zcat -f {hammock_tmp2} | sort -k1,1 -k2,2n {sort_param} | '
            'bgzip -cf > {hammock_gz}'.format(
                hammock_tmp2=hammock_tmp2,
                sort_param=get_gnu_sort_param(mem_gb * 1024 ** 3, ratio=0.5),
                hammock_gz=hammock_gz,
            )
        )
        run_shell_cmd(
            'tabix -f -p bed {hammock_gz}'.format(
                hammock_gz=hammock_gz,
            )
        )

        rm_f([hammock, hammock_tmp, hammock_tmp2])

    return (hammock_gz, hammock_gz_tbi)


def peak_to_starch(peak, out_dir):
    """Convert peak (BED) into starch.
    Required softwares:
        BEDOPS (tested with v2.4.39): sort-bed, starch
    """
    prefix = os.path.join(
        out_dir, os.path.basename(strip_ext_gz(peak))
    )
    starch = '{}.starch'.format(prefix)
    run_shell_cmd(
        'zcat -f {peak} | sort-bed - | starch - > {starch}'.format(
            peak=peak,
            starch=starch,
        )
    )
    return starch

def starch_to_bed_gz(starch, out_dir):
    """Convert starch into gzipped BED.
    Required softwares:
        BEDOPS (tested with v2.4.39): unstarch
    """
    prefix = os.path.join(
        out_dir, os.path.basename(strip_ext(starch))
    )
    bed_gz = '{}.bed.gz'.format(prefix)
    run_shell_cmd(
        'unstarch {starch} | gzip -nc > {bed_gz}'.format(
            starch=starch,
            bed_gz=bed_gz,
        )
    )
    return bed_gz


def peak_to_bigbed(peak, peak_type, chrsz, mem_gb, out_dir):
    prefix = os.path.join(out_dir,
                          os.path.basename(strip_ext(peak)))
    bigbed = '{}.{}.bb'.format(prefix, peak_type)
    as_file = '{}.as'.format(prefix)
    chrsz_tmp = '{}.chrsz.tmp'.format(prefix)
    bigbed_tmp = '{}.bb.tmp'.format(prefix)
    bigbed_tmp2 = '{}.bb.tmp2'.format(prefix)

    if peak_type.lower() == 'narrowpeak' or peak_type.lower() == 'regionpeak':
        as_file_contents = '''table narrowPeak
"BED6+4 Peaks of signal enrichment based on pooled, normalized (interpreted) data."
(
    string chrom;        "Reference sequence chromosome or scaffold"
    uint   chromStart;   "Start position in chromosome"
    uint   chromEnd;     "End position in chromosome"
    string name;     "Name given to a region (preferably unique). Use . if no name is assigned"
    uint   score;        "Indicates how dark the peak will be displayed in the browser (0-1000) "
    char[1]  strand;     "+ or - or . for unknown"
    float  signalValue;  "Measurement of average enrichment for the region"
    float  pValue;       "Statistical significance of signal value (-log10). Set to -1 if not used."
    float  qValue;       "Statistical significance with multiple-test correction applied (FDR -log10). Set to -1 if not used."
    int   peak;         "Point-source called for this peak; 0-based offset from chromStart. Set to -1 if no point-source called."
)
'''
        bed_param = '-type=bed6+4 -as={as_file}'.format(as_file=as_file)
    elif peak_type.lower() == 'broadpeak':
        as_file_contents = '''table broadPeak
"BED6+3 Peaks of signal enrichment based on pooled, normalized (interpreted) data."
(
    string chrom;        "Reference sequence chromosome or scaffold"
    uint   chromStart;   "Start position in chromosome"
    uint   chromEnd;     "End position in chromosome"
    string name;     "Name given to a region (preferably unique). Use . if no name is assigned."
    uint   score;        "Indicates how dark the peak will be displayed in the browser (0-1000)"
    char[1]   strand;     "+ or - or . for unknown"
    float  signalValue;  "Measurement of average enrichment for the region"
    float  pValue;       "Statistical significance of signal value (-log10). Set to -1 if not used."
    float  qValue;       "Statistical significance with multiple-test correction applied (FDR -log10). Set to -1 if not used."
)
'''
        bed_param = '-type=bed6+3 -as={as_file}'.format(as_file=as_file)
    elif peak_type.lower() == 'gappedpeak':
        as_file_contents = '''table gappedPeak
"This format is used to provide called regions of signal enrichment based on pooled, normalized (interpreted) data where the regions may be spliced or incorporate gaps in the genomic sequence. It is a BED12+3 format."
    (
    string chrom;   "Reference sequence chromosome or scaffold"
    uint chromStart;    "Pseudogene alignment start position"
    uint chromEnd;      "Pseudogene alignment end position"
    string name;        "Name of pseudogene"
    uint score;          "Score of pseudogene with gene (0-1000)"
    char[1] strand;     "+ or - or . for unknown"
    uint thickStart;    "Start of where display should be thick (start codon)"
    uint thickEnd;      "End of where display should be thick (stop codon)"
    uint reserved;      "Always zero for now"
    int blockCount;     "Number of blocks"
    int[blockCount] blockSizes; "Comma separated list of block sizes"
    int[blockCount] chromStarts; "Start positions relative to chromStart"
    float  signalValue;  "Measurement of average enrichment for the region"
    float  pValue;       "Statistical significance of signal value (-log10). Set to -1 if not used."
    float  qValue;       "Statistical significance with multiple-test correction applied (FDR). Set to -1 if not used."
)
'''
        bed_param = '-type=bed12+3 -as={as_file}'.format(as_file=as_file)
    else:
        raise Exception('Unsupported peak file type {peak_type}!'.format(peak_type=peak_type))

    # create temporary .as file
    with open(as_file, 'w') as fp:
        fp.write(as_file_contents)

    run_shell_cmd(
        'cat {chrsz} > {chrsz_tmp}'.format(
            chrsz=chrsz,
            chrsz_tmp=chrsz_tmp,
        )
    )

    run_shell_cmd(
        'zcat -f {peak} | LC_COLLATE=C sort -k1,1 -k2,2n {sort_param} | '
        'awk \'BEGIN{{OFS="\\t"}} {{if ($5>1000) $5=1000; '
        'if ($5<0) $5=0; print $0}}\' > {bigbed_tmp}'.format(
            peak=peak,
            sort_param=get_gnu_sort_param(mem_gb * 1024 ** 3, ratio=0.5),
            bigbed_tmp=bigbed_tmp,
        )
    )

    run_shell_cmd(
        'bedClip {bigbed_tmp} {chrsz_tmp} {bigbed_tmp2}'.format(
            bigbed_tmp=bigbed_tmp,
            chrsz_tmp=chrsz_tmp,
            bigbed_tmp2=bigbed_tmp2,
        )
    )

    run_shell_cmd(
        'bedToBigBed {bed_param} {bigbed_tmp2} {chrsz_tmp} {bigbed}'.format(
            bed_param=bed_param,
            bigbed_tmp2=bigbed_tmp2,
            chrsz_tmp=chrsz_tmp,
            bigbed=bigbed,
        )
    )

    # remove temporary files
    rm_f([as_file, chrsz_tmp, bigbed_tmp, bigbed_tmp2])

    return bigbed


def get_read_length(fastq):
    # code extracted from Daniel Kim's ATAQC module
    # https://github.com/kundajelab/ataqc/blob/master/run_ataqc.py
    def getFileHandle(filename, mode="r"):
        if (re.search('.gz$', filename) or re.search('.gzip', filename)):
            if (mode == "r"):
                mode = "rb"
            return gzip.open(filename, mode)
        else:
            return open(filename, mode)
    total_reads_to_consider = 1000000
    line_num = 0
    total_reads_considered = 0
    max_length = 0
    with getFileHandle(fastq, 'rb') as fp:
        for line in fp:
            if line_num % 4 == 1:
                if len(line.strip()) > max_length:
                    max_length = len(line.strip())
                total_reads_considered += 1
            if total_reads_considered >= total_reads_to_consider:
                break
            line_num += 1
    return int(max_length)


def remove_read_group(bam, out_dir='.'):
    basename = os.path.basename(strip_ext_bam(bam))
    prefix = os.path.join(out_dir, basename)
    new_bam = '{}.no_rg.bam'.format(prefix)

    cmd = 'samtools view -h {} | '
    cmd += 'grep -v "^@RG" | sed "s/\\tRG:Z:[^\\t]*//" | '
    cmd += 'samtools view -bo {} -'
    cmd = cmd.format(bam, new_bam)
    run_shell_cmd(cmd)

    return new_bam


def get_region_size_metrics(peak_file, out_dir='.'):
    '''
    From the peak file, return a plot of the region size distribution and
    the quartile metrics (summary from R)
    '''
    import pandas as pd
    import numpy as np
    import matplotlib as mpl
    mpl.use('Agg')
    from matplotlib import pyplot as plt
    from collections import OrderedDict

    basename = os.path.basename(strip_ext_peak(peak_file))
    prefix = os.path.join(out_dir, basename)
    log = '{}.peak_region_size.qc'.format(prefix)
    plot = '{}.peak_region_size.png'.format(prefix)

    # Load peak file. If it fails, return nothing as above
    peak_df = pd.read_table(peak_file, compression='gzip', header=None)

    # Subtract third column from second to get summary
    region_sizes = peak_df.iloc[:, 2] - peak_df.iloc[:, 1]

    # Summarize and store in ordered dict
    peak_summary_stats = region_sizes.describe()

    peak_size_summ = OrderedDict([
        ('Min size', peak_summary_stats['min']),
        ('25 percentile', peak_summary_stats['25%']),
        ('50 percentile (median)', peak_summary_stats['50%']),
        ('75 percentile', peak_summary_stats['75%']),
        ('Max size', peak_summary_stats['max']),
        ('Mean', peak_summary_stats['mean']),
    ])

    # Plot density diagram using matplotlib
    fig = plt.figure()
    ax = fig.add_subplot(111)

    y, binEdges = np.histogram(region_sizes, bins=100)
    bincenters = 0.5 * (binEdges[1:] + binEdges[:-1])

    # write to log file
    with open(log, 'w') as fp:
        for key, val in peak_size_summ.items():
            fp.write(key + '\t' + str(val) + '\n')

    plt.plot(bincenters, y, '-')
    filename = os.path.basename(peak_file)
    ax.set_title('Peak width distribution for {0}'.format(filename))

    # write to plot file
    fig.savefig(plot, format='png')

    return log, plot


def get_num_peaks(peak_file, out_dir='.'):
    '''
    From the peak file, return number of lines in it
    '''
    basename = os.path.basename(strip_ext_peak(peak_file))
    prefix = os.path.join(out_dir, basename)
    log = '{}.num_peak.qc'.format(prefix)

    with open(log, 'w') as fp:
        fp.write(str(get_num_lines(peak_file))+'\n')
    return log


def determine_paired(bam):
    '''
    Quick function to determine if the BAM file is paired end or single end
    '''
    num_paired_reads = int(subprocess.check_output(['samtools',
                                                    'view', '-f', '0x1',
                                                    '-c', bam]).strip())
    if num_paired_reads > 1:
        return True
    else:
        return False


def bed_clip(bed, chrsz, out_clipped_bed, no_gz=False):
    '''
    Make sure that bedClip (in USCS tools) is installed.
    Clip a BED file between 0 and chromSize (taken from 2-col chrsz file).
    bedClip exits with 255 if both start/end coordinates are
    out of valid range (0-chromSize). Otherwise, reads/peaks will be truncated.

    Args:
        no_gz:
            Do not gzip output.
    '''
    tmp_out = out_clipped_bed +  '.clip_tmp'
    cmd = 'bedClip {bed} {chrsz} {tmp_out} -truncate -verbose=2'.format(
        bed=bed,
        chrsz=chrsz,
        tmp_out=out_clipped_bed if no_gz else tmp_out)
    run_shell_cmd(cmd)

    if not no_gz:
        cmd2 = 'cat {tmp_out} | gzip -nc > {out_clipped_bed}'.format(
            tmp_out=tmp_out,
            out_clipped_bed=out_clipped_bed)
        run_shell_cmd(cmd2)

        rm_f(tmp_out)


def bam_to_pbam(bam, ref_fa, out_dir='.'):
    '''Convert BAM into pBAM.

    Requirements:
        - Python package `ptools_bin`
        - `samtools`
    '''
    prefix = os.path.join(
        out_dir,
        os.path.basename(strip_ext_bam(bam))
    )

    pbam_tmp = '{}.sorted.p.bam'.format(prefix)
    pbam = '{}.p.bam'.format(prefix)

    temp_files = []

    if ref_fa.endswith('.gz'):
        gunzipped_ref_fa = '{}.fasta'.format(
            os.path.join(out_dir, os.path.basename(strip_ext_gz(ref_fa)))
        )
        run_shell_cmd(
            'zcat -f {ref_fa} > {gunzipped_ref_fa}'.format(
                ref_fa=ref_fa,
                gunzipped_ref_fa=gunzipped_ref_fa,
            )
        )
        temp_files.append(gunzipped_ref_fa)
    else:
        gunzipped_ref_fa = ref_fa

    run_shell_cmd(
        'makepBAM_genome.sh {bam} {gunzipped_ref_fa}'.format(
            bam=bam,
            gunzipped_ref_fa=gunzipped_ref_fa,
        )
    )

    run_shell_cmd(
        'mv {pbam_tmp} {pbam}'.format(
            pbam_tmp=pbam_tmp,
            pbam=pbam,
        )
    )
    rm_f(temp_files)

    return pbam
