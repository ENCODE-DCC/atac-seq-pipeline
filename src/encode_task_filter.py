#!/usr/bin/env python

# ENCODE DCC filter wrapper
# Author: Jin Lee (leepc12@gmail.com)

import sys
import os
import argparse
from encode_lib_common import (
    copy_f_to_dir, log, ls_l, mkdir_p, rm_f, run_shell_cmd, strip_ext,
    strip_ext_bam)
from encode_lib_genomic import (
    locate_picard, remove_chrs_from_bam, samstat, samtools_index,
    samtools_name_sort, bam_is_empty,
    get_samtools_res_param)


def parse_arguments():
    parser = argparse.ArgumentParser(
        prog='ENCODE DCC filter.')
    parser.add_argument('bam', type=str,
                        help='Path for raw BAM file.')
    parser.add_argument(
        '--dup-marker', type=str, choices=['picard', 'sambamba'],
        default='picard',
        help='Dupe marker for filtering mapped reads in BAM.')
    parser.add_argument('--mapq-thresh', default=30, type=int,
                        help='Threshold for low MAPQ reads removal.')
    parser.add_argument('--no-dup-removal', action="store_true",
                        help='No dupe reads removal when filtering BAM.')
    parser.add_argument('--paired-end', action="store_true",
                        help='Paired-end BAM.')
    parser.add_argument('--multimapping', default=0, type=int,
                        help='Multimapping reads.')
    parser.add_argument(
        '--filter-chrs', nargs='*',
        help='Chromosomes to be filtered for final (nodup/filt) BAM.')
    parser.add_argument('--chrsz', type=str,
                        help='2-col chromosome sizes file.')
    parser.add_argument('--mito-chr-name', default='chrM',
                        help='Mito chromosome name.')
    parser.add_argument('--nth', type=int, default=1,
                        help='Number of threads to parallelize.')
    parser.add_argument('--mem-gb', type=float,
                        help='Max. memory for samtools sort in GB. '
                        'It should be total memory for this task (not memory per thread).')
    parser.add_argument('--picard-java-heap',
                        help='Picard\'s Java max. heap: java -jar picard.jar '
                             '-Xmx[MAX_HEAP]')
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


def rm_unmapped_lowq_reads_se(bam, multimapping, mapq_thresh, nth, mem_gb, out_dir):
    """There are pipes with multiple samtools commands.
    For such pipes, use multiple threads (-@) for only one of them.
    Priority is on sort > index > fixmate > view.
    """
    prefix = os.path.join(out_dir,
                          os.path.basename(strip_ext_bam(bam)))
    filt_bam = '{}.filt.bam'.format(prefix)

    if multimapping:
        qname_sort_bam = samtools_name_sort(bam, nth, mem_gb, out_dir)

        run_shell_cmd(
            'samtools view -h {qname_sort_bam} | '
            '$(which assign_multimappers.py) -k {multimapping} | '
            'samtools view -F 1804 -Su /dev/stdin | '
            'samtools sort /dev/stdin -o {filt_bam} -T {prefix} {res_param}'.format(
                qname_sort_bam=qname_sort_bam,
                multimapping=multimapping,
                filt_bam=filt_bam,
                prefix=prefix,
                res_param=get_samtools_res_param('sort', nth=nth, mem_gb=mem_gb),
            )
        )
        rm_f(qname_sort_bam)  # remove temporary files
    else:
        run_shell_cmd(
            'samtools view -F 1804 -q {mapq_thresh} -u {bam} | '
            'samtools sort /dev/stdin -o {filt_bam} -T {prefix} {res_param}'.format(
                mapq_thresh=mapq_thresh,
                bam=bam,
                filt_bam=filt_bam,
                prefix=prefix,
                res_param=get_samtools_res_param('sort', nth=nth, mem_gb=mem_gb),
            )
        )

    return filt_bam


def rm_unmapped_lowq_reads_pe(bam, multimapping, mapq_thresh, nth, mem_gb, out_dir):
    """There are pipes with multiple samtools commands.
    For such pipes, use multiple threads (-@) for only one of them.
    Priority is on sort > index > fixmate > view.
    """
    prefix = os.path.join(out_dir,
                          os.path.basename(strip_ext_bam(bam)))
    filt_bam = '{}.filt.bam'.format(prefix)
    tmp_filt_bam = '{}.tmp_filt.bam'.format(prefix)
    fixmate_bam = '{}.fixmate.bam'.format(prefix)

    if multimapping:
        run_shell_cmd(
            'samtools view -F 524 -f 2 -u {bam} | '
            'samtools sort -n /dev/stdin -o {tmp_filt_bam} -T {prefix} {res_param} '.format(
                bam=bam,
                tmp_filt_bam=tmp_filt_bam,
                prefix=prefix,
                res_param=get_samtools_res_param('view', nth=nth),
            )
        )

        run_shell_cmd(
            'samtools view -h {tmp_filt_bam} | '
            '$(which assign_multimappers.py) -k {multimapping} --paired-end | '
            'samtools fixmate -r /dev/stdin {fixmate_bam} {res_param}'.format(
                tmp_filt_bam=tmp_filt_bam,
                multimapping=multimapping,
                fixmate_bam=fixmate_bam,
                res_param=get_samtools_res_param('fixmate', nth=nth),
            )
        )
    else:
        run_shell_cmd(
            'samtools view -F 1804 -f 2 -q {mapq_thresh} -u {bam} | '
            'samtools sort -n /dev/stdin -o {tmp_filt_bam} -T {prefix} {res_param}'.format(
                mapq_thresh=mapq_thresh,
                bam=bam,
                tmp_filt_bam=tmp_filt_bam,
                prefix=prefix,
                res_param=get_samtools_res_param('sort', nth=nth, mem_gb=mem_gb),
            )
        )
        run_shell_cmd(
            'samtools fixmate -r {tmp_filt_bam} {fixmate_bam} {res_param}'.format(
                tmp_filt_bam=tmp_filt_bam,
                fixmate_bam=fixmate_bam,
                res_param=get_samtools_res_param('fixmate', nth=nth),
            )
        )

    rm_f(tmp_filt_bam)

    run_shell_cmd(
        'samtools view -F 1804 -f 2 -u {fixmate_bam} | '
        'samtools sort /dev/stdin -o {filt_bam} -T {prefix} {res_param}'.format(
            fixmate_bam=fixmate_bam,
            filt_bam=filt_bam,
            prefix=prefix,
            res_param=get_samtools_res_param('sort', nth=nth, mem_gb=mem_gb),
        )
    )

    rm_f(fixmate_bam)
    return filt_bam


def mark_dup_picard(bam, out_dir, java_heap=None):  # shared by both se and pe
    prefix = os.path.join(out_dir,
                          os.path.basename(strip_ext_bam(bam)))
    # strip extension appended in the previous step
    prefix = strip_ext(prefix, 'filt')
    dupmark_bam = '{}.dupmark.bam'.format(prefix)
    dup_qc = '{}.dup.qc'.format(prefix)
    if java_heap is None:
        java_heap_param = '-Xmx4G'
    else:
        java_heap_param = '-Xmx{}'.format(java_heap)

    run_shell_cmd(
        'java {java_heap_param} -XX:ParallelGCThreads=1 '
        '-jar {picard} MarkDuplicates '
        'INPUT={bam} '
        'OUTPUT={dupmark_bam} '
        'METRICS_FILE={dup_qc} '
        'VALIDATION_STRINGENCY=LENIENT '
        'USE_JDK_DEFLATER=TRUE '
        'USE_JDK_INFLATER=TRUE '
        'ASSUME_SORTED=TRUE '
        'REMOVE_DUPLICATES=FALSE '.format(
            java_heap_param=java_heap_param,
            picard=locate_picard(),
            bam=bam,
            dupmark_bam=dupmark_bam,
            dup_qc=dup_qc,
        )
    )
    return dupmark_bam, dup_qc


def mark_dup_sambamba(bam, nth, out_dir):  # shared by both se and pe
    prefix = os.path.join(out_dir,
                          os.path.basename(strip_ext_bam(bam)))
    # strip extension appended in the previous step
    prefix = strip_ext(prefix, 'filt')
    dupmark_bam = '{}.dupmark.bam'.format(prefix)
    dup_qc = '{}.dup.qc'

    cmd = 'sambamba markdup -t {} --hash-table-size=17592186044416 '
    cmd += '--overflow-list-size=20000000 '
    cmd += '--io-buffer-size=256 {} {} 2> {}'
    cmd = cmd.format(
        nth,
        bam,
        dupmark_bam,
        dup_qc)
    run_shell_cmd(cmd)
    return dupmark_bam, dup_qc


def rm_dup_se(dupmark_bam, nth, out_dir):
    prefix = os.path.join(out_dir,
                          os.path.basename(strip_ext_bam(dupmark_bam)))
    # strip extension appended in the previous step
    prefix = strip_ext(prefix, 'dupmark')
    nodup_bam = '{}.nodup.bam'.format(prefix)

    run_shell_cmd(
        'samtools view -F 1804 -b {dupmark_bam} {res_param} > {nodup_bam}'.format(
            dupmark_bam=dupmark_bam,
            res_param=get_samtools_res_param('view', nth=nth),
            nodup_bam=nodup_bam,
        )
    )
    return nodup_bam


def rm_dup_pe(dupmark_bam, nth, out_dir):
    prefix = os.path.join(out_dir,
                          os.path.basename(strip_ext_bam(dupmark_bam)))
    # strip extension appended in the previous step
    prefix = strip_ext(prefix, 'dupmark')
    nodup_bam = '{}.nodup.bam'.format(prefix)

    run_shell_cmd(
        'samtools view -F 1804 -f 2 -b {dupmark_bam} {res_param} > {nodup_bam}'.format(
            dupmark_bam=dupmark_bam,
            res_param=get_samtools_res_param('view', nth=nth),
            nodup_bam=nodup_bam,
        )
    )
    return nodup_bam


def pbc_qc_se(bam, mito_chr_name, out_dir):
    prefix = os.path.join(out_dir,
                          os.path.basename(strip_ext_bam(bam)))
    # strip extension appended in the previous step
    prefix = strip_ext(prefix, 'dupmark')
    pbc_qc = '{}.lib_complexity.qc'.format(prefix)

    cmd2 = 'bedtools bamtobed -i {} | '
    cmd2 += 'awk \'BEGIN{{OFS="\\t"}}{{print $1,$2,$3,$6}}\' | '
    cmd2 += 'grep -v "^{}\\s" | sort | uniq -c | '
    cmd2 += 'awk \'BEGIN{{mt=0;m0=0;m1=0;m2=0}} ($1==1){{m1=m1+1}} '
    cmd2 += '($1==2){{m2=m2+1}} {{m0=m0+1}} '
    cmd2 += '{{mt=mt+$1}} END{{m1_m2=-1.0; '
    cmd2 += 'if(m2>0) m1_m2=m1/m2; m0_mt=0; '
    cmd2 += 'if (mt>0) m0_mt=m0/mt; m1_m0=0; if (m0>0) m1_m0=m1/m0; '
    cmd2 += 'printf "%d\\t%d\\t%d\\t%d\\t%f\\t%f\\t%f\\n",'
    cmd2 += 'mt,m0,m1,m2,m0_mt,m1_m0,m1_m2}}\' > {}'
    cmd2 = cmd2.format(
        bam,
        mito_chr_name,
        pbc_qc)
    run_shell_cmd(cmd2)
    return pbc_qc


def pbc_qc_pe(bam, mito_chr_name, nth, mem_gb, out_dir):
    prefix = os.path.join(out_dir,
                          os.path.basename(strip_ext_bam(bam)))
    pbc_qc = '{}.lib_complexity.qc'.format(prefix)

    nmsrt_bam = samtools_name_sort(bam, nth, mem_gb, out_dir)
    cmd3 = 'bedtools bamtobed -bedpe -i {} | '
    cmd3 += 'awk \'BEGIN{{OFS="\\t"}}{{print $1,$2,$4,$6,$9,$10}}\' | '
    cmd3 += 'grep -v "^{}\\s" | sort | uniq -c | '
    cmd3 += 'awk \'BEGIN{{mt=0;m0=0;m1=0;m2=0}} ($1==1){{m1=m1+1}} '
    cmd3 += '($1==2){{m2=m2+1}} {{m0=m0+1}} {{mt=mt+$1}} END{{m1_m2=-1.0; '
    cmd3 += 'if(m2>0) m1_m2=m1/m2; m0_mt=0; '
    cmd3 += 'if (mt>0) m0_mt=m0/mt; m1_m0=0; if (m0>0) m1_m0=m1/m0; '
    cmd3 += 'printf "%d\\t%d\\t%d\\t%d\\t%f\\t%f\\t%f\\n"'
    cmd3 += ',mt,m0,m1,m2,m0_mt,m1_m0,m1_m2}}\' > {}'
    cmd3 = cmd3.format(
        nmsrt_bam,
        mito_chr_name,
        pbc_qc)
    run_shell_cmd(cmd3)
    rm_f(nmsrt_bam)
    return pbc_qc


def main():
    # filt_bam - dupmark_bam - nodup_bam
    #          \ dup_qc      \ pbc_qc

    # read params
    args = parse_arguments()

    log.info('Initializing and making output directory...')
    mkdir_p(args.out_dir)

    # declare temp arrays
    temp_files = []  # files to deleted later at the end

    log.info('Removing unmapped/low-quality reads...')
    if args.paired_end:
        filt_bam = rm_unmapped_lowq_reads_pe(
            args.bam, args.multimapping, args.mapq_thresh,
            args.nth, args.mem_gb, args.out_dir)
    else:
        filt_bam = rm_unmapped_lowq_reads_se(
            args.bam, args.multimapping, args.mapq_thresh,
            args.nth, args.mem_gb, args.out_dir)

    log.info('Checking if filtered BAM file is empty...')

    if bam_is_empty(filt_bam, args.nth):
        help_msg = (
            'No reads found in filtered BAM. '
            'Low quality sample? '
            'Or no reads passing criteria "samtools view -F 1804"? '
            'Check samtools flags at '
            'https://broadinstitute.github.io/picard/explain-flags.html. '
        )
        if args.paired_end:
            help_msg += (
                'Or is this truely PE BAM? '
                'All unpaired SE reads could be removed by "samtools view -f 2". '
            )
        raise ValueError(help_msg)

    log.info('Marking dupes with {}...'.format(args.dup_marker))
    if args.dup_marker == 'picard':
        dupmark_bam, dup_qc = mark_dup_picard(
            filt_bam, args.out_dir, args.picard_java_heap)
    elif args.dup_marker == 'sambamba':
        dupmark_bam, dup_qc = mark_dup_sambamba(
            filt_bam, args.nth, args.out_dir)
    else:
        raise argparse.ArgumentTypeError(
            'Unsupported --dup-marker {}'.format(args.dup_marker))

    if args.no_dup_removal:
        nodup_bam = filt_bam
    else:
        temp_files.append(filt_bam)
        log.info('Removing dupes...')
        if args.paired_end:
            nodup_bam = rm_dup_pe(
                dupmark_bam, args.nth, args.out_dir)
        else:
            nodup_bam = rm_dup_se(
                dupmark_bam, args.nth, args.out_dir)
        samtools_index(dupmark_bam)
        temp_files.append(dupmark_bam+'.bai')
    temp_files.append(dupmark_bam)

    if len(args.filter_chrs) > 0:
        final_bam = remove_chrs_from_bam(nodup_bam, args.filter_chrs,
                                         args.chrsz, args.nth,
                                         args.out_dir)
        temp_files.append(nodup_bam)
    else:
        final_bam = nodup_bam

    log.info('Checking if final BAM file is empty...')
    if bam_is_empty(final_bam, args.nth):
        raise ValueError(
            'No reads found in final (filtered/deduped) BAM. '
            'Low quality sample? '
            'Or BAM with duplicates only? '
        )

    log.info('samtools index (final_bam)...')
    samtools_index(final_bam, args.nth, args.out_dir)

    log.info('samstat...')
    samstat(final_bam, args.nth, args.mem_gb, args.out_dir)

    log.info('Generating PBC QC log...')
    if args.paired_end:
        pbc_qc_pe(dupmark_bam, args.mito_chr_name, args.nth, args.mem_gb,
                  args.out_dir)
    else:
        pbc_qc_se(dupmark_bam, args.mito_chr_name, args.out_dir)

    log.info('samtools index (raw bam)...')
    bam = copy_f_to_dir(args.bam, args.out_dir)
    bai = samtools_index(bam, args.nth, args.out_dir)
    temp_files.extend([bam, bai])

    log.info('Removing temporary files...')
    rm_f(temp_files)

    log.info('List all files in output directory...')
    ls_l(args.out_dir)

    log.info('All done.')


if __name__ == '__main__':
    main()
