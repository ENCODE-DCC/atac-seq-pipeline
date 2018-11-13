#!/usr/bin/env python

# ENCODE DCC filter wrapper
# Author: Jin Lee (leepc12@gmail.com)

import sys
import os
import argparse
import multiprocessing
from encode_common_genomic import *

def parse_arguments():
    parser = argparse.ArgumentParser(prog='ENCODE DCC filter.',
                                        description='')
    parser.add_argument('bam', type=str,
                        help='Path for raw BAM file.')
    parser.add_argument('--dup-marker', type=str, choices=['picard','sambamba'],
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
    parser.add_argument('--nth', type=int, default=1,
                        help='Number of threads to parallelize.')
    parser.add_argument('--out-dir', default='', type=str,
                            help='Output directory.')
    parser.add_argument('--log-level', default='INFO', 
                        choices=['NOTSET','DEBUG','INFO',
                            'WARNING','CRITICAL','ERROR','CRITICAL'],
                        help='Log level')
    args = parser.parse_args()

    log.setLevel(args.log_level)
    log.info(sys.argv)
    return args

def rm_unmapped_lowq_reads_se(bam, multimapping, mapq_thresh, nth, out_dir):
    prefix = os.path.join(out_dir,
        os.path.basename(strip_ext_bam(bam)))
    filt_bam = '{}.filt.bam'.format(prefix)

    if multimapping:        
        # qname_sort_bam = samtools_name_sort(bam, nth, out_dir)
        qname_sort_bam = sambamba_name_sort(bam, nth, out_dir)

        cmd2 = 'samtools view -h {} | '
        cmd2 += '$(which assign_multimappers.py) -k {} | '
        cmd2 += 'samtools view -F 1804 -Su /dev/stdin | '

        # cmd2 += 'samtools sort /dev/stdin -o {} -T {} -@ {}'
        # cmd2 = cmd2.format(
        #     qname_sort_bam,
        #     multimapping,
        #     filt_bam,
        #     prefix,
        #     nth)
        cmd2 += 'sambamba sort /dev/stdin -o {} -t {}'
        cmd2 = cmd2.format(
            qname_sort_bam,
            multimapping,
            filt_bam,
            nth)
        run_shell_cmd(cmd2)
        rm_f(qname_sort_bam) # remove temporary files
    else:
        cmd = 'samtools view -F 1804 -q {} -u {} | '
        cmd += 'samtools sort /dev/stdin -o {} -T {} -@ {}'
        cmd = cmd.format(
            mapq_thresh,
            bam,
            filt_bam,
            prefix,
            nth)
        run_shell_cmd(cmd)

    return filt_bam

def rm_unmapped_lowq_reads_pe(bam, multimapping, mapq_thresh, nth, out_dir):
    prefix = os.path.join(out_dir,
        os.path.basename(strip_ext_bam(bam)))
    filt_bam = '{}.filt.bam'.format(prefix)
    tmp_filt_bam = '{}.tmp_filt.bam'.format(prefix)
    fixmate_bam = '{}.fixmate.bam'.format(prefix)

    if multimapping:
        # cmd1 = 'samtools view -F 524 -f 2 -u {} | '
        # cmd1 += 'samtools sort -n /dev/stdin -o {} -T {} -@ {} '
        # cmd1 = cmd1.format(
        #     bam,
        #     tmp_filt_bam,
        #     prefix,
        #     nth)
        # run_shell_cmd(cmd1)
        cmd1 = 'samtools view -F 524 -f 2 -u {} | '
        cmd1 += 'sambamba sort -n /dev/stdin -o {} -t {} '
        cmd1 = cmd1.format(
            bam,
            tmp_filt_bam,
            nth)
        run_shell_cmd(cmd1)

        cmd2 = 'samtools view -h {} -@ {} | '
        cmd2 += '$(which assign_multimappers.py) -k {} --paired-end | '
        cmd2 += 'samtools fixmate -r /dev/stdin {}'
        cmd2 = cmd2.format(
            tmp_filt_bam,
            nth,
            multimapping,
            fixmate_bam)
        run_shell_cmd(cmd2)
    else:
        # cmd1 = 'samtools view -F 1804 -f 2 -q {} -u {} | '
        # cmd1 += 'samtools sort -n /dev/stdin -o {} -T {} -@ {}'
        # cmd1 = cmd1.format(
        #     mapq_thresh,
        #     bam,
        #     tmp_filt_bam,
        #     prefix,
        #     nth)
        # run_shell_cmd(cmd1)
        cmd1 = 'samtools view -F 1804 -f 2 -q {} -u {} | '
        cmd1 += 'sambamba sort -n /dev/stdin -o {} -t {}'
        cmd1 = cmd1.format(
            mapq_thresh,
            bam,
            tmp_filt_bam,
            nth)
        run_shell_cmd(cmd1)

        cmd2 = 'samtools fixmate -r {} {}'
        cmd2 = cmd2.format(
            tmp_filt_bam,
            fixmate_bam)
        run_shell_cmd(cmd2)
    rm_f(tmp_filt_bam)

    # cmd = 'samtools view -F 1804 -f 2 -u {} | '
    # cmd += 'samtools sort /dev/stdin -o {} -T {} -@ {}'
    # cmd = cmd.format(
    #     fixmate_bam,
    #     filt_bam,
    #     prefix,
    #     nth)
    cmd = 'samtools view -F 1804 -f 2 -u {} | '
    cmd += 'sambamba sort /dev/stdin -o {} -t {}'
    cmd = cmd.format(
        fixmate_bam,
        filt_bam,
        nth)
    run_shell_cmd(cmd)
    rm_f(fixmate_bam)
    return filt_bam

def mark_dup_picard(bam, out_dir): # shared by both se and pe
    prefix = os.path.join(out_dir,
        os.path.basename(strip_ext_bam(bam)))
    # strip extension appended in the previous step
    prefix = strip_ext(prefix,'filt') 
    dupmark_bam = '{}.dupmark.bam'.format(prefix)
    dup_qc = '{}.dup.qc'.format(prefix)

    cmd = 'java -Xmx4G -jar '
    cmd += locate_picard()
    cmd += ' MarkDuplicates '
    # cmd = 'picard MarkDuplicates '
    cmd += 'INPUT={} OUTPUT={} '
    cmd += 'METRICS_FILE={} VALIDATION_STRINGENCY=LENIENT '
    cmd += 'ASSUME_SORTED=true REMOVE_DUPLICATES=false'
    cmd = cmd.format(
        bam,
        dupmark_bam,
        dup_qc)
    run_shell_cmd(cmd)
    return dupmark_bam, dup_qc

def mark_dup_sambamba(bam, nth, out_dir): # shared by both se and pe
    prefix = os.path.join(out_dir,
        os.path.basename(strip_ext_bam(bam)))
    # strip extension appended in the previous step
    prefix = strip_ext(prefix,'filt') 
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
    prefix = strip_ext(prefix,'dupmark') 
    nodup_bam = '{}.nodup.bam'.format(prefix)

    cmd1 = 'samtools view -@ {} -F 1804 -b {} > {}'
    cmd1 = cmd1.format(
        nth,
        dupmark_bam,
        nodup_bam)
    run_shell_cmd(cmd1)
    return nodup_bam

def rm_dup_pe(dupmark_bam, nth, out_dir):
    prefix = os.path.join(out_dir,
        os.path.basename(strip_ext_bam(dupmark_bam)))
    # strip extension appended in the previous step
    prefix = strip_ext(prefix,'dupmark') 
    nodup_bam = '{}.nodup.bam'.format(prefix)

    cmd1 = 'samtools view -@ {} -F 1804 -f 2 -b {} > {}'
    cmd1 = cmd1.format(
        nth,
        dupmark_bam,
        nodup_bam)
    run_shell_cmd(cmd1)
    return nodup_bam

def pbc_qc_se(bam, out_dir):
    prefix = os.path.join(out_dir,
        os.path.basename(strip_ext_bam(bam)))
    # strip extension appended in the previous step
    prefix = strip_ext(prefix,'dupmark') 
    pbc_qc = '{}.pbc.qc'.format(prefix)

    cmd2 = 'bedtools bamtobed -i {} | '
    cmd2 += 'awk \'BEGIN{{OFS="\\t"}}{{print $1,$2,$3,$6}}\' | '
    cmd2 += 'grep -v "chrM" | sort | uniq -c | '
    cmd2 += 'awk \'BEGIN{{mt=0;m0=0;m1=0;m2=0}} ($1==1){{m1=m1+1}} '
    cmd2 += '($1==2){{m2=m2+1}} {{m0=m0+1}} {{mt=mt+$1}} END{{m1_m2=-1.0; '
    cmd2 += 'if(m2>0) m1_m2=m1/m2; m0_mt=0; if (mt>0) m0_mt=m0/mt; m1_m0=0; if (m0>0) m1_m0=m1/m0; '
    cmd2 += 'printf "%d\\t%d\\t%d\\t%d\\t%f\\t%f\\t%f\\n",'
    cmd2 += 'mt,m0,m1,m2,m0_mt,m1_m0,m1_m2}}\' > {}'
    cmd2 = cmd2.format(
        bam,
        pbc_qc)
    run_shell_cmd(cmd2)
    return pbc_qc

def pbc_qc_pe(bam, nth, out_dir):
    prefix = os.path.join(out_dir,
        os.path.basename(strip_ext_bam(bam)))
    pbc_qc = '{}.pbc.qc'.format(prefix)

    # nmsrt_bam = samtools_name_sort(bam, nth, out_dir)
    nmsrt_bam = sambamba_name_sort(bam, nth, out_dir)
    cmd3 = 'bedtools bamtobed -bedpe -i {} | '
    cmd3 += 'awk \'BEGIN{{OFS="\\t"}}{{print $1,$2,$4,$6,$9,$10}}\' | '
    cmd3 += 'grep -v "chrM" | sort | uniq -c | '
    cmd3 += 'awk \'BEGIN{{mt=0;m0=0;m1=0;m2=0}} ($1==1){{m1=m1+1}} '
    cmd3 += '($1==2){{m2=m2+1}} {{m0=m0+1}} {{mt=mt+$1}} END{{m1_m2=-1.0; '
    cmd3 += 'if(m2>0) m1_m2=m1/m2; m0_mt=0; if (mt>0) m0_mt=m0/mt; m1_m0=0; if (m0>0) m1_m0=m1/m0; '
    cmd3 += 'printf "%d\\t%d\\t%d\\t%d\\t%f\\t%f\\t%f\\n"'
    cmd3 += ',mt,m0,m1,m2,m0_mt,m1_m0,m1_m2}}\' > {}'
    cmd3 = cmd3.format(
        nmsrt_bam,
        pbc_qc)
    run_shell_cmd(cmd3)
    rm_f(nmsrt_bam)
    return pbc_qc

# if --no-dup-removal is on, 
# Cromwell/WDL wants to have a empty file 
# for output { File pbc_qc, File dup_qc }

def make_mito_dup_log(dupmark_bam, out_dir):
    prefix = os.path.join(out_dir, os.path.basename(strip_ext_bam(dupmark_bam)))    
    mito_dup_log = '{}.mito_dup.txt'.format(prefix)
    
    # Get the mitochondrial reads that are marked duplicates
    # cmd1 = 'echo -e "mito_dups\\t$(samtools view -f 1024 -c {} chrM)" > {}'.format(dupmark_bam, mito_dup_log)
    cmd1 = 'printf "mito_dups\\t$(samtools view -f 1024 -c {} chrM)\\n" > {}'.format(dupmark_bam, mito_dup_log)
    run_shell_cmd(cmd1)

    # cmd2 = 'echo -e "total_dups\\t$(samtools view -f 1024 -c {})" >> {}'.format(dupmark_bam, mito_dup_log)
    cmd2 = 'printf "total_dups\\t$(samtools view -f 1024 -c {})\\n" >> {}'.format(dupmark_bam, mito_dup_log)
    run_shell_cmd(cmd2)

    return mito_dup_log

def main():
    # filt_bam - dupmark_bam - nodup_bam
    #          \ dup_qc      \ pbc_qc

    # read params
    args = parse_arguments()

    log.info('Initializing and making output directory...')
    mkdir_p(args.out_dir)

    # declare temp arrays
    temp_files = [] # files to deleted later at the end

    log.info('Removing unmapped/low-quality reads...')
    if args.paired_end:
        filt_bam = rm_unmapped_lowq_reads_pe(
                args.bam, args.multimapping, args.mapq_thresh, 
                args.nth, args.out_dir)
    else:
        filt_bam = rm_unmapped_lowq_reads_se(
                args.bam, args.multimapping, args.mapq_thresh, 
                args.nth, args.out_dir)

    if args.no_dup_removal:
        nodup_bam = filt_bam        
    else:
        log.info('Marking dupes with {}...'.format(args.dup_marker))
        if args.dup_marker=='picard':
            dupmark_bam, dup_qc = mark_dup_picard(
                                filt_bam, args.out_dir)
        elif args.dup_marker=='sambamba':
            dupmark_bam, dup_qc = mark_dup_sambamba(
                                filt_bam, args.nth, args.out_dir)
        else:
            raise argparse.ArgumentTypeError(
            'Unsupported --dup-marker {}'.format(args.dup_marker))
        temp_files.append(filt_bam)

        log.info('Removing dupes...')
        if args.paired_end:
            nodup_bam = rm_dup_pe(
                        dupmark_bam, args.nth, args.out_dir)
        else:
            nodup_bam = rm_dup_se(
                        dupmark_bam, args.nth, args.out_dir)
        samtools_index(dupmark_bam)
        temp_files.append(dupmark_bam)
        temp_files.append(dupmark_bam+'.bai')

    # initialize multithreading
    log.info('Initializing multi-threading...')
    num_process = min(3,args.nth)
    log.info('Number of threads={}.'.format(num_process))
    pool = multiprocessing.Pool(num_process)

    # log.info('samtools index...')
    # ret_val_1 = pool.apply_async(samtools_index, 
    #                             (nodup_bam, args.out_dir))
    # log.info('samtools flagstat...')
    # ret_val_2 = pool.apply_async(samtools_flagstat,
    #                             (nodup_bam, args.out_dir))
    log.info('sambamba index...')
    ret_val_1 = pool.apply_async(sambamba_index, 
                                (nodup_bam, args.nth, args.out_dir))
    log.info('sambamba flagstat...')
    ret_val_2 = pool.apply_async(sambamba_flagstat,
                                (nodup_bam, args.nth, args.out_dir))

    log.info('Generating PBC QC log...')
    if not args.no_dup_removal:
        if args.paired_end:
            ret_val_3 = pool.apply_async(pbc_qc_pe,
                            (dupmark_bam,
                                max(1,args.nth-2),
                                args.out_dir))
        else:
            ret_val_3 = pool.apply_async(pbc_qc_se,
                            (dupmark_bam, args.out_dir))
            
    # gather
    nodup_bai = ret_val_1.get(BIG_INT)
    nodup_flagstat_qc = ret_val_2.get(BIG_INT)

    if not args.no_dup_removal:
        pbc_qc = ret_val_3.get(BIG_INT)

        log.info('Making mito dup log...')
        mito_dup_log = make_mito_dup_log(dupmark_bam, args.out_dir)

    log.info('Closing multi-threading...')
    pool.close()
    pool.join()

    log.info('Removing temporary files...')
    rm_f(temp_files)

    log.info('List all files in output directory...')
    ls_l(args.out_dir)

    log.info('All done.')

if __name__=='__main__':
    main()
