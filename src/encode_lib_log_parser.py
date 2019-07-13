#!/usr/bin/env python
"""
ENCODE QC log parser wrapper which converts a log file into a dict

Author: Jin Lee (leepc12@gmail.com)
"""

from collections import OrderedDict


MAP_KEY_DESC_FRAC_MITO_QC = {
    'total_reads': 'Total reads',
    'mito_reads': 'Mito. reads',
    'frac_mito': 'Frac. of mito. reads',
    'total_nodup_reads': 'Total deduped reads',
    'mito_nodup_reads': 'Mito. deduped reads',
    'frac_nodup_mito': 'Frac. of deduped mito. reads',
    'total_dup_reads': 'Total dup reads',
    'mito_dup_reads': 'Mito. dup reads',
    'frac_dup_mito': 'Frac. of mito. dup reads',
}

def parse_frac_mito_qc(txt):
    result = OrderedDict()
    with open(txt, 'r') as fp:
        for line in fp.read().strip('\n').split('\n'):
            k, v = line.split('\t')
            result[k] = v
    return result

MAP_KEY_DESC_FLAGSTAT_QC = {
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
}

def parse_flagstat_qc(txt):
    result = OrderedDict()
    if not txt: return result
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
            if ' total ' in line:
                if ' in total ' in line:
                    tmp1 = line.split(' in total ')
                else:
                    tmp1 = line.split(' total ')
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
                mapped_pct = tmp3_2[0]  #.replace('%','')
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
                paired_properly_pct = tmp3_2[0]  #.replace('%','')
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
                singletons_pct = tmp3_2[0]  #.replace('%','')       
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
        if 'nan' not in mapped_pct and 'N/A' not in mapped_pct \
                and 'NA' not in mapped_pct:
            if '%' in mapped_pct:
                mapped_pct = mapped_pct.replace('%','')
                result['mapped_pct'] = float(mapped_pct)
            else:
                result['mapped_pct'] = 100.0 * float(mapped_pct)
        else:
            result['mapped_pct'] = 0.0
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
    if paired_properly_pct:
        if 'nan' not in paired_properly_pct and 'N/A' not in paired_properly_pct \
                and 'NA' not in paired_properly_pct:
            if '%' in paired_properly_pct:
                paired_properly_pct = paired_properly_pct.replace('%','')
                result['paired_properly_pct'] = float(paired_properly_pct)
            else:
                result['paired_properly_pct'] = 100.0 * float(paired_properly_pct)
        else:
            result['paired_properly_pct'] = 0.0
    if with_itself:
        result['with_itself'] = int(with_itself)
    if with_itself_qc_failed:
        result['with_itself_qc_failed'] = int(with_itself_qc_failed)
    if singletons:
        result['singletons'] = int(singletons)
    if singletons_qc_failed:
        result['singletons_qc_failed'] = int(singletons_qc_failed)
    if singletons_pct:
        if 'nan' not in singletons_pct and 'N/A' not in singletons_pct \
                and 'NA' not in singletons_pct:
            if '%' in singletons_pct:
                singletons_pct = singletons_pct.replace('%','')
                result['singletons_pct'] = float(singletons_pct)
            else:
                result['singletons_pct'] = 100.0 * float(singletons_pct)
        else:
            result['singletons_pct'] = 0.0
    if diff_chroms:
        result['diff_chroms'] = int(diff_chroms)
    if diff_chroms_qc_failed:
        result['diff_chroms_qc_failed'] = int(diff_chroms_qc_failed)
    return result

MAP_KEY_DESC_DUP_QC = {
    'unpaired_reads' : 'Unpaired Reads',
    'paired_reads' : 'Paired Reads',
    'unmapped_reads' : 'Unmapped Reads',
    'unpaired_dupes' : 'Unpaired Dupes',
    'paired_dupes' : 'Paired Dupes',
    'paired_opt_dupes' : 'Paired Opt. Dupes',
    'dupes_pct' : '% Dupes/100',
}

def parse_dup_qc(txt):
    result = OrderedDict()
    if not txt: return result
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
                content = line.replace(',','.')
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
        if 'PERCENT_DUPLICATION' in m:
            dupes_pct = m['PERCENT_DUPLICATION']
        else:
            dupes_pct = '0'
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
                elif paired_reads:
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

MAP_KEY_DESC_PBC_QC = {
    'total_read_pairs' : 'Total Reads (Pairs)',
    'distinct_read_pairs' : 'Distinct Reads (Pairs)',
    'one_read_pair' : 'One Read (Pair)',
    'two_read_pair' : 'Two Reads (Pairs)',
    'NRF' : 'NRF = Distinct/Total',
    'PBC1' : 'PBC1 = OnePair/Distinct',
    'PBC2' : 'PBC2 = OnePair/TwoPair'
}

def parse_pbc_qc(txt):
    result = OrderedDict()
    if not txt: return result
    with open(txt, 'r') as f:
        for line in f:
            arr = line.strip().split('\t')
            break
    result['total_read_pairs'] = to_int(arr[0])
    result['distinct_read_pairs'] = to_int(arr[1])
    result['one_read_pair'] = to_int(arr[2])
    result['two_read_pair'] = to_int(arr[3])
    result['NRF'] = to_float(arr[4])
    result['PBC1'] = to_float(arr[5])
    result['PBC2'] = to_float(arr[6])
    return result

MAP_KEY_DESC_XCOR_SCORE = {
    'num_reads' : 'Reads',
    'est_frag_len' : 'Est. Fragment Len.',
    'corr_est_frag_len' : 'Corr. Est. Fragment Len.',
    'phantom_peak' : 'Phantom Peak',
    'corr_phantom_peak' : 'Corr. Phantom Peak',
    'argmin_corr' : 'Argmin. Corr.',
    'min_corr' : 'Min. Corr.',
    'NSC' : 'NSC',
    'RSC' : 'RSC',
}

def parse_xcor_score(txt):
    result = OrderedDict()
    if not txt: return result
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

MAP_KEY_DESC_JSD_QC = {
    'pct_gen_enrich' : '% genome enriched',
    'auc' : 'AUC',
    'ch_div' : 'CHANCE divergence',
    'elbow_pt' : 'Elbow Point',
    'jsd' : 'JS Distance',
    'syn_auc' : 'Synthetic AUC',
    'syn_elbow_pt' : 'Synthetic Elbow Point',
    'syn_jsd' : 'Synthetic JS Distance',
    'syn_x_intcpt' : 'Synthetic X-intercept',
    'x_intcpt' : 'X-intercept',
    'diff_enrich' : 'diff. enrichment',
}

def parse_jsd_qc(txt):
    result = OrderedDict()
    if not txt: return result
    with open(txt, 'r') as f:
        arr = f.readlines()[0].strip().split('\t')
    result['pct_gen_enrich'] = float(arr[0])
    result['auc'] = float(arr[1])
    result['ch_div'] = float(arr[2])
    result['elbow_pt'] = float(arr[3])
    result['jsd'] = float(arr[4])
    result['syn_auc'] = float(arr[5])
    result['syn_elbow_pt'] = float(arr[6])
    result['syn_jsd'] = float(arr[7])
    # result['syn_x_intcpt'] = float(arr[8])
    # result['x_intcpt'] = float(arr[9])
    # result['diff_enrich'] = float(arr[10])
    return result

MAP_KEY_DESC_REPRODUCIBILITY_QC = {
    'Np' : 'Np',
    'Nt' : 'Nt',
    'N1' : 'N1',
    'N2' : 'N2',
    'N3' : 'N3',
    'N4' : 'N4',
    'N5' : 'N5',
    'N6' : 'N6',
    'N7' : 'N7',
    'N8' : 'N8',
    'N9' : 'N9',
    'N10' : 'N10',
    'N_opt' : 'N optimal',
    'N_consv' : 'N conservative',
    'opt_set' : 'Optimal Set',
    'consv_set' : 'Conservative Set',
    'rescue_ratio' : 'Rescue Ratio',
    'self_consistency_ratio' : 'Self Consistency Ratio',
    'reproducibility' : 'Reproducibility Test',
}

def parse_reproducibility_qc(txt):
    if not txt: return OrderedDict()
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

MAP_KEY_DESC_FRIP_QC = {
    'frip' : 'Fraction of Reads in Peaks',
}

def parse_frip_qc(txt):
    result = OrderedDict()
    if not txt: return result
    with open(txt, 'r') as f:
        frip = f.readlines()[0].strip()
    result['frip'] = float(frip)
    return result

MAP_KEY_DESC_ANNOT_ENRICH_QC = {
    'fri_dhs' : 'Fraction of Reads in universal DHS regions',
    'fri_blacklist' : 'Fraction of Reads in blacklist regions',
    'fri_prom' : 'Fraction of Reads in promoter regions',
    'fri_enh' : 'Fraction of Reads in enhancer regions',
}

def parse_annot_enrich_qc(txt):
    result = OrderedDict()
    if not txt: return result
    with open(txt, 'r') as fp:
        lines = fp.read().strip('\n').split('\n')
    for line in lines:
        key, reads, frac = line.split('\t')
        frac = to_float(frac)
        if key == 'fraction_of_reads_in_universal_DHS_regions':
            result['fri_dhs'] = frac
        elif key == 'fraction_of_reads_in_blacklist_regions':
            result['fri_blacklist'] = frac
        elif key == 'fraction_of_reads_in_promoter_regions':
            result['fri_prom'] = frac
        elif key == 'fraction_of_reads_in_enhancer_regions':
            result['fri_enh'] = frac
        else:
            raise ValueError(
                'Wrong line in annot_enrich QC file')
    return result

MAP_KEY_DESC_PICARD_EST_LIB_SIZE_QC = {
    'picard_est_lib_size' : 'Estimated library size by Picard tools',
}

def parse_picard_est_lib_size_qc(txt):
    result = OrderedDict()
    if not txt: return result
    with open(txt, 'r') as f:
        val = f.readlines()[0].strip()
    result['picard_est_lib_size'] = float(val)
    return result

MAP_KEY_DESC_TSS_ENRICH_QC = {
    'tss_enrich' : 'TSS enrichment',
}

def parse_tss_enrich_qc(txt):
    result = OrderedDict()
    if not txt: return result
    with open(txt, 'r') as f:
        val = f.readlines()[0].strip()
    result['tss_enrich'] = float(val)
    return result

MAP_KEY_DESC_NUCLEOSOMAL_QC = {
    'frac_reads_in_nfr' : 'Fraction of reads in NFR',
    'frac_reads_in_nfr_qc_pass' : 'Fraction of reads in NFR (QC pass)',
    'frac_reads_in_nfr_qc_reason' : 'Fraction of reads in NFR (QC reason)',
    'nfr_over_mono_nun_reads' : 'NFR / mono-nuc reads',
    'nfr_over_mono_nun_reads_qc_pass' : 'NFR / mono-nuc reads (QC pass)',
    'nfr_over_mono_nun_reads_qc_reason' : 'NFR / mono-nuc reads (QC reason)',
    'nfr_peak_exists' : 'Presence of NFR peak',
    'mono_nuc_peak_exists' : 'Presence of Mono-Nuc peak',
    'di_nuc_peak_exists' : 'Presence of Di-Nuc peak',
}

def parse_nucleosomal_qc(txt):
    result = OrderedDict()
    if not txt: return result
    with open(txt, 'r') as fp:
        lines = fp.read().strip('\n').split('\n')
    for line in lines:
        arr = line.split('\t')
        key = arr[0]
        if key == 'Fraction of reads in NFR':
            result['frac_reads_in_nfr'] = to_float(arr[2])
            result['frac_reads_in_nfr_qc_pass'] = to_bool(arr[1])
            result['frac_reads_in_nfr_qc_reason'] = arr[3]

        elif key == 'NFR / mono-nuc reads':
            result['nfr_over_mono_nun_reads'] = to_float(arr[2])
            result['nfr_over_mono_nun_reads_qc_pass'] = to_bool(arr[1])
            result['nfr_over_mono_nun_reads_qc_reason'] = arr[3]

        elif key == 'Presence of NFR peak':
            result['nfr_peak_exists'] = to_bool(arr[1])

        elif key == 'Presence of Mono-Nuc peak':
            result['mono_nuc_peak_exists'] = to_bool(arr[1])

        elif key == 'Presence of Di-Nuc peak':
            result['di_nuc_peak_exists'] = to_bool(arr[1])

        else:
            raise ValueError(
                'Wrong line in nucleosomal QC file')
    return result

MAP_KEY_DESC_PEAK_REGION_SIZE_QC = {
    'min_size' : 'Min size',
    '25_pct' : '25 percentile',
    '50_pct' : '50 percentile (median)',
    '75_pct' : '75 percentile',
    'max_size' : 'Max size',
    'mean' : 'Mean',
}

def parse_peak_region_size_qc(txt):
    result = OrderedDict()
    if not txt: return result
    with open(txt, 'r') as fp:
        lines = fp.read().strip('\n').split('\n')
    for line in lines:
        key, val = line.split('\t')
        if key == 'Min size':
            result['min_size'] = to_float(val)
        elif key == '25 percentile':
            result['25_pct'] = to_float(val)
        elif key == '50 percentile (median)':
            result['50_pct'] = to_float(val)
        elif key == '75 percentile':
            result['75_pct'] = to_float(val)
        elif key == 'Max size':
            result['max_size'] = to_float(val)
        elif key == 'Mean':
            result['mean'] = to_float(val)
        else:
            raise ValueError(
                'Wrong line in peak region size log file')
    return result

MAP_KEY_DESC_NUM_PEAK_QC = {
    'num_peak' : 'Number of peaks',
}

def parse_num_peak_qc(txt):
    result = OrderedDict()
    if not txt: return result
    with open(txt, 'r') as f:
        val = f.readlines()[0].strip()
    result['num_peak'] = int(val)
    return result

def to_int(var):
    try:
        return int(var)
    except ValueError:
        return None

def to_float(var):
    try:
        return float(var)
    except ValueError:
        return None

def to_number(var):
    """Convert to number or return None
    """
    try:
        if '.' in var:
            raise ValueError
        return int(var)
    except ValueError:
        try:
            return float(var)
        except ValueError:
            return None

def to_bool(var):
    return var.lower() in ('true', 't', 'ok', 'yes', '1')

