#!/usr/bin/env python2

# written by Jin Lee, 2016

import os
import glob
import sys
import re
import argparse
import json
import csv
from collections import OrderedDict, defaultdict

def parse_arguments():
    parser = argparse.ArgumentParser(prog='qc.json parser for ENCODE ATAC/Chip-Seq pipelines',
                                        description='Recursively find qc.json, '
                                        'parse it and make a TSV spreadsheet of all quality metrics.')
    parser.add_argument('--search-dir', type=str, default='.',
                            help='Root directory to be recursively searched for qc.json (or --qc-json-file-basename).')
    parser.add_argument('--out-file', type=argparse.FileType('w'), default=sys.stdout, \
                            help='(Optional) Output TSV filename. Prints to STDOUT if not defined.')
    parser.add_argument('--criteria-def-json-file', type=str,
                            help='(Optional but important) Specify criteria definition JSON file. '
                            '"criteria" category will be added to the output file/STDOUT.')
    parser.add_argument('--qc-json-file-basename', type=str, default='qc.json',
                            help='(Optional) Specify QC JSON file basename to be parsed. '
                            'Files with this name will be recursively found and parsed.')
    parser.add_argument('--tsv-mapping-qc-json-path-to-title', type=str,
                            help='(Optional) Two-column TSV (ABSOLUTE path for qc.json (or --qc-json-file-basename) [TAB] title).'
                            'This is useful when you forgot to define titles for your pipeline runs.')
    parser.add_argument('--tsv-mapping-workflow-id-to-title', type=str,
                            help='(Optional) Two-column TSV (cromwell workflow ID [TAB] title).'
                            'This is useful when you forgot to define titles for your pipeline runs.')
    args = parser.parse_args()
    return args

def find_workflow_id_from_path(path):
    abspath = os.path.abspath(path)
    # *cromwell-executions/*/WORKFLOW_ID/call-*
    pattern = 'cromwell-executions\/.+\/(.+)\/call-.+'
    m = re.findall(pattern, abspath)
    return m[0] if m else None

def recursively_read_qc_jsons(search_dir, qc_json_file_basename, 
        map_qc_json_path_to_title=None, map_workflow_id_to_title=None):
    # find all qc.json recursively
    json_files = [y for x in os.walk(search_dir) \
        for y in glob.glob(os.path.join(x[0], qc_json_file_basename))]

    # read all qc.json files
    qc_jsons = []
    for json_file in json_files:        
        with open(json_file,'r') as fp:
            qc = json.load(fp, object_pairs_hook=OrderedDict)

            # backward compatibility: change qc.json format (v1.1.4->v1.1.5)
            if 'name' in qc: # < v1.1.5
                if 'ataqc' in qc:
                    paired_end = (qc['ataqc'][0]['Paired/Single-ended']=="Paired-ended")
                elif 'flagstat_qc' in qc:
                    paired_end = (qc['flagstat_qc'][0]['paired']>0)
                elif 'nodup_flagstat_qc' in qc:
                    paired_end = (qc['nodup_flagstat_qc'][0]['paired']>0)
                else:
                    paired_end = None
                qc['general'] = OrderedDict()
                qc['general']['genome'] = qc['ataqc'][0]['Genome'] if 'ataqc' in qc else None
                qc['general']['description'] = qc.pop('desc')
                qc['general']['title'] = qc.pop('name')
                qc['general']['paired_end'] = paired_end

                for category_name in qc:
                    if type(qc[category_name])==list:
                        tmp_dict = OrderedDict()
                        for i, category_item in enumerate(qc[category_name]):
                            tmp_dict['rep{}'.format(i+1)] = category_item
                        qc[category_name] = tmp_dict
                        # print(json.dumps(tmp_dict["rep1"], sort_keys=False, indent=4, separators=(',',': ')))                        

            qc['general']['rep_id'] = None

            # find workflow_id (hash string) if possible
            workflow_id = find_workflow_id_from_path(json_file)

            if workflow_id and map_workflow_id_to_title:
                qc['general']['title'] = map_workflow_id_to_title[workflow_id]
            elif map_qc_json_path_to_title:
                qc['general']['title'] = map_qc_json_path_to_title[os.path.abspath(json_file)]

            qc_jsons.append(qc)
    return qc_jsons

def read_2col_tsv(tsv): # tsv (key \t val) -> map (key:val)
    if tsv:
        with open(tsv) as fp:
            tsv = csv.reader(fp, delimiter='\t')
            return {row : tsv[row] for row in tsv}
    else:
        return None

def read_json_file(json_file):
    if json_file:        
        with open(json_file,'r') as fp:
            return json.load(fp, object_pairs_hook=OrderedDict)
    else:
        return None

def parse_rep_key_name(rep_key_name):
    # repX-pr, repX-pr1/2 goes to repX, others (ppr, rep1-rep2, ...) go to rep1
    # 'suffix' will be suffixed to the quality_metric_name
    repX = re.findall('^rep(\d+)$',rep_key_name)
    repX_repY = re.findall('^rep(\d+)-rep(\d+)$',rep_key_name)
    repX_pr = re.findall('^rep(\d+)-pr$',rep_key_name)
    repX_pr1 = re.findall('^rep(\d+)-pr1$',rep_key_name)
    repX_pr2 = re.findall('^rep(\d+)-pr2$',rep_key_name)
    pooled = re.findall('^pooled$',rep_key_name)
    ppr = re.findall('^ppr$',rep_key_name)
    ppr1 = re.findall('^ppr1$',rep_key_name)
    ppr2 = re.findall('^ppr2$',rep_key_name)
    if repX:
        rep = 'rep{}'.format(repX[0])
        suffix = ''
    elif repX_repY:
        rep = 'rep1'
        suffix = ' (rep{}-rep{})'.format(repX_repY[0][0], repX_repY[0][1])
    elif repX_pr:
        rep = 'rep{}'.format(repX_pr[0])
        suffix = ' (rep{}-pr)'.format(repX_pr[0])
    elif repX_pr1:
        rep = 'rep{}'.format(repX_pr1[0])
        suffix = ' (rep{}-pr1)'.format(repX_pr1[0])
    elif repX_pr2:
        rep = 'rep{}'.format(repX_pr2[0])
        suffix = ' (rep{}-pr2)'.format(repX_pr2[0])
    elif pooled:
        rep = 'rep1'
        suffix = ' (pooled)'
    elif ppr:
        rep = 'rep1'
        suffix = ' (ppr)'
    elif ppr1:
        rep = 'rep1'
        suffix = ' (ppr1)'
    elif ppr2:
        rep = 'rep1'
        suffix = ' (ppr2)'
    else:
        rep = 'rep1'
        suffix = ''
    return rep, suffix

def make_a_sorted_qc_json(qc):
    result = OrderedDict()
    result['general'] = qc['general']
    if 'criteria' in qc:
        result['criteria'] = qc['criteria']
    for category in qc:
        if category in ['general','criteria']:
            continue
        result[category] = qc[category]
    return result

def pretty_print_json(d):
    print(json.dumps(d, sort_keys=False, indent=4, separators=(',',': ')))

def add_criteria_category_to_qc_json(qc, criteria):
    # count number of replicates
    tmp_dict_rep = OrderedDict()
    for category in qc: # for each category
         for rep_key_name in qc[category]: # for each rep
            rep, suffix = parse_rep_key_name(rep_key_name)            
            tmp_dict_rep[rep] = None

    # add criteria if criterie definition file is given
    if criteria:
        qc['criteria'] = OrderedDict()

        for rep in tmp_dict_rep:
            qc['criteria'][rep] = OrderedDict()
            for c in criteria:
                # two quality metric items (val, pass/fail) per condition
                # read condition
                if 'condition' in criteria[c]:
                    condition = criteria[c]['condition']
                elif qc['general']['paired_end']!=None:
                    if 'condition_pe' in criteria[c] and qc['general']['paired_end']:
                        condition = criteria[c]['condition_pe']
                    elif 'condition_se' in criteria[c] and not qc['general']['paired_end']:
                        condition = criteria[c]['condition_se']
                    else:
                        condition = None
                else:
                    condition = None

                val1 = 'N/A'
                if condition:
                    val2 = 'N/A'
                try:
                    val1 = eval(criteria[c]['eval'].replace('rep?',rep))
                    if condition:
                        for key in condition:
                            cond_met = eval('{} {}'.format(val1, condition[key]))
                            if cond_met:
                                val2 = key
                                break

                except:
                    sys.stderr.write('Failed to evaluate condition ({}) of criterion ({})\n'.format(condition, c))

                qc['criteria'][rep][c] = val1
                if condition:
                    qc['criteria'][rep]['{} (QC)'.format(c)] = val2
            # pretty_print_json(qc['criteria'])
    return qc

def main():
    args = parse_arguments()

    map_qc_json_path_to_title = read_2col_tsv(args.tsv_mapping_qc_json_path_to_title)
    map_workflow_id_to_title = read_2col_tsv(args.tsv_mapping_workflow_id_to_title)
    criteria = read_json_file(args.criteria_def_json_file)

    qc_jsons = recursively_read_qc_jsons(args.search_dir, args.qc_json_file_basename, 
        map_qc_json_path_to_title, map_workflow_id_to_title)

    # parse each qc_json and add criteria category if criteria definition file is given
    for qc_json in qc_jsons:
        add_criteria_category_to_qc_json(qc_json, criteria)

    # sort qc json (general, criteria, ...)
    sorted_qc_jsons = []
    for qc_json in qc_jsons:
        sorted_qc_jsons.append(make_a_sorted_qc_json(qc_json))

    flattened_qc_jsons = []
    for qc_json in sorted_qc_jsons:
        flat_qc_json = OrderedDict()
        for cat in qc_json:
            flat_qc_json[cat] = OrderedDict()
            for rep_key_name in qc_json[cat]:
                if type(qc_json[cat][rep_key_name])==OrderedDict:                    
                    rep, suffix = parse_rep_key_name(rep_key_name)
                    for key in qc_json[cat][rep_key_name]:
                        val = qc_json[cat][rep_key_name][key]                    
                        if not rep in flat_qc_json[cat]:
                            flat_qc_json[cat][rep] = OrderedDict()
                        flat_qc_json[cat][rep][key+suffix] = val
                else:
                    flat_qc_json[cat]['rep1'] = qc_json[cat]
                    break
        # pretty_print_json(flat_qc_json)
        flattened_qc_jsons.append(flat_qc_json)

    # dict with all category and quality metric name
    all_metric_names = OrderedDict()

    for qc_json in flattened_qc_jsons:
        for category in qc_json:
            if not category in all_metric_names:
                all_metric_names[category] = OrderedDict()
            if 'rep1' in qc_json[category]:                
                for key in qc_json[category]['rep1']:
                    all_metric_names[category][key] = None
    # TSV format
    # header
    #   layer1: category
    #   layer2: quality metric name
    # cells
    #   quality metric value

    # layer1
    header_layer1 = '\t'.join([cat if i==0 else '' \
                        for cat in all_metric_names \
                            for i, metric_name in enumerate(all_metric_names[cat])])
    # layer2
    header_layer2 = '\t'.join([metric_name \
                        for cat in all_metric_names \
                            for i, metric_name in enumerate(all_metric_names[cat])])
    
    args.out_file.write(header_layer1+'\n')
    args.out_file.write(header_layer2+'\n')

    # cells = []
    for qc_json in flattened_qc_jsons:
        num_rep = max([len(qc_json[cat]) for cat in qc_json])
        for rep_id in range(1,num_rep+1):
            line = []
            for cat in all_metric_names:                
                for i, metric_name in enumerate(all_metric_names[cat]):
                    if metric_name=='rep_id':
                        line.append(str(rep_id))
                        break
                    found_elem = False
                    if cat in qc_json:                        
                        for rep in qc_json[cat]: # for each rep
                            if rep=='rep{}'.format(rep_id):
                                for key in qc_json[cat][rep]:
                                    if key==metric_name:
                                        line.append(str(qc_json[cat][rep][key]))
                                        found_elem = True
                                        break
                    if not found_elem:
                        line.append('')
            args.out_file.write('\t'.join(line)+'\n')

if __name__=='__main__':
    main()
