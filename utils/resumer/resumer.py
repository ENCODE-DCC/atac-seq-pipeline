#!/usr/bin/env python2

# written by Jin Lee, 2016

import argparse
import json
from collections import OrderedDict

def parse_arguments():
    parser = argparse.ArgumentParser(prog='Resumer for ENCODE ATAC/Chip-Seq pipelines',
                                        description='Parse cromwell\'s metadata JSON file and generate a new input JSON file '
                                        'to resume a pipeline from where it left off.')
    parser.add_argument('metadata_json_file', type=str, help='Cromwell metadata JSON file from a previous failed run.')
    parser.add_argument('output_def_json_file', type=str, help='Output definition JSON file for your pipeline. '
                        'Use atac.json/chip.json for ATAC-Seq/ChIP-Seq pipelines. '
                        'You can also use your own JSON file for your pipeline. '
                        'Entries in "Array[Object]" is for Array[Object] in an input JSON. This is useful to take outputs from a scatter block. '
                        'For example, the 1st entry of "Array[Object]" in chip.json is "chip.bwa" : {"bam" : "chip.bams", "flagstat_qc" : "chip.flagstat_qcs"}. '
                        'chip.flagstat_qcs : [...(taken from an output of chip.bwa.flagstat_qc)...] will be added to your new input JSON. '
                        'For example, the 1st entry of "Object" in chip.json is "chip.pool_ta" : {"ta_pooled" : "chip.ta_pooled"}. '
                        'chip.ta_pooled : "(taken from an output of chip.pool_ta.ta_pooled)" will be added to your new input JSON. ')
    args = parser.parse_args()
    return args

def read_json_file(json_file):
    with open(json_file,'r') as fp:
        return json.load(fp, object_pairs_hook=OrderedDict)

def parse_cromwell_metadata_json_file(json_file):    
    metadata_json = read_json_file(json_file)

    workflow_id = metadata_json['labels']['cromwell-workflow-id'].replace('cromwell-','')
    org_input_json = json.loads(metadata_json['submittedFiles']['inputs'], object_pairs_hook=OrderedDict)
    calls = metadata_json['calls']

    return workflow_id, org_input_json, calls

def find_output_of_successful_calls(calls, output_def_json):
    result = OrderedDict()

    if 'Array[Object]' in output_def_json:
        for call_name in output_def_json['Array[Object]']:
            if call_name in calls:
                call = calls[call_name] # call is a list of the same task for multiple replicates
                failed = False
                for i, c in enumerate(call): # i = 0-based replicate id
                    if c['executionStatus']!='Done':
                        failed = True
                        break
                if not failed:
                    for key in output_def_json['Array[Object]'][call_name]:
                        wdl_var_name = output_def_json['Array[Object]'][call_name][key]
                        result[wdl_var_name] = [call[i]['outputs'][key] for i, _ in enumerate(call)]

    if 'Object' in output_def_json:
        for call_name in output_def_json['Object']:
            if call_name in calls:
                call = calls[call_name] # call is a list of the same task for multiple replicates
                failed = False
                for i, c in enumerate(call): # i = 0-based replicate id
                    if c['executionStatus']!='Done':
                        failed = True
                        break
                if not failed:
                    assert(len(call)==1)
                    for key in output_def_json['Object'][call_name]:
                        wdl_var_name = output_def_json['Object'][call_name][key]
                        result[wdl_var_name] = call[0]['outputs'][key]

    return result

def main():
    args = parse_arguments()

    workflow_id, org_input_json, calls = parse_cromwell_metadata_json_file(args.metadata_json_file)

    output_def_json = read_json_file(args.output_def_json_file)

    new_input_json = find_output_of_successful_calls(calls, output_def_json)

    # merge new input json over original input json    
    for key in new_input_json:
        org_input_json[key] = new_input_json[key]

    with open('resume.{}.json'.format(workflow_id),'w') as fp:
        fp.write(json.dumps(org_input_json, indent=4))

if __name__=='__main__':
    main()
