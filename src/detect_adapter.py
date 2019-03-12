# written by Nathan Boley, from https://github.com/nboley/GGR_code

import sys
import gzip

VERBOSE = False

adapters = {
    'Illumina': b'AGATCGGAAGAGC',
    'Nextera ': b'CTGTCTCTTATA',
    'smallRNA': b'TGGAATTCTCGG'
}

def open_gz(fname):
    return gzip.open(fname) if fname.endswith('.gz') else open(fname,'r')

def detect_adapters_and_cnts(fname, max_n_lines=1000000):
    adapter_cnts = {
        'Illumina': 0,
        'Nextera ': 0,
        'smallRNA': 0
    }

    with open_gz(fname) as fp:
        # read the first million sequences or to the end of the while -- whichever
        # comes first, and then use the adapter for trimming which was found to
        # occur most often
        for seq_index, line in enumerate(fp):
            if seq_index >= max_n_lines: break
            if seq_index%4 != 1: continue
            for key in adapters:
                if line.find(adapters[key]) > -1:
                    adapter_cnts[key] += 1

    observed_adapters = [
        adapter for adapter, cnt in sorted(
            adapter_cnts.items(), key=lambda x: -x[1])
        if cnt > 0
    ]
    return observed_adapters, adapter_cnts, seq_index//4

def detect_most_likely_adapter(fname):
    observed_adapters, adapter_cnts, n_obs_adapters = detect_adapters_and_cnts(fname)
    if observed_adapters:
        best_adapter = observed_adapters[0]
    else:
        best_adapter = ""

    if VERBOSE:
        print("\n\nAUTO-DETECTING ADAPTER TYPE\n===========================")
        print("Attempting to auto-detect adapter type from the first 1 million sequences of the first file (>> {} <<)\n".format(
            fname)
        )
        print("Found perfect matches for the following adapter sequences:")
        print("Adapter type\tCount\tSequence\tSequences analysed\tPercentage")
        for adapter in observed_adapters:
            print("{}\t{}\t{}\t{}\t\t\t{:.2%}".format(
                adapter,
                adapter_cnts[adapter],
                adapters[adapter].decode(),
                n_obs_adapters,
                adapter_cnts[adapter]/n_obs_adapters)
            )
    if best_adapter:
        return adapters[best_adapter].decode()
    else:
        return ""

def main():
    global VERBOSE
    VERBOSE = False
    best_adapter = detect_most_likely_adapter(sys.argv[1])
    print(best_adapter)

if __name__ == '__main__':
    main()
