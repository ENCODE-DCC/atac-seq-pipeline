import re
import csv
import gzip

def strip_ext_fastq(fastq):
    return re.sub(r'\.(fastq|fq|Fastq|Fq)\.gz$','',fastq)

def strip_ext_bam(bam):
    return re.sub(r'\.(bam|Bam)$','',bam)

def get_read_length(fastq):
    # code extracted from Daniel Kim's ATAQC module
    # https://github.com/kundajelab/ataqc/blob/master/run_ataqc.py
    def getFileHandle(filename, mode="r"):
        if (re.search('.gz$',filename) or re.search('.gzip',filename)):
            if (mode=="r"):
                mode="rb";
            return gzip.open(filename,mode)
        else:
            return open(filename,mode)
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

def read_tsv(tsv):
    result = []
    with open(tsv,'r') as fp:        
        for row in csv.reader(fp,delimiter='\t'):
            result.append(row)
    return result

def main():
    pass

if __name__=='__main__':
    main()
    #pass