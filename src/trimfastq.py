##################################
#                                #
# Last modified 2017/11/08       # 
#                                #
# Georgi Marinov                 #
#                                # 
##################################

import sys
import os

READS_LOG_INTERVAL=5000000
MILLION=1000000

def run():

    if len(sys.argv) < 2:
        print('usage: python %s <inputfilename> <bpToKeep | max> [-trim5 bp] [-flowcellID flowcell] [-addEnd 1 | 2] [-replace string newstring | blank] [-renameIDs prefix] [-stdout]' % sys.argv[0])
        print('\tthe -trim5 option will trim additional bp from the 5 end, i.e. if you want the middle 36bp of 38bp reads, use 36 as bp to keep and 1 as the trim5 argument')
        print('\tUse - to specify standard input, the script will print(to standard output by default')
        print('\tThe script can read compressed files as long as they have the correct suffix - .bz2 or .gz')
        sys.exit(1)

    inputfilename = sys.argv[1]
    doMax=False
    if sys.argv[2] == 'max':
        doMax=True
        trim='max'
    else: 
        trim = int(sys.argv[2])
    outputfilename = inputfilename.split('/')[-1].split('.fastq')[0] + '.' +str(trim)+'mers.fastq'
    doFlowcellID=False

    doStdOut=True

    if '-flowcellID' in sys.argv:
        doFlowcellID=True
        flowcellID=sys.argv[sys.argv.index('-flowcellID')+1]
        if doStdOut:
            pass
        else:
            print('will include flowcell ID', flowcellID, 'in reads headers')

    doRenameIDs = False
    if '-renameIDs' in sys.argv:
        doRenameIDs = True
        RID = '@' + sys.argv[sys.argv.index('-renameIDs') + 1]

    dotrim5=False
    if '-trim5' in sys.argv:
        dotrim5=True
        trim5=int(sys.argv[sys.argv.index('-trim5')+1])
        if doStdOut:
            pass
        else:
            print('will trim ', trim5, 'bp from the 5-end')
        outputfilename = inputfilename.split('.fastq')[0] + '.' +str(trim)+'bp-5prim-trim.fastq'

    doAddEnd=False
    if '-addEnd' in sys.argv:
        doAddEnd=True
        END=sys.argv[sys.argv.index('-addEnd')+1]
        if doStdOut:
            pass
        else:
            print('will add',  '/'+END, 'to read IDs')

    doReplace=False
    if '-replace' in sys.argv:
        doReplace=True
        oldstring=sys.argv[sys.argv.index('-replace')+1]
        newstring=sys.argv[sys.argv.index('-replace')+2]
        if newstring == 'blank':
            newstring=''
        if doStdOut:
            pass
        else:
            print('will replace',  oldstring, 'with', newstring, 'in read IDs')

    if doStdOut:
        pass
    else:
        outfile = open(outputfilename, 'w')

    doStdIn = False
    if inputfilename != '-':
        if inputfilename.endswith('.bz2'):
            cmd = 'bzip2 -cd ' + inputfilename
        elif inputfilename.endswith('.gz'):
            cmd = 'gunzip -c ' + inputfilename
        else:
            cmd = 'cat ' + inputfilename
        p = os.popen(cmd, "r")
    else:
        doStdIn = True

    line = 'line'

    i=0
    shorter=0

    if dotrim5:
        i=1
        j=0
        while line != '':
            if doStdIn:
                line = sys.stdin.readline()
            else:
                line = p.readline()
            if line == '':
                break
            if i==1 and line[0]=='@':
                if doFlowcellID and flowcellID not in line:
                    ID='@'+flowcellID+'_'+line.replace(' ','_')[1:-1]+'\n'
                else:
                    ID=line.replace(' ','_')
                if doReplace:
                    ID=ID.replace(oldstring,newstring)
                if doRenameIDs:
                    ID = RID + str(j)
                if doAddEnd:
                    ID=ID.strip()+'/'+END+'\n'
                i=2
                continue
            if i==2:
                i=3
                sequence=line[trim5:len(line)].strip()
                continue
            if i==3 and line[0]=='+':
                plus='+\n'
                i=4
                continue
            if i==4:
                scores=line
                i=1
                scores=line[trim5:len(line)].strip()
                scores=scores[0:trim]
                j=j+1
                if j % READS_LOG_INTERVAL == 0:
                    if doStdOut:
                        pass
                    else:
                        print(str(j/MILLION) + 'M reads processed')
                if doMax: 
                    sequence=sequence.replace('.','N')
                else:
                    sequence=sequence[0:trim].replace('.','N')+'\n'
                if doStdOut:
                    print(ID.strip())
                    print(sequence.strip())
                    print(plus.strip())
                    print(scores)
                else:
                    outfile.write(ID.strip()+'\n')
                    outfile.write(sequence.strip()+'\n')
                    outfile.write(plus.strip()+'\n')
                    outfile.write(scores + '\n')
                continue
    else:
        i=1
        j=0
        while line != '':
            if doStdIn:
                line = sys.stdin.readline()
            else:
                line = p.readline()
            if line == '':
                break
            if i==1 and line[0]=='@':
                if doFlowcellID and flowcellID not in line:
                    ID='@'+flowcellID+'_'+line.replace(' ','_')[1:-1]+'\n'
                else:
                    ID=line.replace(' ','_')
                if doReplace:
                    ID=ID.replace(oldstring,newstring)
                if doRenameIDs:
                    ID = RID + str(j)
                if doAddEnd:
                    ID=ID.strip()+'/'+END+'\n'
                i=2
                continue
            if i==2:
                i=3
                j=j+1
                if j % READS_LOG_INTERVAL == 0:
                    if doStdOut:
                        pass
                    else:
                        print(str(j/MILLION) + 'M reads processed')
                if doMax: 
                    sequence=line
                else:
                    if len(line.strip())<trim:
                        shorter+=1
                        sequence=line.strip().replace('.','N')+'\n'
                    else:
                        sequence=line[0:trim].replace('.','N')+'\n'
                continue
            if i==3 and line[0]=='+':
                plus='+\n'
                i=4
                continue
            if i==4:
                i=1
                if doMax: 
                    scores=line
                    if doStdOut:
                        print(ID.strip())
                        print(sequence.strip())
                        print(plus.strip())
                        print(line.strip())
                    else:
                        outfile.write(ID)
                        outfile.write(sequence)
                        outfile.write(plus)
                        outfile.write(line)
                else:
                    if len(line.strip())<trim:
                        continue
                    scores=line[0:trim]+'\n'
                    if doStdOut:
                        print(ID.strip())
                        print(sequence.strip())
                        print(plus.strip())
                        print(scores.strip())
                    else:
                        outfile.write(ID)
                        outfile.write(sequence)
                        outfile.write(plus)
                        outfile.write(scores)
                continue

    if doStdOut:
        pass
    else:
        outfile.close()

    if shorter>0:
        print(shorter, 'sequences shorter than desired length')
run()

