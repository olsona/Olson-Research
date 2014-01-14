#! /usr/bin/env python

'''bootstrap.py - wrapper class for my MS project.'''

import sys, getopt, string, os, re
from bootstrapConstants import *
from bootstrapUtils import *

def main(argv):
    # get inputs, check validity
    inputFile = ''
    outputFile = ''
    coolingSchedule = []
    raiPath = ''
    # Command line arguments
    try:
        opts, args = getopt.getopt(argv,"hi:o:r:c:p:",\
            ["ifile=","ofile=","reference=","db=","schedule=",\
            "path="])
    except getopt.GetoptError:
        print usageString
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print usageString
            print bigUsageString
            sys.exit()
        elif opt in ("-i", "--ifile"):
            inputFile = arg
        elif opt in ("-o", "--ofile"):
            outputFile = arg
        elif opt in ("-c", "--schedule"):
            coolingSchedule = [int(n) for n in arg.lstrip()[1:-1].split(',')]
        elif opt in ("-p", "--path:"):
            raiPath = arg
    if len(inputFile) == 0:
        print 'Missing argument: -i'
        print usageString
        sys.exit(2)
    elif len(outputFile) == 0:
        print 'Missing argument: -o'
        print usageString
        sys.exit(2)
    elif len(raiPath) == 0:
        print 'Missing argument: -p'
        print usageString
        sys.exit(2)
    if len(coolingSchedule) == 0:
        coolingSchedule = defaultSchedule
    # Checking validity of inputs
    if raiPath[-1] != "/":
        raiPath = raiPath + "/"
    try:
        temp = open(inputFile, 'r')
        temp.close()
    except IOError:
        print inputFile + "cannot be opened."
        sys.exit(1)
    try:
        temp = open(outputFile, 'w')
        temp.close()
    except IOError:
        print outputFile + "cannot be opened."
        sys.exit(1)
    try:
        temp = open(raiPath+"rait",'r')
        temp.close()
    except IOError:
        print raiPath+"rait cannot be opened."
        sys.exit(1)
   
        
    # Separate out all sizes
    fileSep = inputFile
    baseName = fileSep.split(".")[-2]
    myFiles = []
    rangeList = []
    for i in range(len(coolingSchedule)):
        num = coolingSchedule[i]
        f = open(fileSep, 'r')
        ln = f.readline()
        if string.find(ln,"\t") == -1:
            # convert to contiguous line AND tabbed format
            newName = baseName+"_TAB.fa"
            fN = open(newName,'w')
            s = 0
            while ln:
                if ln[0] == '>': # deal with name lines
                    m = re.search('[A-Za-z]',ln).start()
                    if s == 0: # start file
                        fN.write(ln.rstrip()[m:]+'\t')
                        s = 1
                    else: # make new line
                        fN.write('\n'+ln.rstrip()[m:]+'\t')
                else: # genetic lines
                    fN.write(ln.rstrip())
                ln = f.readline()
        if i == 0:
            bgr = "{!s}-gt{!s}k-LIST".format(baseName, num)
            rangeList.append("gt{!s}k".format(num))
        elif num == 0:
            bgr = "{!s}-lt{!s}k-LIST".format(baseName, coolingSchedule[i-1])
            rangeList.append("lt{!s}k".format(coolingSchedule[i-1]))
        else:
            bgr = "{!s}-gt{!s}k-lt{!s}k-LIST".format(baseName,\
                num, coolingSchedule[i-1])
            rangeList.append("gt{!s}k-lt{!s}k".format(num, coolingSchedule[i-1]))
        smlr = baseName + "-lt" + str(num) + "k.fa"
        pth = fileSep.rsplit("/",1)[0]+"/"
        # Produce RAI input lists for bigger contigs and fasta file of smaller contigs
        os.system("perl sepSizeListTopDown.pl {!s} {!s} {!s} {!s} {!s}".\
            format(1000*num, pth, fileSep, smlr, bgr))
        fileSep = smlr
        myFiles.append(bgr)
    myFiles.append("{!s}-lt{!s}k-LIST".format(baseName, coolingSchedule[-1]))
    rangeList.append("lt{!s}k".format(coolingSchedule[-1]))
    
    #matches = {}
    #firstSeeds = set()
    # Main grouping loop
    #for i in range(len(rangeList)-1):
        # Seed this round
    #    os.system("{!s}rait -new -o {!s}{!s}DB -i {!s}-2 >/dev/null 2>&1".format(\
    #        raiPath, pth, rangeList[i], myFiles[i]))
        # Match round of smaller contigs to database of longer contigs
    #    os.system("{!s}rai -d {!s}{!s}DB -I {!s}-1".format(raiPath, pth,\
    #        rangeList[i], myFiles[i+1]))
    #    myFileShort = format(myFiles[i+1].split("/")[-1])
    #    os.system("cp {!s}/{!s}-1.bin {!s}{!s}-1.bin".format(os.getcwd(),\
    #        myFileShort, pth, myFileShort)) # moves results to results folder
    #    os.system("rm {!s}/{!s}-1.bin".format(os.getcwd(), myFileShort))
        
        # Keep track of who was attached to what larger contig
    #    fmatch = open("{!s}{!s}-1.bin".format(pth, myFileShort),'r')
    #    for l in fmatch.readlines():
    #        [u1,u2] = l.rstrip().split(" ")
    #        if i == 0:
    #            firstSeeds.add(u2)
    #        if u2 in matches:
    #            matches[u2].append(u1)
    #        else:
    #            matches[u2] = [u1]
    
    #finalOut = open(outputFile,'w')
    #for fs in firstSeeds:
    #    finalOut.write("{!s}\n".format(fs) + "\n  ".join(str(x) for x in matches[fs])+"\n")
    #finalOut.close()
    

if __name__ == "__main__":
    main(sys.argv[1:])