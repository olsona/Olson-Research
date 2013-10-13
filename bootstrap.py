#! /usr/bin/env python

'''bootstrap.py - wrapper class for my MS project.'''

import sys, getopt, string, os
from bootstrapConstants import *

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
    baseName = fileSep.split(".")[0]
    myFiles = []
    for i in range(len(coolingSchedule)):
        num = coolingSchedule[i]
        f = open(fileSep, 'r')
        ln = f.readline()
        if string.find(ln,"\t") == -1:
            # convert to tabbed format using tab2fasta.pl
            newName = fileSep.split(".")[0]+"_TAB.fa"
            os.system("cat {!s} | perl -pe's/[\r\n]+$/\t/ if $i = !$i' > {!s}"\
                .format(fileSep, newName))
            fileSep = newName
        if i == 0:
            smlr = "{!s}_lt{!s}k.LIST".format(baseName, num)
        else:
            smlr = "{!s}_gt{!s}k_lt{!s}k.LIST".format(baseName,\
                coolingSchedule[i-1], num)
        bgr = baseName + "_gt" + str(num) + "k.fa"
        pth = fileSep.rsplit("/",1)[0]+"/"
        os.system("perl separateBySizeListFormat.pl {!s} {!s} {!s} {!s} {!s}".\
            format(1000*num, pth, fileSep, smlr, bgr))
        fileSep = bgr
        myFiles.append(smlr)
    myFiles.append(fileSep)
    
    # De novo clustering of first files
    print "{!s}rait -o {!s}{!s}kDB -i {!s}.2".format(raiPath, pth,\
        coolingSchedule[0], myFiles[0])
    os.system("{!s}rait -o {!s}lt{!s}kDB -i {!s}.2".format(raiPath, pth,\
        coolingSchedule[0], myFiles[0]))
    

if __name__ == "__main__":
    main(sys.argv[1:])