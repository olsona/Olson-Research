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


    # properly format input file
    f = open(inputFile, 'r')
    baseName = inputFile.rsplit(".",1)[0]
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
                    fN.write('>: '+ln.rstrip()[m:]+'\t')
                    s = 1
                else: # make new line
                    fN.write('\n>: '+ln.rstrip()[m:]+'\t')
            else: # genetic lines
                fN.write(ln.rstrip())
            ln = f.readline()

    # separate out files by size, using sepSizeListDownUp.pl
    workingFile = newName
    print workingFile
    genePath = newName.rsplit("/",1)[0]+"/"
    #for i in range(len(coolingSchedule)):
    i = 0
    thr = int(coolingSchedule[i])
    bgr = baseName + "_next"
    smlr = baseName + "_{!s}".format(thr)
    os.system("perl sepSizeListDownUp.pl {!s} {!s} {!s} {!s} {!s}".format(thr*1000, genePath, workingFile, smlr, bgr))
    workingFile = bgr
    print workingFile

    # main loop: iterate through cooling schedule, and once matches are made, concatenate each seed (pseudo)contig with matched contigs to make next round

    # process results from main loop to get initial "trees"


if __name__ == "__main__":
    main(sys.argv[1:])