#! /usr/bin/env python

'''bootstrap.py - wrapper class for my MS project.'''

import sys, getopt, string, os, re, pprint
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
                    myStr = ln.rstrip()[m:]
                    fN.write('>'+string.replace(myStr,' ','_')+'\t')
                    s = 1
                else: # make new line
                    myStr = ln.rstrip()[m:]
                    fN.write('\n>'+string.replace(myStr,' ','_')+'\t')
            else: # genetic lines
                fN.write(ln.rstrip())
            ln = f.readline()
    f.close()
    fN.close()

    # separate out files by size, using sepSizeListDownUp.pl
    fNext = newName
    genePath = newName.rsplit("/",1)[0]+"/contigs/"
    ensureDir(genePath)
    leng = len(coolingSchedule)
    for i in range(leng):
        #for i in [0]:
        workingFile = fNext
        thr = int(coolingSchedule[i])
        bgr = "{!s}_{!s}_next".format(baseName,i)
        smlr = "{!s}_{!s}".format(baseName, i)
        os.system("perl sepSizeListDownUp.pl {!s} {!s} {!s} {!s} {!s}".format(thr*1000, genePath, workingFile, smlr, bgr))
        fNext = bgr

    # Make initial seed file
    fSeed = baseName+"_seed"
    os.system("perl processSeedFile.pl {!s} {!s} {!s}".format(genePath, fNext, fSeed))
    
    masterDict = {}
    roots = set()
    ct = 0
    fOut = open(outputFile,'a')

    # Main loop: iterate through cooling schedule, creating databases, making matches, and once matches are made, concatenate each seed (pseudo)contig with matched contigs to make next round
    DB = baseName + "_DB"
    matches = baseName + "_matches"
    for i in range(leng-1,-1,-1):
    #for i in [leng-1]:
        # Make DB out of fSeed, whatever it is right now
        os.system("{!s}rait -new -i {!s}-2 -o {!s}-{!s} >/dev/null 2>&1".format(raiPath, fSeed, DB, i))
        # Match ith contigs to DB
        toMatch = "{!s}_{!s}".format(baseName,i)
        os.system("{!s}rai -I {!s}-1 -d {!s}-{!s} >/dev/null 2>&1".format(raiPath, toMatch, DB, i))
        short = toMatch.rsplit("/",1)[1]
        os.system("cp {!s}/{!s}-1.bin {!s}".format(os.getcwd(), short, matches)) # moves results to results folder
        os.system("rm {!s}/{!s}-1.bin".format(os.getcwd(), short))
        
        # Construct matching dictionary
        matchDict = {}
        fMatch = open(matches,'r')
        for l in fMatch.readlines():
            [u1,u2] = l.rstrip().split(" ")
            if u2 in matches:
                matchDict[u2].append(u1)
            else:
                matchDict[u2] = [u1]
    
        # Make concatenated seeds for next DB
        l2 = open(fSeed + "-2",'w')
        for j in matchDict.keys():
            newContig = "pseudocontig_"+"{!s}".format(ct).zfill(3)
            roots.add(newContig)
            masterDict[newContig] = [j]
            fpc = open("{!s}{!s}.fna".format(genePath,newContig),'w')
            fpc.write(">{!s}\n".format(newContig))
            _, seq = readSequence("{!s}{!s}.fna".format(genePath, j))
            fpc.write(seq)
            os.system("rm {!s}{!s}.fna".format(genePath,j)) # clear up space
            for v in matchDict[j]:
                _, seq = readSequence("{!s}{!s}.fna".format(genePath, v))
                fpc.write(seq)
                os.system("rm {!s}{!s}.fna".format(genePath,v)) # clear up space
                masterDict[newContig].append(v)
                if v in roots:
                    roots.remove(v)
            fpc.write("\n")
            fpc.close()
            l2.write("{!s}\t{!s}{!s}.fna\n".format(newContig,genePath,newContig))
            ct += 1
        l2.close()


    # process results from main loop to get initial clusters
    for r in roots:
        fOut.write("{!s}: {!s}\n".format(r,getLeaves(masterDict,r)))


if __name__ == "__main__":
    main(sys.argv[1:])