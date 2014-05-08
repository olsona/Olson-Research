#! /usr/bin/env python

# imports necessary for functioning
import sys, getopt, string, re, os
import horatioConstants as hcon
import horatioUtils as hutil
import horatioFunctions as hfun
from horatioClasses import Cluster, Contig
#import cPickle as pickle

# imports necessary for debugging
#import matplotlib as mpl
#mpl.use('Agg')
import pprint

def main(argv):
    #---PREPROCESSING---#
    # initialize values
    inputFile = ''
    outputFile = ''
    cutSchedule = []
    computePath = ''
    scoreFunction = ''
    nameFile = ''
    names = []
    jtArg = ''
    joinThreshold = 0.5
    ntArg = ''
    neighborThreshold = 0.1
    matchLevel = 'genus' #***
	
    # get inputs
    try:
	opts, args = getopt.getopt(argv,"hi:o:c:p:s:f:j:n:",["ifile=","ofile=","cut=","path=","score=","namefile=","jointhreshold=","neighborthreshold="])
    except getopt.GetoptError:
	print hcon.usageString
	sys.exit(2)
    for opt, arg in opts:
	if opt == '-h': #user needs help
	    print hcon.usageString
	    print hcon.bigUsageString
	    sys.exit()
        elif opt in ("-i", "--ifile"):     # input file
	   inputFile = arg
	elif opt in ("-o", "--ofile"):     # output file
	   outputFile = arg
	elif opt in ("-c", "--cut"):       # cut schedule
	   cutSchedule = [int(n) for n in arg.lstrip()[1:-1].split(',')]
	elif opt in ("-s", "--score"):     # score function
	   scoreFunction = arg.lower()
	elif opt in ("-p", "--path"):      # compute path (only used for RAIphy scoring)
	   computePath = arg
	elif opt in ("-j", "--joiningthreshold"):  # joinThreshold, for cluster joining
	   jtArg = arg
	elif opt in ("-n", "--neighborthreshold"): # neighborThreshold, for counting a cluster as a neighbor
	   ntArg = arg
	# ***
	elif opt in ("-f", "--namefile"):  # *** FOR CORRECTNESS TESTING ONLY
	   nameFile = arg
	# ***

    # check inputs for validity
    if len(inputFile) == 0: # no input file provided
        print 'Missing argument: -i'
        print hcon.usageString
        sys.exit(2)
    elif len(outputFile) == 0:
	print 'Missing argument: -o'
	print hcon.usageString
	sys.exit(2)
    if scoreFunction not in ["raiphy","tetra","tacoa"]:
	scoreFunction = "raiphy"
    elif len(computePath) == 0 and scoreFunction == "raiphy":
	print 'Missing argument: -p'
	print hcon.usageString
	sys.exit(2)
    if len(cutSchedule) == 0:
	cutSchedule = hcon.defaultSchedule
    if jtArg:
        try:
            joinThreshold = float(jtArg)
        except ValueError:
            print "Cannot parse {!s} as a float.".format(jtArg)
            sys.exit(2)
    if ntArg:
        try:
            neighborThreshold = float(ntArg)
        except ValueError:
            print "Cannot parse {!s} as a float.".format(ntArg)
            sys.exit(2)
    # ***
    if nameFile:
        try:
            nf = open(nameFile,'r')
            for l in nf.readlines():
                n = l.rstrip()
                if len(n) > 0:
                    names.append(n)
        except IOError:
            print 'NameFile {!s} could not be open/read.'.format(nameFile)
            sys.exit(2)
        except:
            print "Unexpected error:", sys.exc_info()[0]
            sys.exit(2)
    # ***

    # properly format input metagenome file
    f = open(inputFile, 'r')
    baseName = inputFile.rsplit(".",1)[0]+"_working_"+scoreFunction
    newName = inputFile
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
	fN.close()		
    f.close()
	
    # separate out contigs by size, according to user-supplied cut schedule
    fNext = newName
    genePath = newName.rsplit("/",1)[0]+"/contigs/"
    hutil.ensureDir(genePath)
    leng = len(cutSchedule)
    for i in range(leng):
        workingFile = fNext
	thr = cutSchedule[i]*1000
	bgr = "{!s}_{!s}_next".format(baseName,i)
	if i == 0:
	    os.system("perl discardSmlr.pl {!s} {!s} {!s} {!s}".format(thr, genePath, workingFile, bgr))
	else:
	    smlr = "{!s}_{!s}".format(baseName, i)
	    os.system("perl sepSizeListDownUp.pl {!s} {!s} {!s} {!s} {!s}".format(thr, genePath, workingFile, smlr, bgr))
	fNext = bgr
	
    # make initial seed file
    fSeed = "{!s}_{!s}_seed".format(baseName, leng)
    os.system("perl processSeedFile.pl {!s} {!s} {!s}".format(genePath, fNext, fSeed))
    # initialize clusters and contigs
    allClusters = {}      # a dict matching cluster seeds to clusters themselves
    allContigs = {}       # a dict matching contig names to contigs themselves
    contigs2Clusters = {} # matches contig names to cluster names
    clusters2Contigs = {} # matches cluster seeds (aka names) to contig names
    f = open(fSeed+"-2",'r')
    for l in f.readlines():
        sp = l.rstrip().split("\t")
   	nm = sp[0]      # nm is the name of the contig
	cl = Cluster(nm) # starting a new cluster cl with nm as its seed
	co = Contig(nm)
	allClusters[nm] = cl
	allContigs[nm] = co
	contigs2Clusters[nm] = cl
	clusters2Contigs[nm] = [co]
    newContigCount = 0

    #---MAIN LOOP---#
    for i in range(leng-1, 0, -1):
        iterString = "{!s}-{!s}".format(str(cutSchedule[i-1]).zfill(2),str(cutSchedule[i]).zfill(2))
        
        # create DB and query files, apply user-supplied scoring function to them
        DB = "{!s}_{!s}_DB".format(baseName,i)
	matches = "{!s}_{!s}_matches".format(baseName,i)
	toMatch = "{!s}_{!s}".format(baseName,i)
	if scoreFunction == "tacoa":
	    hfun.scoreTACOA(DB, fSeed, matches, toMatch, allContigs)				
	elif scoreFunction == "tetra":
	    hfun.scoreTETRA(DB, fSeed, matches, toMatch, allContigs)
	elif scoreFunction == "raiphy":
	    hfun.scoreRAIphy(DB, computePath, fSeed, matches, toMatch, allContigs)
	    
	# construct matching dictionary for internal use
	matchDict = {}
	fMatching = open(matches,'r')
	lns = fMatching.readlines()
	dbNames = lns[0].rstrip().split(",")
	queryNames = lns[1].rstrip().split(",")
	for d in dbNames:
	    matchDict[d] = []
	for row in range(2,len(lns)):
	    line = lns[row]
	    bestMatch = line.rstrip().split(", ")[0].split(":")
	    bestIndex = int(bestMatch[1])
	    #bestScore = float(bestMatch[0])
	    dbItem = dbNames[bestIndex]
	    queryItem = queryNames[row-2]
	    matchDict[dbItem].append(queryItem)
	    
	fMatching.close()
	    
	for co in sorted(contigs2Clusters):
            print "{!s}: {!s}".format(co, contigs2Clusters[co].seed)
            
	# create new seeds through concatenation and prepare for next DB creation
	fSeed = "{!s}_{!s}_seed".format(baseName, i)
	l2 = open(fSeed + "-2",'w')
	# Make concatenated seeds, add to DB index list
	for seed in matchDict.keys():
	    cl = contigs2Clusters[seed]
	    newContigName = "pseudocontig_"+"{!s}".format(newContigCount).zfill(4)
	    fNewContig = open("{!s}{!s}.fna".format(genePath,newContigName),'w')
	    fNewContig.write(">{!s}\n".format(newContigName))
	    _, seq = hutil.readSequence("{!s}{!s}.fna".format(genePath, seed))
	    fNewContig.write(seq)
	    #os.system("rm {!s}{!s}.fna".format(genePath,seed)) # clear up space
	    nCo = Contig(newContigName)
	    allContigs[newContigName] = nCo
	    contigs2Clusters[newContigName] = cl
	    clusters2Contigs[cl.seed].append(nCo)
	    cl.root = newContigName
	    cl.addNode(newContigName, seed)
	    # concatenate all sequences from all matching children
	    for child in matchDict[seed]:
	        co = allContigs[child]
		cl.addNode(newContigName, child)
		# os.system("rm {!s}{!s}.fna".format(genePath,child)) # clear up space
	    fNewContig.write("\n")
	    fNewContig.close()
	    newContigCount += 1
	    l2.write("{!s}\t{!s}{!s}.fna\n".format(newContigName,genePath,newContigName))
        
        print iterString + " done"
	l2.close()
			
    #---POSTPROCESSING---#

if __name__ == "__main__":
    main(sys.argv[1:])