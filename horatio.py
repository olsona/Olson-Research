#! /usr/bin/env python

# imports necessary for functioning
import sys, getopt, string, re, os
import horatioConstants as hcon
import horatioUtils as hutil
import horatioFunctions as hfun
from horatioClasses import Cluster, Contig
#import cPickle as pickle

# imports necessary for debugging and correctness
import matplotlib as mpl
mpl.use('Agg')
import horatioCorrectness as hcorr
#import pprint

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
    splitThreshold = []
    stArg = ''	
    #matchLevel = 'genus' #***
    # get inputs
    try:
	opts, args = getopt.getopt(argv,"hi:o:c:p:s:f:j:n:l:",["ifile=","ofile=",\
	   "cut=","path=","score=","namefile=","jointhreshold=",\
	   "neighborthreshold=","splitthreshold="])
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
	elif opt in ("-l", "--splitthreshold:"):   # splitThreshold, for starting a new cluster
	    stArg = arg
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
    if stArg:
        try:
            splitThreshold = [float(n) for n in arg.lstrip()[1:-1].split(',')]
            if len(splitThreshold) != len(cutSchedule)-1:
                print "Incorrect number of thresholds."
                sys.exit(2)
        except ValueError, IOError:
            print "Cannot parse {!s} as a list of floats.".format(stArg)
            sys.exit(2)
    else:
        splitThreshold = [-1000.0] * (len(cutSchedule)-1)
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
	    os.system("perl discardSmlr.pl {!s} {!s} {!s} {!s}".\
	       format(thr, genePath, workingFile, bgr))
	else:
	    smlr = "{!s}_{!s}".format(baseName, i)
	    os.system("perl sepSizeListDownUp.pl {!s} {!s} {!s} {!s} {!s}".\
	       format(thr, genePath, workingFile, smlr, bgr))
	fNext = bgr
	
    # make initial seed file
    fSeed = "{!s}_{!s}_seed".format(baseName, leng)
    os.system("perl processSeedFile.pl {!s} {!s} {!s}".format(genePath, fNext, fSeed))
    # initialize clusters and contigs
    allClusters = {}      # a dict matching cluster seeds to clusters themselves
    allContigs = {}       # a dict matching contig names to contigs themselves
    contigs2Clusters = {} # matches contig names to clusters
    clusters2Contigs = {} # matches cluster seeds (aka names) to a list of contigs
    f = open(fSeed+"-2",'r')
    for l in f.readlines():
        sp = l.rstrip().split("\t")
   	nm = sp[0]        # nm is the name of the contig
	cl = Cluster(nm)  # starting a new cluster cl with nm as its seed
	co = Contig(nm)
	allClusters[nm] = cl
	allContigs[nm] = co
	contigs2Clusters[nm] = cl
	clusters2Contigs[nm] = [co]
    newContigCount = 0

    # ***
    rightDists = {"{!s}-{!s}".format(str(cutSchedule[i]).zfill(2),str(cutSchedule[i+1]).zfill(2)):[] for i in range(leng-1)}
    wrongDists = {"{!s}-{!s}".format(str(cutSchedule[i]).zfill(2),str(cutSchedule[i+1]).zfill(2)):[] for i in range(leng-1)}

    #---MAIN LOOP---#
    for i in range(leng-1, 0, -1):
        iterString = "{!s}-{!s}".format(str(cutSchedule[i-1]).zfill(2),\
            str(cutSchedule[i]).zfill(2))
        
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
	newSeeds = []
	for row in range(2,len(lns)):
	    line = lns[row]
	    bestMatch = line.rstrip().split(", ")[0].split(":")
	    bestIndex = int(bestMatch[1])
	    bestScore = float(bestMatch[0])
	    dbItem = dbNames[bestIndex]
	    queryItem = queryNames[row-2]
	    if bestScore < splitThreshold[i-1]:
	        newSeeds.append(queryItem)
	    else:
	       if dbItem in matchDict:
	           matchDict[dbItem].append(queryItem)
	       else:
	           matchDict[dbItem] = [queryItem]
	    # check for other good matches to queryItem
	    myThresh = 0.0
	    queryContig = allContigs[queryItem]
	    # matchCluster = contigs2Clusters[dbItem]  # ***
	    if bestScore > 0:
	        myThresh = (1.0-neighborThreshold)*bestScore
	    else:
	        myThresh = (1.0+neighborThreshold)*bestScore
	    # check for neighbors
	    for l in line.rstrip().split(", ")[1:]:
	        entry = l.split(":")
	        score = float(entry[0])
	        if score >= myThresh:
	            neighborInd = int(entry[1])
	            neighborName = dbNames[neighborInd]
	            neighborCluster = contigs2Clusters[neighborName]
	            queryContig.goodMatches.append([neighborCluster.seed,score])
	    # ***
	    #print "{!s} matched: {!s}\n\tNear: {!s}\n".format(queryContig,\
	    #   matchCluster.seed, [m[0] for m in queryContig.goodMatches])
	    # ***
	            
	    # ***
            clName = contigs2Clusters[dbItem].seed
            isCorrect = hcorr.checkCorrectGenusOlsonFormat(clName,queryItem)
            if isCorrect:
                rightDists[iterString].append(bestScore)
            else:
                wrongDists[iterString].append(bestScore)
            # ***
	
	fMatching.close()
            
	# create new seeds through concatenation and prepare for next DB creation
	fSeed = "{!s}_{!s}_seed".format(baseName, i)
	l2 = open(fSeed + "-2",'w')
	# Make concatenated seeds
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
		_, seq = hutil.readSequence("{!s}{!s}.fna".format(genePath,child))
		fNewContig.write(seq)
		# add neighbors
		for m in co.goodMatches:
		    cl.addMatch(m)
		# os.system("rm {!s}{!s}.fna".format(genePath,child)) # clear up space
	    fNewContig.write("\n")
	    fNewContig.close()
	    newContigCount += 1
	    l2.write("{!s}\t{!s}{!s}.fna\n".format(newContigName,genePath,newContigName))
	# add split seeds to DB
	for nSeed in newSeeds:
	    print nSeed
	    co = allContigs[nSeed]
	    cl = Cluster(nSeed)
	    allClusters[nSeed] = cl
	    contigs2Clusters[nSeed] = cl
	    clusters2Contigs[nSeed] = [co]
	    l2.write("{!s}\t{!s}{!s}.fna\n".format(nSeed, genePath, nSeed))
	    
	# go through clusters again, evaluate what clusters should be joined
	# identify edges (make simplifying assumption that if A is a neighbor of B, then B is a neighbor of A)
	partition = []
	sortedClusterIDs = sorted(allClusters.keys())
	for a in range(len(sortedClusterIDs)):
	    nameA = sortedClusterIDs[a]
	    clustA = allClusters[sortedClusterIDs[a]]
	    ratio, nameB = clustA.getMostCommonNeighbor()
	    if ratio >= joinThreshold:
	        print nameA, nameB, ratio
	        indB = sortedClusterIDs.index(nameB)
	        if indB > a:
	            print "Joining {!s} and {!s}".format(nameA,nameB)
	            # search through partition to see if nameA or nameB are already there
	            for p in partition:
	                if nameA in p or nameB in p:
	                    p.add(nameA)
	                    p.add(nameB)
	                    break
	            newSet = set([nameA,nameB])
	            partition.append(newSet)
	print "\nPartition:"
	for p in partition:
	    print "\t* "+", ".join(c for c in p)
	print     
	                                           
	# iterate through partition of clusters and merge as appropriate
	#for p in partition:
	#    pList = list(p)
	#    mainClID = pList[0]
	#    mainClust = allClusters[mainClID]
	#    restClust = [allClusters[ID] for ID in pList[1:]]
	#    # TODO: make ubercontig
	#    newContigName = "pseudocontig_"+"{!s}".format(newContigCount).zfill(4)
	#    fNewContig = open("{!s}{!s}.fna".format(genePath,newContigName),'w')
	#    fNewContig.write(">{!s}\n".format(newContigName))
	#    _, seq = hutil.readSequence("{!s}{!s}.fna".format(genePath, mainClID))
	#    fNewContig.write(seq)
	#    for rCl in restClust:
	#        co = allContigs[rCl.root]
	#	_, seq = hutil.readSequence("{!s}{!s}.fna".format(genePath,co))
	#	fNewContig.write(seq)
	#	for m in co.goodMatches:
	#	    mainClust.add(m)    # it's entirely possible that the root is not something I made above in the initial matching loop, and so would have neighbors
	#	# os.system("rm {!s}{!s}.fna".format(genePath,child)) # clear up space
	#    fNewContig.write("\n")
	#    fNewContig.close()
	#    newContigCount += 1
	#    # add clusters
	#    mainClust.addClusters(restClust, newContigName)
	#    # remove restClust from allClusters, update all entries in clusters2Contigs and contigs2Clusters
	#    for rCl in restClust:
	#        contigNames = clusters2Contigs[rCl.seed]
	#        for con in contigNames:
	#            clusters2Contigs[mainClust.seed].append(con)
	#            contigs2Clusters[con] = mainClust
	#        clusters2Contigs.pop(rCl.seed)
	#        allClusters.pop(rCl)
	#
	## finally write root contigs to db index file
	#for clID in allClusters:
	#    clust = allClusters[clID]
	#    root = clust.root
	#    l2.write("{!s}\t{!s}{!s}.fna\n".format(root, genePath, root))
	 
	# *** compute correctness distributions
	#rdata = rightDists[iterString]
	#wdata = wrongDists[iterString]
	#if rdata:
	#    print "Right ({!s}): {:03.2f}-{:03.2f}".format(len(rdata),min(rdata),max(rdata))
	#if wdata:
	#    print "Wrong ({!s}): {:03.2f}-{:03.2f}".format(len(wdata),min(wdata),max(wdata))
	#hcorr.comparisonPlot(rdata, wdata, iterString, outputFile, "seed", "Correct distances", "Incorrect Distances")
	# ***   
        
        print iterString + " done"
	l2.close()
			
    #---POSTPROCESSING---#

if __name__ == "__main__":
    main(sys.argv[1:])