#! /usr/bin/env python

# imports necessary for functioning
import sys, string, re, os
import horatioConstants as hcon
import horatioUtils as hutil
import horatioFunctions as hfun
from horatioClasses import Cluster, Contig
import argparse
import cPickle as pickle
import sklearn.cluster
import numpy

# imports necessary for debugging and correctness
#import matplotlib as mpl
#mpl.use('Agg')
import horatioCorrectness as hcorr
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
	joinThreshold = 0.5
	neighborThreshold = 0.1
	splitThreshold = []
	#matchLevel = 'genus' #***
	# get inputs
	parser = argparse.ArgumentParser(description="Read in arguments for horatio.py")
	parser.add_argument("-i","--input", help="Input .fasta file", required=True)
	parser.add_argument("-o","--output", help="Prefix for all output files", required=True)
	parser.add_argument("-c","--cut", help="Cut schedule", required=True)
	parser.add_argument("-s","--score", help="Scoring function", choices=['raiphy','tetra','tacoa'], default='tetra')
	parser.add_argument("-p","--path", help="Computation path (necessary only for RAIphy scoring)")
	parser.add_argument("-n","--neighbor", help="Neighborhood threshold",type=float, default=0.01)
	parser.add_argument("-j","--join", help="Joining threshold",type=float, default=0.5)
	parser.add_argument("-a","--ap", help="AP preference", choices=['min','median','mean','max','40','60','70','80','90','len','len2'],default="max")
	parser.add_argument("-d","--doAP", help="Do AP or not", type=int, choices=[0,1],default=0)
	parser.add_argument("-k","--clusterLimit", help="Minimum size for a cluster to be included in AP", type=int, default=2)
	parser.add_argument("-l","--split", help="Split threshold")
	parser.add_argument("-f","--names", help="Name file")
	args = parser.parse_args()
	inputFile = args.input
	outputFile = args.output
	scoreFunction = args.score
	prefFun = args.ap
	doAP = args.doAP
	k = args.clusterLimit
	cutSchedule = [int(n) for n in args.cut.lstrip()[1:-1].split(',')]
	if args.split:
		splitThreshold = [float(n) for n in args.split.lstrip()[1:-1].split(',')]
	else:
		splitThreshold = [-1000.0]*(len(cutSchedule)-1)
	computePath = args.path
	joinThreshold = args.join
	neighborThreshold = args.neighbor
	if args.names:
		nameFile = args.names

	# check inputs for validity
	if scoreFunction == "raiphy" and computePath is None:
		print 'Missing argument: -p'
		print hcon.usageString
		sys.exit(2)
	if len(splitThreshold) != len(cutSchedule) -1:
		print "Incorrect number of split thresholds."
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

	#print "Starting"

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
	
	#print "Separating"
	# separate out contigs by size, according to user-supplied cut schedule
	fNext = newName
	#print fNext
	genePath = newName.rsplit("/",1)[0]+"/contigs/"
	hutil.ensureDir(genePath)
	leng = len(cutSchedule)
	for i in range(leng):
		workingFile = fNext
		#print i, workingFile
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
	
	#print "Seeding"
	# make initial seed file
	fSeed = "{!s}_{!s}_seed".format(baseName, leng)
	os.system("perl processSeedFile.pl {!s} {!s} {!s}".format(genePath, fNext, fSeed))
	# initialize clusters and contigs
	allClusters = {}	  # a dict matching cluster seeds to clusters themselves
	allContigs = {}		  # a dict matching contig names to contigs themselves
	contigs2Clusters = {} # matches contig names to clusters
	clusters2Contigs = {} # matches cluster seeds (aka names) to a list of contigs
	f = open(fSeed+"-2",'r')
	#os.system("head -1 {!s}-2".format(fSeed))
	for l in f.readlines():
		sp = l.rstrip().split("\t")
		nm = sp[0]		  # nm is the name of the contig
		cl = Cluster(nm)  # starting a new cluster cl with nm as its seed
		co = Contig(nm)
		allClusters[nm] = cl
		allContigs[nm] = co
		contigs2Clusters[nm] = cl
		clusters2Contigs[nm] = [co]
	newContigCount = 0

	# ***
	#rightDists = {"{!s}-{!s}".format(str(cutSchedule[i]).zfill(2),str(cutSchedule[i+1]).zfill(2)):[] for i in range(leng-1)}
	#wrongDists = {"{!s}-{!s}".format(str(cutSchedule[i]).zfill(2),str(cutSchedule[i+1]).zfill(2)):[] for i in range(leng-1)}
	#distLog = open("{!s}_distLog".format(outputFile),"w")
	# ***

	numLeaves = 0

	#---MAIN LOOP---#
	for i in range(leng-1, 0, -1):
		#print i
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
		#print fMatching
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
			#print "{!s} -> {!s}, {:01.3f}".format(queryItem,dbItem,bestScore)
			if bestScore < splitThreshold[i-1]:
				newSeeds.append(queryItem)
			else:
				if dbItem in matchDict:
					matchDict[dbItem].append([queryItem,bestScore])
				else:
					matchDict[dbItem] = [[queryItem,bestScore]]
			# check for other good matches to queryItem
			myThresh = 0.0
			queryContig = allContigs[queryItem]
			# matchCluster = contigs2Clusters[dbItem]  # ***
			if bestScore >= splitThreshold[i-1]:
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
			#clName = contigs2Clusters[dbItem].seed
			#isCorrect = hcorr.checkCorrectGenusOlsonFormat(clName,queryItem)
			#if isCorrect:
			#	rightDists[iterString].append(bestScore)
			#else:
			#	wrongDists[iterString].append(bestScore)
			# ***
		fMatching.close()
		#pprint.pprint(matchDict)
			
		# create new seeds through concatenation and prepare for next DB creation
		fSeed = "{!s}_{!s}_seed".format(baseName, i)
		l2 = open(fSeed + "-2",'w')
		# Make concatenated seeds
		for seed in matchDict.keys():
			#print seed, matchDict[seed]
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
			cl.addNode(newContigName, seed, 0.0)
			# concatenate all sequences from all matching children
			for [child,score] in matchDict[seed]:
				numLeaves += 1
				#print child
				co = allContigs[child]
				#print child, co.goodMatches
				cl.addNode(newContigName, child, score)
				_, seq = hutil.readSequence("{!s}{!s}.fna".format(genePath,child))
				fNewContig.write(seq)
				# add neighbors
				for m in co.goodMatches:
					cl.addMatch(m)
				# os.system("rm {!s}{!s}.fna".format(genePath,child)) # clear up space
			#print cl.root	
			#pprint.pprint(cl.dict)
			#print cl.getLeaves()
			fNewContig.write("\n")
			fNewContig.close()
			newContigCount += 1
			#l2.write("{!s}\t{!s}{!s}.fna\n".format(newContigName,genePath,newContigName))

		# add split seeds to DB
		for nSeed in newSeeds:
			#print "{!s} is a new seed".format(nSeed)
			co = allContigs[nSeed]
			cl = Cluster(nSeed)
			allClusters[nSeed] = cl
			contigs2Clusters[nSeed] = cl
			clusters2Contigs[nSeed] = [co]
			numLeaves += 1
			#l2.write("{!s}\t{!s}{!s}.fna\n".format(nSeed, genePath, nSeed))
			
		# go through clusters again, evaluate what clusters should be joined
		# identify edges (make simplifying assumption that if A is a neighbor of B, then B is a neighbor of A)
		graph = {}
		for clID in allClusters:
			clust = allClusters[clID]
			ratio, nameB = clust.getMostCommonNeighbor()
			if ratio >= joinThreshold:
				#print "Edge between {!s} and {!s}".format(clID,nameB)
				if clID in graph:
					graph[clID].add(nameB)
				else:
					graph[clID] = set([nameB])
				if nameB in graph:
					graph[nameB].add(clID)
				else:
					graph[nameB] = set([clID])
		#print "Graph"			 
		#pprint.pprint(graph)
		#print
		
		# DFS on clusters to find all connected components
		partition = []
		metaVisited = set()
		for root in graph:
			if root not in metaVisited:
				component = hutil.dfs(graph,root)
				if len(component) > 1:
					partition.append(component)
				metaVisited = metaVisited | component
										  
		# iterate through partition of clusters and merge as appropriate
		for p in partition:
			pList = list(p)
			mainClID = pList[0]
			mainClust = allClusters[mainClID]
			restClust = [allClusters[ID] for ID in pList[1:]]
			# make ubercontig
			newContigName = "pseudocontig_"+"{!s}".format(newContigCount).zfill(4)
			newContig = Contig(newContigName)
			allContigs[newContigName] = newContig
			contigs2Clusters[newContigName] = mainClust
			clusters2Contigs[mainClID].append(newContig)
			fNewContig = open("{!s}{!s}.fna".format(genePath,newContigName),'w')
			fNewContig.write(">{!s}\n".format(newContigName))
			_, seq = hutil.readSequence("{!s}{!s}.fna".format(genePath, mainClID))
			fNewContig.write(seq)
			for rCl in restClust:
				co = allContigs[rCl.root]
				_, seq = hutil.readSequence("{!s}{!s}.fna".format(genePath,co.name))
				fNewContig.write(seq)
				for m in co.goodMatches:
					mainClust.addMatch(m)	  # it's entirely possible that the root is not something I made above in the initial matching loop, and so would have neighbors
				# os.system("rm {!s}{!s}.fna".format(genePath,child)) # clear up space
			fNewContig.write("\n")
			fNewContig.close()
			newContigCount += 1
			# add clusters
			#print mainClust.getLeaves()
			mainClust.addClusters(restClust, newContigName)
			#print mainClust.getLeaves()
			# remove restClust from allClusters, update all entries in clusters2Contigs and contigs2Clusters
			for rCl in restClust:
				contigNames = clusters2Contigs[rCl.seed]
				for con in contigNames:
					clusters2Contigs[mainClust.seed].append(con)
					contigs2Clusters[con] = mainClust
				clusters2Contigs.pop(rCl.seed)
				allClusters.pop(rCl.seed)
				#print "Popped {!s}".format(rCl.seed)
			#print sorted(mainClust.getLeaves())
			#print ", ".join([c for c in sorted(allClusters.keys())])
	
		# finally write root contigs to db index file
		for clID in allClusters:
			clust = allClusters[clID]
			root = clust.root
			leaves = clust.getLeaves()
			
			l2.write("{!s}\t{!s}{!s}.fna\n".format(root, genePath, root))
			clust.closeDict = {}
		
		#for cl in sorted(allClusters.keys()):
		#	 print cl
	 
		# *** compute correctness distributions
		#rdata = rightDists[iterString]
		#wdata = wrongDists[iterString]
		#distLog.write(iterString + "\n")
		#if rdata:
		#	distLog.write("Right ({!s}): {:03.2f} - {:03.2f}\n".format(len(rdata),min(rdata),max(rdata)))
		#if wdata:
		#	distLog.write("Wrong ({!s}): {:03.2f} - {:03.2f}\n".format(len(wdata),min(wdata),max(wdata)))
		#distLog.write("\n")
		#hcorr.comparisonPlot(rdata, wdata, iterString, outputFile, "seed", "Correct distances", "Incorrect Distances")
		# ***	
		
		l2.close()
		#print i, len(allContigs)
			
	#---POSTPROCESSING---#
	
	
	#print allClusters.keys()
	actualClusterList1 = open("{!s}_actualclusters-1".format(outputFile),'w')
	actualClusterList2 = open("{!s}_actualclusters-2".format(outputFile),'w')
	clusterLengths = []
	totalLength = 0
	clusters = []
	numLeaves = 0
	for c in allClusters:
		r = allClusters[c].root
		l = allClusters[c].getLeaves()
		numLeaves = len(l)
		if numLeaves >= k:
			actualClusterList2.write("{!s}\t{!s}/{!s}.fna\n".format(r,genePath,r))
			actualClusterList1.write("{!s}/{!s}.fna\n".format(genePath,r))
			clusters.append(allClusters[c])
			sz = os.path.getsize("{!s}/{!s}.fna".format(genePath,r))
			clusterLengths.append(sz)
	actualClusterList1.close()
	actualClusterList2.close()
	print "{!s},{!s}".format(len(allClusters),len(clusters))
	#print numLeaves
	fSeed = "{!s}_actualclusters".format(outputFile)

	#final distances
	DB = "{!s}_final_DB".format(baseName)
	if scoreFunction == "tacoa":
		hfun.scoreTACOAFinal(DB, fSeed, outputFile)			
	elif scoreFunction == "tetra":
		hfun.scoreTETRAFinal(DB, fSeed, outputFile)
	elif scoreFunction == "raiphy":
		hfun.scoreRAIphyFinal(DB, fSeed, computePath, outputFile)
		
	finalDists = hutil.makeDistanceMatrix("{!s}_dists_sorted".format(outputFile))
	if doAP:	
		#os.system("rm {!s}_dists_sorted".format(outputFile))
		os.system("rm {!s}_actualclusters*".format(outputFile))
		if prefFun == "len":
			minDist = numpy.min(finalDists)
			maxDist = numpy.max(finalDists)
			minLen = numpy.min(clusterLengths)
			maxLen = numpy.max(clusterLengths)
			m = (maxDist - minDist)/(maxLen - minLen)
			b = minDist - m*minLen
			pref = [m*c + b for c in clusterLengths]
		elif prefFun == "len2":
			minDist = numpy.min(finalDists)
			maxDist = numpy.max(finalDists)
			minLen = numpy.min(clusterLengths)**2
			maxLen = numpy.max(clusterLengths)**2
			m = m = (maxDist - minDist)/(maxLen - minLen)
			b = minDist - m*minLen
			pref = [m*(c**2) + b for c in clusterLengths]
		else:
			pref = hcon.apPreferences[prefFun](finalDists)
		_, labels = sklearn.cluster.affinity_propagation(finalDists,preference=pref)
		metaClustering = hutil.processAPLabels(labels, clusters)
		#print metaClustering
		# iterate through partition of clusters and merge as appropriate
		finalClusters = []
		fOutC = open("{!s}_clusters".format(outputFile),'w')
		for p in metaClustering:
			pList = list(p)
			mainClust = pList[0]
			mainClID = mainClust.seed
			restClust = pList[1:]
			# make ubercontig
			newContigName = "pseudocontig_"+"{!s}".format(newContigCount).zfill(4)
			newContig = Contig(newContigName)
			allContigs[newContigName] = newContig
			contigs2Clusters[newContigName] = mainClust
			clusters2Contigs[mainClID].append(newContig)
			fNewContig = open("{!s}{!s}.fna".format(genePath,newContigName),'w')
			fNewContig.write(">{!s}\n".format(newContigName))
			_, seq = hutil.readSequence("{!s}{!s}.fna".format(genePath, mainClID))
			fNewContig.write(seq)
			for rCl in restClust:
				co = allContigs[rCl.root]
				_, seq = hutil.readSequence("{!s}{!s}.fna".format(genePath,co.name))
				fNewContig.write(seq)
				# os.system("rm {!s}{!s}.fna".format(genePath,child)) # clear up space
			fNewContig.write("\n")
			fNewContig.close()
			#newContigCount += 1
			# add clusters
			mainClust.addClusters(restClust, newContigName)
			finalClusters.append(mainClust)
			leaves = mainClust.getLeaves()
			fOutC.write(str(leaves) + "\n")
			# remove restClust from allClusters, update all entries in clusters2Contigs and contigs2Clusters
			for rCl in restClust:
				contigNames = clusters2Contigs[rCl.seed]
				for con in contigNames:
					clusters2Contigs[mainClust.seed].append(con)
					contigs2Clusters[con] = mainClust
				clusters2Contigs.pop(rCl.seed)
				allClusters.pop(rCl.seed)
		pickle.dump(finalClusters,open("{!s}_clusters_pickle".format(outputFile),"wb")) 
	else:
		clusters = []
		for c in allClusters:
			r = allClusters[c].root
			l = allClusters[c].getLeaves()
			clusters.append(allClusters[c])
		pickle.dump(clusters, open("{!s}_clusters_pickle".format(outputFile),'wb'))
	
	# Get rid of files we're not using any more
	os.system("rm -r {!s} >/dev/null 2>&1".format(genePath))
	os.system("rm {!s} >/dev/null 2>&1".format(DB))
	os.system("rm {!s} >/dev/null 2>&1".format(toMatch))
	os.system("rm {!s} >/dev/null 2>&1".format(fSeed))
	for i in range(leng+1):
		os.system("rm {!s}_{!s}* >/dev/null 2>&1".format(baseName,i))
	os.system("rm {!s}_dists_sorted >/dev/null 2>&1".format(outputFile))
	os.system("rm {!s}_actualclusters*".format(outputFile))
	
	#distLog.close()


if __name__ == "__main__":
	main(sys.argv[1:])