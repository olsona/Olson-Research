#! /usr/bin/env python

'''bootstrap.py - wrapper class for my MS project.'''

import sys, getopt, string, os, re, pprint, numpy
import cPickle as pickle
# http://stackoverflow.com/questions/2801882/generating-a-png-with-matplotlib-when-display-is-undefined
import matplotlib
matplotlib.use('Agg')
from bootstrapConstants import *
from bootstrapUtils import *
from bootstrapClasses import *
from bootstrapFunctions import *
from bootstrapCorrectness import *

def main(argv):
	# get inputs, check validity
	inputFile = ''
	outputFile = ''
	coolingSchedule = []
	computePath = ''
	matchLevel = ''
	scoreFunction = ''
	nameFile = ''
	names = []
	# Command line arguments
	try:
		opts, args = getopt.getopt(argv,"hi:o:c:p:s:n:",["ifile=","ofile=","cut=","path=","score=","namefile="])
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
		elif opt in ("-c", "--cut"):
			coolingSchedule = [int(n) for n in arg.lstrip()[1:-1].split(',')]
		elif opt in ("-p", "--path"):
			computePath = arg
		elif opt in ("-s", "--score"):
			scoreFunction = arg.lower()
		elif opt in ("-n", "--namefile"):
			nameFile = arg
			nf = open(nameFile,'r')
			names = [l.rstrip() for l in nf.readlines()]
	if len(inputFile) == 0:
		print 'Missing argument: -i'
		print usageString
		sys.exit(2)
	elif len(outputFile) == 0:
		print 'Missing argument: -o'
		print usageString
		sys.exit(2)
	if scoreFunction not in ["raiphy","tetra","tacoa"]:
		scoreFunction = "raiphy"
	elif len(computePath) == 0 and scoreFunction == "raiphy":
		print 'Missing argument: -p'
		print usageString
		sys.exit(2)
	if len(coolingSchedule) == 0:
		coolingSchedule = defaultSchedule
	
	matchLevel = 'genus'

	# Checking validity of inputs
	if scoreFunction == "raiphy":
		if computePath[-1] != "/":
			computePath = computePath + "/"
		if scoreFunction == "raiphy":
			try:
				temp = open(computePath+"rait",'r')
				temp.close()
			except IOError:
				print computePath+"rait cannot be opened."
				sys.exit(1)
	try:
		temp = open(inputFile, 'r')			
		temp.close()
	except IOError:
		print inputFile + "cannot be opened."			
		sys.exit(1)
	try:
		temp = open(outputFile+"_clusters", 'w')
		temp.close()
	except IOError:
		print outputFile + "cannot be opened."
		sys.exit(1)

	# properly format input file
	f = open(inputFile, 'r')
	baseName = inputFile.rsplit(".",1)[0]+"_working"
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

	# separate out files by size, using sepSizeListDownUp.pl
	fNext = newName
	genePath = newName.rsplit("/",1)[0]+"/contigs/"
	ensureDir(genePath)
	leng = len(coolingSchedule)
	for i in range(leng):
		workingFile = fNext
		thr = coolingSchedule[i]*1000
		bgr = "{!s}_{!s}_next".format(baseName,i)
		if i == 0:
			os.system("perl discardSmlr.pl {!s} {!s} {!s} {!s}".format(thr, genePath, workingFile, bgr))
		else:
			smlr = "{!s}_{!s}".format(baseName, i)
			os.system("perl sepSizeListDownUp.pl {!s} {!s} {!s} {!s} {!s}".format(thr, genePath, workingFile, smlr, bgr))
		fNext = bgr

	# Make initial seed file
	fSeed = "{!s}_{!s}_seed".format(baseName, leng)
	os.system("perl processSeedFile.pl {!s} {!s} {!s}".format(genePath, fNext, fSeed))
	# initialize clusters and contigs
	allClusters = {}
	allContigs = {}
	f = open(fSeed+"-2",'r')
	for l in f.readlines():
		sp = l.rstrip().split("\t")
		nm = sp[0]
		cl = Cluster(nm)
		allClusters[nm] = cl
		co = Contig(nm,myCluster=cl)
		allContigs[nm] = co
	ct = 0
	
	# ***
	rightDistsSeed = {"{!s}-{!s}".format(str(coolingSchedule[i]).zfill(2),str(coolingSchedule[i+1]).zfill(2)):[] for i in range(leng-1)}
	wrongDistsSeed = {"{!s}-{!s}".format(str(coolingSchedule[i]).zfill(2),str(coolingSchedule[i+1]).zfill(2)):[] for i in range(leng-1)}
	
	rightNeighborsDistsSeed = {"{!s}-{!s}".format(str(coolingSchedule[i]).zfill(2),str(coolingSchedule[i+1]).zfill(2)):[] for i in range(leng-1)}
	wrongNeighborsDistsSeed = {"{!s}-{!s}".format(str(coolingSchedule[i]).zfill(2),str(coolingSchedule[i+1]).zfill(2)):[] for i in range(leng-1)}
	
	rightDistsMax = {"{!s}-{!s}".format(str(coolingSchedule[i]).zfill(2),str(coolingSchedule[i+1]).zfill(2)):[] for i in range(leng-1)}
	wrongDistsMax = {"{!s}-{!s}".format(str(coolingSchedule[i]).zfill(2),str(coolingSchedule[i+1]).zfill(2)):[] for i in range(leng-1)}
	
	rightNeighborsDistsMax = {"{!s}-{!s}".format(str(coolingSchedule[i]).zfill(2),str(coolingSchedule[i+1]).zfill(2)):[] for i in range(leng-1)}
	wrongNeighborsDistsMax = {"{!s}-{!s}".format(str(coolingSchedule[i]).zfill(2),str(coolingSchedule[i+1]).zfill(2)):[] for i in range(leng-1)}
	
	log = open("{!s}_log".format(outputFile),'w')
	neighborLog = open("{!s}_neighborlog".format(outputFile),'w')
	# ***
	
	close = closeThreshold
	
	# Main loop: iterate through cooling schedule, creating databases, making matches, and once matches are made, concatenate each seed (pseudo)contig with matched contigs to make next round
	for i in range(leng-1,0,-1):
		iterString = "{!s}-{!s}".format(str(coolingSchedule[i-1]).zfill(2),str(coolingSchedule[i]).zfill(2))
		
		DB = "{!s}_{!s}_DB".format(baseName,i)
		matches = "{!s}_{!s}_matches".format(baseName,i)
		toMatch = "{!s}_{!s}".format(baseName,i)
		if scoreFunction == "raiphy":
			scoreRAIphy(DB, computePath, fSeed, matches, toMatch, allContigs)				
		else:
			scoringMethod[scoreFunction](DB, fSeed, matches, toMatch, allContigs)
		
		# Construct matching dictionary for internal use
		matchDict = {}
		fMatch = open(matches,'r')
		lns = fMatch.readlines()
		dbNames = lns[0].rstrip().split(",")
		contigNames = lns[1].rstrip().split(",")
		for d in dbNames:
			matchDict[d] = []
		for row in range(2,len(lns)):
			line = lns[row]
			bestMatch = line.rstrip().split(", ")[0].split(":")
			bestIndex = int(bestMatch[1])
			bestScore = float(bestMatch[0])
			parent = dbNames[bestIndex]
			child = contigNames[row-2]
			co = allContigs[child]
			matchDict[parent].append(child)
			mythresh = 0.0
			if bestScore > 0:
				mythresh = (1.0-close)*bestScore
			else:
				mythresh = (1.0+close)*bestScore
			for l in line.rstrip().split(", ")[0:]:
				entry = l.split(":")
				score = float(entry[0])
				if score >= mythresh:
					ind = int(entry[1])
					name = dbNames[ind]
					mco = allContigs[name]
					mcl = mco.myCluster
					clnm = mcl.seed
					co.goodMatches.append([clnm,score])
					# *** check quality of neighbors 
					corrSeed = correctnessDictSeed[matchLevel](clnm, child)
					if corrSeed == 1:
						rightNeighborsDistsSeed[iterString].append(score)
					else:
						wrongNeighborsDistsSeed[iterString].append(score)
					corrMax = correctnessDictMax[matchLevel](mcl, child, names)
					if corrMax == 1:
						rightNeighborsDistsMax[iterString].append(score)
					else:
						wrongNeighborsDistsMax[iterString].append(score)
					# *** 
					
			# *** check correctness of match
			cl = allContigs[parent].myCluster
			correctS = correctnessDictSeed[matchLevel](cl.seed, child)
			if correctS == 1:
				rightDistsSeed[iterString].append(bestScore)
			else:
				wrongDistsSeed[iterString].append(bestScore)
				log.write("Wrong match (seed): {!s} to {!s}\n".format(child, cl.seed))
				log.write("Original score: {!s}\n".format(bestScore))
				log.write("Matches within {!s}%: {!s}\n\n".format(close*100,co.goodMatches))
				cl = allContigs[parent].myCluster
				
			correctM = correctnessDictMax[matchLevel](cl, child, names)
			if correctM == 1:
				rightDistsMax[iterString].append(bestScore)
			else:
				wrongDistsMax[iterString].append(bestScore)
				log.write("Wrong match (max): {!s} to {!s}\n".format(child, cl.seed))
				log.write("Original score: {!s}\n".format(bestScore))
				log.write("Matches within {!s}%: {!s}\n\n".format(close*100,co.goodMatches))
			# ***

		# *** compute correctness distributions
		rdata = rightDistsSeed[iterString]
		wdata = wrongDistsSeed[iterString]
		comparisonPlot(rdata, wdata, iterString, outputFile, "_seed", "Correct distances", "Incorrect Distances")
		
		rdata = rightDistsMax[iterString]
		wdata = wrongDistsMax[iterString]
		comparisonPlot(rdata, wdata, iterString, outputFile, "_max", "Correct distances", "Incorrect Distances")
		
		rdata = rightNeighborsDistsSeed[iterString]
		wdata = wrongNeighborsDistsSeed[iterString]
		comparisonPlot(rdata, wdata, iterString, outputFile, "_neighbors_seed", "Distances between correct neighbors", "Distances between incorrect neighbors")

		rdata = rightNeighborsDistsMax[iterString]
		wdata = wrongNeighborsDistsMax[iterString]
		comparisonPlot(rdata, wdata, iterString, outputFile, "_neighbors_max", "Distances between correct neighbors", "Distances between incorrect neighbors")
		# ***
		
		fMatch.close()

		# Prepare for next DB creation
		fSeed = "{!s}_{!s}_seed".format(baseName, i)
		l2 = open(fSeed + "-2",'w')
		mergeLogSeed = []
		mergeLogMax = []
		# Make concatenated seeds
		for j in matchDict.keys():
			#print j
			cl = allContigs[j].myCluster
			newContig = "pseudocontig_"+"{!s}".format(ct).zfill(3)
			fpc = open("{!s}{!s}.fna".format(genePath,newContig),'w')
			fpc.write(">{!s}\n".format(newContig))
			_, seq = readSequence("{!s}{!s}.fna".format(genePath, j))
			fpc.write(seq)
			#os.system("rm {!s}{!s}.fna".format(genePath,j)) # clear up space
			nCo = Contig(newContig, myCluster = cl)
			cl.root = newContig
			cl.addNode(newContig, j)
			for v in matchDict[j]:
				cl.addNode(newContig, v)
				co = allContigs[v]
				for m in co.goodMatches:
					cl.addMatch(m)
					#print "Match added"
				_, seq = readSequence("{!s}{!s}.fna".format(genePath, v))
				fpc.write(seq)
				#os.system("rm {!s}{!s}.fna".format(genePath,v)) # clear up space
			allContigs[newContig] = nCo
			fpc.write("\n")
			fpc.close()
			l2.write("{!s}\t{!s}{!s}.fna\n".format(newContig,genePath,newContig))
			ct += 1
			# get info on cluster closeness
			goodMatchList, total = cl.getMatches(names)
			bestMatch = goodMatchList[0]
			mergeLogSeed.append([cl.seed, bestMatch[0], float(bestMatch[1])/float(total)])
			_, mymax = cl.purityMax(names)
			mergeLogMax.append([mymax, bestMatch[0], float(bestMatch[1])/float(total)])
			
		print iterString + " done"
		l2.close()
		
		neighborLog.write(iterstring+"\n")
		for a in range(len(mergeLogSeed)):
			neighborLog.write("{!s}, {:03.2f}%:\t{!s}\t{!s}\n".format(mergeLogMax[a][1], mergeLogMax[a][2]*100.0, mergeLogMax[a][0], mergeLogSeed[a][0]))
		neighborLog.write("\n")

	log.close()
	neighborLog.close()

	# process results from main loop to get clusters and distances
	totalCluster = {}
	fOutC = open("{!s}_clusters".format(outputFile),'w')
	for c in allClusters:
		#print "-----"
		r = allClusters[c].root
		l = allClusters[c].getLeaves()
		cl = [li for li in l]
		cl.append(r)
		fOutC.write("{!s}\n".format(cl))
		totalCluster[r] = allClusters[c]
	fOutC.close()
	
	# get right/wrong distance distributions ***
	dists={"right":rightDists,"wrong":wrongDists}
	pickle.dump(dists,open("{!s}_right_wrong_distances_pickle".format(outputFile),"wb"))
	# ***
	
	# get right/wrong distance distributions for neighbors ***
	neighborDists = {"right":rightNeighborsDists, "wrong":wrongNeighborsDists}
	pickle.dump(neighborDists,open("{!s}_neighbor_right_wrong_distances_pickle".format(outputFile),"wb"))
	# ***

	pickle.dump(totalCluster,open("{!s}_clusters_pickle".format(outputFile),"wb"))
	pickle.dump(allContigs,open("{!s}_contigs_pickle".format(outputFile),"wb"))

	# get distances between extant clusters
	#toMatch = baseName.rsplit("/",1)[0]+"/l1"
	#fSeed = baseName.rsplit("/",1)[0]+"/l2"
	#DB = baseName + "_finalDB"
	#os.system("ls {!s}* > {!s}".format(genePath,toMatch))
	#os.system("bash ./ListScript.sh {!s} > {!s}".format(genePath[:-1],fSeed))
	#if scoreFunction == "raiphy":
	#	scoreRAIphyFinal(DB, fSeed, toMatch, computePath, outputFile)
	#else:
	#	scoringMethodFinal[scoreFunction](DB, fSeed, toMatch, outputFile)
	#fOutD = open("{!s}_distances".format(outputFile),'w')
	#fDists = makeDistanceMatrix("{!s}".format(outputFile+"_dists_sorted"))
	#for row in fDists:
	#	 fOutD.write(",".join(str(r) for r in row)+"\n")
	#fOutD.close()

	# Get rid of files we're not using any more
	os.system("rm -r {!s}".format(genePath))
	os.system("rm {!s}".format(DB))
	os.system("rm {!s}".format(toMatch))
	os.system("rm {!s}".format(fSeed))
	for i in range(leng+1):
		os.system("rm {!s}_{!s}*".format(baseName,i))


if __name__ == "__main__":
	main(sys.argv[1:])