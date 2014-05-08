#! /usr/bin/env python

'''bootstrap.py - wrapper class for my MS project.'''

import sys, getopt, string, os, re, pprint
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
	joinThreshold = 0.5
	neighborThreshold = 0.1
	# Command line arguments
	try:
		opts, args = getopt.getopt(argv,"hi:o:c:p:s:f:j:n:",["ifile=","ofile=","cut=","path=","score=","namefile=","jointhreshold=","neighborthreshold="])
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
		elif opt in ("-j", "--joiningthreshold"):
			joinThreshold = float(arg)
		elif opt in ("-n", "--neighborthreshold"):
			neighborThreshold = float(arg)
		elif opt in ("-f", "--namefile"):
			nameFile = arg
			nf = open(nameFile,'r')
			names = []
			for l in nf.readlines():
			    n = l.rstrip()
			    if len(n) > 0:
			        names.append(n)
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
	
	rightDistsMax = {"{!s}-{!s}".format(str(coolingSchedule[i]).zfill(2),str(coolingSchedule[i+1]).zfill(2)):[] for i in range(leng-2)}
	wrongDistsMax = {"{!s}-{!s}".format(str(coolingSchedule[i]).zfill(2),str(coolingSchedule[i+1]).zfill(2)):[] for i in range(leng-2)}
	
	rightNeighborsDistsMax = {"{!s}-{!s}".format(str(coolingSchedule[i]).zfill(2),str(coolingSchedule[i+1]).zfill(2)):[] for i in range(leng-2)}
	wrongNeighborsDistsMax = {"{!s}-{!s}".format(str(coolingSchedule[i]).zfill(2),str(coolingSchedule[i+1]).zfill(2)):[] for i in range(leng-2)}
	
	log = open("{!s}_log".format(outputFile),'w')
	neighborLog = open("{!s}_neighborlog".format(outputFile),'w')
	# ***
	
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
		print matches
		for d in dbNames:
			matchDict[d] = []
		for row in range(2,len(lns)):
			line = lns[row]
			bestMatch = line.rstrip().split(", ")[0].split(":")
			if len(bestMatch) < 2:
				print row, line.rstrip()
			bestIndex = int(bestMatch[1])
			bestScore = float(bestMatch[0])
			parent = dbNames[bestIndex]
			child = contigNames[row-2]
			co = allContigs[child]
			matchDict[parent].append(child)
			mythresh = 0.0
			if bestScore > 0:
				mythresh = (1.0-neighborThreshold)*bestScore
			else:
				mythresh = (1.0+neighborThreshold)*bestScore
			for l in line.rstrip().split(", ")[0:]:
				entry = l.split(":")
				score = float(entry[0])
				if score >= mythresh:
					ind = int(entry[1])
					name = dbNames[ind]
					mco = allContigs[name]
					mcl = mco.myCluster
					co.goodMatches.append([mcl.seed, score])
					#print str(mcl)
					#_, mclMax = mcl.purityMax(names)
					#co.goodMatchesMax.append([mclMax, score])
					#print "\n"
					
					# *** check quality of neighbors 
					#corrSeed = correctnessDictSeed[matchLevel](child, mcl.seed)
					#if corrSeed == 1:
					#	rightNeighborsDistsSeed[iterString].append(score)
					#else:
					#	wrongNeighborsDistsSeed[iterString].append(score)
					#if i < leng-1:
					#	corrMax = correctnessDictMax[matchLevel](child, mcl, names)
					#	if corrMax == 1:
					#		rightNeighborsDistsMax[iterString].append(score)
					#	elif corrMax == 0:
					#		wrongNeighborsDistsMax[iterString].append(score)
					# *** 
			# *** check correctness of match
			cl = allContigs[parent].myCluster
			correctS = correctnessDictSeed[matchLevel](child, cl.seed)
			if correctS == 1:
				rightDistsSeed[iterString].append(bestScore)
			else:
				wrongDistsSeed[iterString].append(bestScore)
				log.write("Wrong match (seed): {!s} to {!s}\n".format(child, cl.seed))
				log.write("Original score: {!s}\n".format(bestScore))
				log.write("Matches within {!s}%: {!s}\n\n".format(neighborThreshold*100, co.goodMatches))
			
			#if i < leng-1:	
			#	correctM = correctnessDictMax[matchLevel](child, cl, names)
			#	if correctM == 1:
			#		rightDistsMax[iterString].append(bestScore)
			#	elif correctM == 0:
			#		wrongDistsMax[iterString].append(bestScore)
			#		log.write("Wrong match (max): {!s} to {!s}\n".format(child, cl.seed))
			#		log.write("Original score: {!s}\n".format(bestScore))
			#		log.write("Matches within {!s}%: {!s}\n\n".format(neighborThreshold*100, co.goodMatchesMax))
			# ***
		
		fMatch.close()
		
		# Prepare for next DB creation
		fSeed = "{!s}_{!s}_seed".format(baseName, i)
		l2 = open(fSeed + "-2",'w')
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
				co = allContigs[v]
				#matchingScore = sorted(cl.closeList[co])[-1]
				#if matchingScore > 
				cl.addNode(newContig, v)
				for m in co.goodMatches:
					cl.addMatch(m)
				_, seq = readSequence("{!s}{!s}.fna".format(genePath, v))
				fpc.write(seq)
				#os.system("rm {!s}{!s}.fna".format(genePath,v)) # clear up space
			allContigs[newContig] = nCo
			fpc.write("\n")
			fpc.close()
			ct += 1
			
		toWrite = set(allClusters.keys())
		#print "keys:", toWrite
		toPop = set()
		# merge clusters as appropriate	
		for clust in allClusters:
			cl = allClusters[clust]
			# get info on cluster closeness
			if i < leng-1 and len(cl.closeList) > 0 and cl not in toPop:
				ratioS, bestS = cl.getMostCommonNeighbor()
				if ratioS > joinThreshold and bestS not in toPop:
					myRoot = cl.root
					bestSCL = allClusters[bestS]
					bestRoot = bestSCL.root
					newContig = "pseudocontig_"+"{!s}".format(ct).zfill(3)
					fpc = open("{!s}{!s}.fna".format(genePath,newContig),'w')
					fpc.write(">{!s}\n".format(newContig))
					_, seq1 = readSequence("{!s}{!s}.fna".format(genePath, myRoot))
					fpc.write(seq1)
					_, seq2 = readSequence("{!s}{!s}.fna".format(genePath, bestRoot))
					fpc.write(seq2)
					fpc.write("\n")
					fpc.close()
					nCo = Contig(newContig, myCluster = cl)
					allContigs[newContig] = nCo
					cl.addClusters([bestSCL],newContig)
					#print "Adding {!s} and {!s}".format(cl.seed.split("_")[0], bestSCL.seed.split("_")[0])
					toPop.add(bestS)
					toWrite.remove(bestS)
					ct += 1
				else:
					newContig = cl.root
			
		print "toWrite:", i, ":", toWrite	
		for item in toWrite:
			cl = allClusters[item]
			newContig = cl.root
			l2.write("{!s}\t{!s}{!s}.fna\n".format(newContig,genePath,newContig))
		
		for item in toPop:
			#print "Popped {!s}".format(item)
			allClusters.pop(item)
					
		# *** compute correctness distributions
		rdata = rightDistsSeed[iterString]
		wdata = wrongDistsSeed[iterString]
		comparisonPlot(rdata, wdata, iterString, outputFile, "seed", "Correct distances", "Incorrect Distances")
		
		#if i < leng-1:
		#	rdata = rightNeighborsDistsSeed[iterString]
		#	wdata = wrongNeighborsDistsSeed[iterString]
		#	if rdata:
		#		print "Right range: {:03.2f}%-{:03.2f}%".format(min(rdata)*100.0, max(rdata)*100.0)
		#	if wdata:
		#		print "Wrong range: {:03.2f}%-{:03.2f}%".format(min(wdata)*100.0, max(wdata)*100.0)
		#	print ",".join("{:03.2f}%".format(d*100.0) for d in sorted(rdata))
		#	print ",".join("{:03.2f}%".format(d*100.0) for d in sorted(wdata))
		#	print len(rdata), len(wdata)
		#	comparisonPlot(rdata, wdata, iterString, outputFile, "neighbors_seed", "Distances between correct neighbors", "Distances between incorrect neighbors")

		#if i < leng-1:
		#	rdata = rightDistsMax[iterString]
		#	wdata = wrongDistsMax[iterString]
		#	comparisonPlot(rdata, wdata, iterString, outputFile, "max", "Correct distances", "Incorrect Distances")

		#	rdata = rightNeighborsDistsMax[iterString]
		#	wdata = wrongNeighborsDistsMax[iterString]
		#	comparisonPlot(rdata, wdata, iterString, outputFile, "neighbors_max", "Distances between correct neighbors", "Distances between incorrect neighbors")
		# ***
		
		for n in allClusters:
			cl = allClusters[n]
			cl.closeList = {}
			
		print iterString + " done"
		l2.close()
		
		#if i < leng -1:
		#	neighborLog.write(iterString+"\n")
		#	for a in range(len(mergeLog)):
		#		item = mergeLog[a]
		#		corr = checkCorrectGenusOlsonFormat(item[0],item[2])
		#		if corr == 1:
		#			neighborLog.write("{!s} (max {!s}):\n    {:03.2f}%, {!s} (max {!s})\n".format(item[0],item[1],item[4]*100.0,item[2],item[3]))
		#		else:
		#			neighborLog.write("{!s} (max {!s}):\nXXXX{:03.2f}%, {!s} (max {!s})\n".format(item[0],item[1],item[4]*100.0,item[2],item[3]))
		#	neighborLog.write("\n")

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
	dists={"rightMax":rightDistsMax,"wrongMax":wrongDistsMax, "rightSeed":rightDistsSeed, "wrongSeed":wrongDistsSeed}
	pickle.dump(dists,open("{!s}_right_wrong_distances_pickle".format(outputFile),"wb"))
	# ***
	
	# get right/wrong distance distributions for neighbors ***
	neighborDists = {"rightMax":rightNeighborsDistsMax, "wrongMax":wrongNeighborsDistsMax, "rightSeed":rightNeighborsDistsSeed, "wrongSee":rightNeighborsDistsSeed}
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
	#os.system("rm -r {!s}".format(genePath))
	#os.system("rm {!s}".format(DB))
	#os.system("rm {!s}".format(toMatch))
	#os.system("rm {!s}".format(fSeed))
	#for i in range(leng+1):
	#	os.system("rm {!s}_{!s}*".format(baseName,i))


if __name__ == "__main__":
	main(sys.argv[1:])