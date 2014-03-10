#! /usr/bin/env python

'''bootstrap.py - wrapper class for my MS project.'''

import sys, getopt, string, os, re, pprint
from bootstrapConstants import *
from bootstrapUtils import *
from bootstrapClasses import *

def main(argv):
	# get inputs, check validity
	inputFile = ''
	outputFile = ''
	coolingSchedule = []
	raiPath = ''
	# Command line arguments
	try:
		opts, args = getopt.getopt(argv,"hi:o:r:c:p:",["ifile=","ofile=","reference=","db=","schedule=","path="])
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
		temp = open(outputFile+"_clusters", 'w')
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
		workingFile = fNext
		thr = int(float(coolingSchedule[i])*1000.0) # so you can have a threshold of 1500
		bgr = "{!s}_{!s}_next".format(baseName,i)
		smlr = "{!s}_{!s}".format(baseName, i)
		#print("Thr = {!s}".format(thr))
		#print("Bgr = {!s}".format(bgr))
		#print("Smlr = {!s}".format(smlr))
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

	pprint.pprint(allContigs)

	ct = 0
	rightDists = []
	wrongDists = []
	thresh = matchThreshold
	close = closeThreshold
	newClusters = []
	
	# Main loop: iterate through cooling schedule, creating databases, making matches, and once matches are made, concatenate each seed (pseudo)contig with matched contigs to make next round
	#for i in range(leng-1,-1,-1):
	for i in [leng-1]:
		# Make DB out of fSeed, whatever it is right now
		print coolingSchedule[i]*1000
		print fSeed
		DB = "{!s}_{!s}_DB".format(baseName,i)
		os.system("{!s}rait -new -i {!s}-2 -o {!s} >/dev/null 2>&1".format(raiPath, fSeed, DB))
		# Process contigs to match
		matches = "{!s}_{!s}_matches".format(baseName,i)
		toMatch = "{!s}_{!s}".format(baseName,i)
		f = open(toMatch+"-2",'r')
		for l in f.readlines():
			sp = l.rstrip().split("\t")
			nm = sp[0]
			fi = sp[1]
			co = Contig(nm,fi)
			allContigs[nm] = co
		# Match ith contigs to DB
		os.system("{!s}rai -I {!s}-1 -d {!s} >/dev/null 2>&1".format(raiPath, toMatch, DB))
		short = toMatch.rsplit("/",1)[1]
		os.system("cp {!s}/{!s}-1.bin {!s}".format(os.getcwd(), short, matches)) # moves results to results folder
		os.system("rm {!s}/{!s}-1.bin".format(os.getcwd(), short))
		
		pprint.pprint(allContigs)
		
		# Construct matching dictionary for internal use
		matchDict = {}
		newClusters = []
		fMatch = open(matches,'r')
		lns = fMatch.readlines()
		dbNames = lns[0].rstrip().split(",")
		contigNames = lns[1].rstrip().split(",")
		for row in range(2,len(lns)):
			line = lns[row]
			bestMatch = line.rstrip().split(", ")[0].split(":")
			bestIndex = int(bestMatch[1])
			bestScore = float(bestMatch[0])
			parent = dbNames[bestIndex]
			child = contigNames[row-2]
			co = allContigs[child]
			if bestScore > thresh: # check if contig is close enough to add
				if parent in matchDict:
					matchDict[parent].append(child)
				else:
					matchDict[parent] = [child]
			else:
				newClusters.append(child)
			for l in line.rstrip().split(", ")[1:]:
				entry = l.split(":")
				score = float(entry[0])
				if score >= (1.0-close)*bestScore:
					ind = int(entry[1])
					name = dbNames[ind]
					co.goodMatches.append([name,score])

			# *** check correctness of match
			cl = allContigs[parent].myCluster
			if checkCorrectSpeciesOlsonFormat(cl.seed, child) == 1:
					rightDists.append(bestDist)
			else:
				wrongDists.append(bestDist)
			# ***

		fMatch.close()
	
		# Prepare next DB
		fSeed = "{!s}_{!s}_seed".format(baseName, i)
		l2 = open(fSeed + "-2",'w')
		# Make concatenated seeds
		for j in matchDict.keys():
			cl = allContigs[j].myCluster
			newContig = "pseudocontig_"+"{!s}".format(ct).zfill(3)
			fpc = open("{!s}{!s}.fna".format(genePath,newContig),'w')
			fpc.write(">{!s}\n".format(newContig))
			_, seq = readSequence("{!s}{!s}.fna".format(genePath, j))
			fpc.write(seq)
			os.system("rm {!s}{!s}.fna".format(genePath,j)) # clear up space
			nCo = Contig(newContig, myCluster = cl)
			cl.root = newContig
			cl.addNode(newContig, j)
			for v in matchDict[j]:
				cl.addNode(newContig, v)
				co = allContigs[v]
				for m in co.goodMatches:
                    cl.addMatch(m)
				_, seq = readSequence("{!s}{!s}.fna".format(genePath, v))
				fpc.write(seq)
				os.system("rm {!s}{!s}.fna".format(genePath,v)) # clear up space
			allContigs[newContig] = nCo
			fpc.write("\n")
			fpc.close()
			l2.write("{!s}\t{!s}{!s}.fna\n".format(newContig,genePath,newContig))
			ct += 1
		for j in newClusters:
			cl = Cluster(j)
			allClusters[j] = cl
			co = allContigs[j]
			co.myCluster = cl
			l2.write("{!s}\t{!s}{!s}.fna\n".format(j,genePath,j))
		l2.close()
        for cl in allContigs:
            print allContigs[cl].seed
            print allContigs[cl].closeList
            print

	# process results from main loop to get clusters and distances
	fOutC = open("{!s}_clusters".format(outputFile),'w')
	for c in allClusters:
		fOutC.write("{!s}\n".format(allClusters[c].get_All()))
	fOutC.close()

	# get distances between extant clusters
	#toMatch = baseName.rsplit("/",1)[0]+"/l1"
	#fSeed = baseName.rsplit("/",1)[0]+"/l2"
	#DB = baseName + "_finalDB"
	#print toMatch
	#print fSeed
	#print DB
	#os.system("ls {!s}* > {!s}".format(genePath,toMatch))
	#os.system("bash ./ListScript.sh {!s} > {!s}".format(genePath[:-1],fSeed))
	#os.system("{!s}rait -new -i {!s} -o {!s} >/dev/null 2>&1".format(raiPath, fSeed, DB))
	#os.system("{!s}rai -I {!s} -d {!s} >/dev/null 2>&1".format(raiPath, toMatch, DB))
	#short = toMatch.rsplit("/",1)[1]
	#os.system("cp {!s}/{!s}.bin {!s}".format(os.getcwd(), short, outputFile+"_dists_sorted")) # moves results to results folder
	#os.system("rm {!s}/{!s}.bin".format(os.getcwd(), short))
	#fOutD = open("{!s}_distances".format(outputFile),'w')
	#fDists = makeDistanceMatrix("{!s}".format(outputFile+"_dists_sorted"))
	#for row in fDists:
	#	 fOutD.write(",".join(str(r) for r in row)+"\n")
	#fOutD.close()

	# get right/wrong distance distributions ***
	#fOutDiff = open("{!s}_right_wrong_distances".format(outputFile),'w')
	#fOutDiff.write("Correct:\n" + ",".join(str(r) for r in rightDists) + "\n")
	#fOutDiff.write("Incorrect:\n" + ",".join(str(w) for w in wrongDists) + "\n")
	#fOutDiff.close()
	# ***

	# Get rid of files we're not using any more
	#os.system("rm -r {!s}".format(genePath))
	#os.system("rm {!s}".format(DB))
	#os.system("rm {!s}".format(toMatch))
	#os.system("rm {!s}".format(fSeed))
	#for i in range(leng+1):
	#	 os.system("rm {!s}_{!s}*".format(baseName,i))


if __name__ == "__main__":
	main(sys.argv[1:])